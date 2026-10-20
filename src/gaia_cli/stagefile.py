"""``gaia stage`` – stage several datasets described by one ``stage.yaml``.

Every run goes through :func:`execute`: a stage file, a built-in preset filled
in with ``-i``/``-s``/``-e``, or a single dataset's flags (``gaia stage prism``),
which build the same configuration. See dev/stage-file.md for the file format.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import difflib
import inspect
import json
import logging
import re
import sys
import typing
from enum import Enum
from pathlib import Path
from types import ModuleType
from typing import Annotated, Any

import geopandas as gpd
import pandas as pd
import shapely
import xarray as xr
import yaml
from cyclopts import Parameter
from xarray.backends.file_manager import FILE_CACHE

from . import (
    burn_severity,
    dem,
    hrrr,
    io,
    nlcd,
    prism,
    provenance,
    regrid,
    soil,
    synoptic,
)

log = logging.getLogger(__name__)

TOP_KEYS = (
    "name",
    "time",
    "origin",
    "attrs",
    "aoi",
    "aois",
    "window",
    "output",
    "datasets",
)
# Each module has load(), EXTENT (a STAC Collection extent), KIND and TIMED, and
# optionally VARIANT: a load() option named by the key (`dem/usgs10m`), each
# value its own top-level group; and REGRID_METHOD: `gaia align`'s default for
# its raster variables (else nearest)
DATASETS: dict[str, ModuleType] = {
    "synoptic": synoptic,
    "prism": prism,
    "hrrr": hrrr,
    "nlcd": nlcd,
    "dem": dem,
    "soil": soil,
    "burn_severity": burn_severity,
}
# KIND -> group prefix in the output store
GROUPS = {"raster": "rasters", "stations": "stations"}
# Built-in stage files without AOIs or windows, used as --config NAME -i -s -e
PRESETS = Path(__file__).parent / "presets"
WINDOW_ARGS = ("start_date", "end_date")
WINDOW_KEYS = ("start", "end", "before", "after", "select", "count")
NAMED_WINDOWS = ("pre", "post")
SELECT = ("all", "first", "last", "nearest")
DURATION = re.compile(r"(\d+(?:\.\d+)?)(s|min|h|d|w)")
SECONDS = {"s": 1, "min": 60, "h": 3600, "d": 86400, "w": 604800}
DATE_ONLY = re.compile(r"\d{4}-\d{2}-\d{2}")
# `name` names the store and STAC Item in the output folder
NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
# The resolved configuration, written to the output folder next to the store
CONFIG_FILE = "stage.yaml"
AOI_FORMS = ({"bbox"}, {"intersects"}, {"latitude", "longitude", "radius_km"})
# Timed datasets are read this far beyond each end of their window, so that
# periods straddling the window's edges (e.g. PRISM's 12 UTC day) are complete
# when `gaia align` aggregates one product into another's time steps
PAD = pd.Timedelta(days=1)


class StageFileError(ValueError):
    """Every problem found in a stage file, reported together."""

    def __init__(self, path: str, problems: list[str]):
        self.problems = problems
        lines = "\n".join(f"  - {p}" for p in problems)
        super().__init__(f"{path} has {len(problems)} problem(s):\n{lines}")


@dataclasses.dataclass
class Window:
    """A half-open UTC interval ``[start, end)``; ``None`` is an open end."""

    start: pd.Timestamp | None
    end: pd.Timestamp | None
    name: str | None = None
    select: str | None = None
    count: int | None = None

    def __str__(self) -> str:
        return "/".join(
            ".." if t is None else t.isoformat() for t in (self.start, self.end)
        )

    def padded(self) -> Window:
        """The interval actually read: ``PAD`` beyond each closed end."""
        return dataclasses.replace(
            self,
            start=None if self.start is None else self.start - PAD,
            end=None if self.end is None else self.end + PAD,
        )


@dataclasses.dataclass
class Request:
    """One dataset to stage: its AOI, windows and options."""

    name: str  # the datasets key, `<dataset>` or `<dataset>/<variant>`
    aoi: str | dict
    windows: list[Window]
    options: dict[str, Any]

    @property
    def dataset(self) -> str:
        return self.name.partition("/")[0]


@dataclasses.dataclass
class StageConfig:
    name: str
    output: str | None
    time: pd.Timestamp | None
    origin: dict | None
    attrs: dict
    requests: list[Request]


# ---- parsing helpers: each raises ValueError with a message for the user ----


def _unknown(keys, allowed, what: str) -> list[str]:
    out = []
    for key in keys:
        if key not in allowed:
            close = difflib.get_close_matches(str(key), allowed, n=1)
            hint = f"; did you mean '{close[0]}'?" if close else ""
            out.append(f"unknown {what} '{key}'{hint}")
    return out


def parse_time(value, *, end: bool = False):
    """An RFC 3339 timestamp or date as naive UTC.

    A date-only *end* includes that whole day, so it becomes the next midnight.
    """
    date_only = (isinstance(value, dt.date) and not isinstance(value, dt.datetime)) or (
        isinstance(value, str) and DATE_ONLY.fullmatch(value) is not None
    )
    if not isinstance(value, (str, dt.date)):
        raise TypeError(f"{value!r} is not a date or timestamp")
    try:
        t = pd.Timestamp(value)
    except ValueError:
        raise ValueError(f"{value!r} is not a date or timestamp") from None
    if t.tzinfo is not None:
        t = t.tz_convert("UTC").tz_localize(None)
    return t + pd.Timedelta(days=1) if end and date_only else t


def parse_duration(value):
    """A non-negative duration such as ``30d``, ``17min`` or ``PT17M``."""
    s = str(value).strip()
    if m := DURATION.fullmatch(s):
        return pd.Timedelta(seconds=float(m[1]) * SECONDS[m[2]])
    if s.startswith("P"):
        try:
            return pd.Timedelta(s)
        except ValueError:
            pass
    raise ValueError(
        f"{value!r} is not a duration; use e.g. 30d, 17min or 6h (units s, min, h, d, w)"
    )


def parse_window(
    spec, time: pd.Timestamp | None, name: str | None = None
) -> list[Window]:
    """Resolve a window spec (see dev/stage-file.md) into one or more windows."""
    if isinstance(spec, str):
        parts = spec.split("/")
        if len(parts) != 2:
            raise ValueError(f"{spec!r} is not an interval; use 'start/end'")
        start, end = (None if p == ".." else p for p in parts)
        spec = {"start": start, "end": end}
    if not isinstance(spec, dict):
        raise TypeError(f"{spec!r} is not a window")
    if name is None and spec and set(spec) <= set(NAMED_WINDOWS):
        return [w for k, sub in spec.items() for w in parse_window(sub, time, k)]
    if problems := _unknown(spec, WINDOW_KEYS, "window key"):
        raise ValueError("; ".join(problems))
    if "start" in spec and "before" in spec or "end" in spec and "after" in spec:
        raise ValueError("give either start or before, and either end or after")
    before, after = (
        parse_duration(spec[k]) if k in spec else None for k in ("before", "after")
    )
    relative = before is not None or after is not None
    if (relative or spec.get("select") == "nearest") and time is None:
        raise ValueError("relative windows need a top-level 'time'")
    anchor = time if relative else None
    start = (
        parse_time(spec["start"])
        if spec.get("start") is not None
        else (anchor - before if before is not None and anchor is not None else anchor)
    )
    end = (
        parse_time(spec["end"], end=True)
        if spec.get("end") is not None
        else (anchor + after if after is not None and anchor is not None else anchor)
    )
    if start is None and end is None:
        raise ValueError("window has neither a start nor an end")
    if start is not None and end is not None and start >= end:
        raise ValueError(f"window starts at {start} but ends at {end}")
    select, count = spec.get("select"), spec.get("count")
    if select is not None and select not in SELECT:
        raise ValueError(f"select must be one of {', '.join(SELECT)}, not {select!r}")
    if count is not None:
        if select is None:
            raise ValueError("count needs select")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError(f"count must be a positive integer, not {count!r}")
    return [Window(start, end, name, select, count)]


def parse_aoi(spec, base_dir: Path) -> str | dict:
    """Check an AOI path, URL or inline geometry; resolve relative paths."""
    if isinstance(spec, str):
        if "://" in spec or Path(spec).is_absolute():
            return spec
        path = base_dir / spec
        if not path.exists():
            raise ValueError(f"AOI file {path} not found")
        return str(path)
    if not isinstance(spec, dict) or set(spec) not in AOI_FORMS:
        raise ValueError(
            f"{spec!r} is not an AOI; use a path or URL, {{bbox}}, {{intersects}} "
            "or {latitude, longitude, radius_km}"
        )
    if "bbox" in spec:
        b = spec["bbox"]
        if not (isinstance(b, list) and len(b) == 4 and all(_is_number(v) for v in b)):
            raise ValueError("bbox must be [west, south, east, north]")
        _check_lonlat(b[0], b[1])
        _check_lonlat(b[2], b[3])
        if b[0] >= b[2] or b[1] >= b[3]:
            raise ValueError("bbox must be [west, south, east, north] with west < east")
    elif "intersects" in spec:
        try:
            geom = shapely.from_geojson(json.dumps(spec["intersects"]))
        except (shapely.errors.ShapelyError, TypeError) as e:
            raise ValueError(f"intersects is not a GeoJSON geometry ({e})") from None
        if geom.geom_type != "Polygon":
            raise ValueError(f"intersects must be a Polygon, not {geom.geom_type}")
    else:
        _check_lonlat(spec["longitude"], spec["latitude"])
        r = spec["radius_km"]
        if not _is_number(r) or r <= 0:
            raise ValueError("radius_km must be a positive number")
    return spec


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _check_lonlat(lon, lat):
    if not (
        _is_number(lon) and _is_number(lat) and -180 <= lon <= 180 and -90 <= lat <= 90
    ):
        raise ValueError(f"({lon}, {lat}) is not a longitude, latitude pair")


def load_aoi(spec: str | dict) -> gpd.GeoDataFrame:
    """Read a path or URL, or build an inline AOI, in EPSG:4326 for inline forms."""
    if isinstance(spec, str):
        return io.load_aoi(spec)
    if "bbox" in spec:
        geom = shapely.box(*spec["bbox"])
    elif "intersects" in spec:
        geom = shapely.from_geojson(json.dumps(spec["intersects"]))
    else:
        lon, lat = spec["longitude"], spec["latitude"]
        # Buffer in a local azimuthal equidistant projection so the radius is true
        aeqd = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +units=m"
        circle = gpd.GeoSeries([shapely.Point(0, 0)], crs=aeqd).buffer(
            spec["radius_km"] * 1000, quad_segs=32
        )
        geom = circle.to_crs("EPSG:4326").iloc[0]
    return gpd.GeoDataFrame(geometry=[geom], crs="EPSG:4326")


def options(load: typing.Callable) -> dict[str, tuple[Any, Any]]:
    """``{option: (type, default)}`` from a dataset's ``load()`` signature."""
    hints = typing.get_type_hints(load)
    return {
        n: (hints.get(n), p.default)
        for n, p in inspect.signature(load).parameters.items()
        if n not in ("aoi", *WINDOW_ARGS)
    }


def coerce(value, hint):
    """Check a YAML *value* against a ``load()`` type hint and convert it."""
    if isinstance(hint, type) and issubclass(hint, Enum):
        try:
            return hint(value)
        except ValueError:
            choices = ", ".join(str(m.value) for m in hint)
            raise ValueError(f"{value!r} is not one of {choices}") from None
    kinds = [typing.get_origin(h) for h in (hint, *typing.get_args(hint))]
    if list in kinds or tuple in kinds:
        items = [value] if isinstance(value, str) else value
        if not (
            isinstance(items, (list, tuple)) and all(isinstance(v, str) for v in items)
        ):
            raise ValueError(f"{value!r} is not a list of strings")
        # A list of Literal strings, e.g. prism's variables, takes only those
        allowed = [
            a
            for h in (hint, *typing.get_args(hint))
            for e in typing.get_args(h)
            if typing.get_origin(e) is typing.Literal
            for a in typing.get_args(e)
        ]
        for v in items:
            if allowed and v not in allowed:
                choices = ", ".join(dict.fromkeys(allowed))
                raise ValueError(f"{v!r} is not one of {choices}")
        return tuple(items)
    if hint in (int, str) and (not isinstance(value, hint) or isinstance(value, bool)):
        raise ValueError(f"{value!r} is not of type {hint.__name__}")
    return value


# ---- loading a file ----------------------------------------------------------


def presets() -> list[str]:
    return sorted(p.stem for p in PRESETS.glob("*.yaml"))


def read_config(path: str) -> dict:
    """The YAML document in *path*; raise :class:`StageFileError` if unreadable."""
    try:
        doc = yaml.safe_load(Path(path).read_text())
    except (OSError, yaml.YAMLError) as e:
        raise StageFileError(path, [str(e)]) from None
    if not isinstance(doc, dict):
        raise StageFileError(path, ["the file must be a mapping of keys to values"])
    return doc


def load_config(path: str) -> StageConfig:
    """Read and validate a stage file; raise :class:`StageFileError` listing
    every problem found. No network calls are made."""
    return parse_config(read_config(path), path, Path(path).resolve().parent)


def parse_config(doc: dict, path: str, base_dir: Path) -> StageConfig:
    """Validate a stage file's contents; relative AOI paths resolve against
    *base_dir*, and *path* names the source in errors."""
    problems = _unknown(doc, TOP_KEYS, "key")

    def check(where: str, fn, *args):
        try:
            return fn(*args)
        except (ValueError, TypeError) as e:
            problems.append(f"{where}: {e}")

    name = doc.get("name")
    if not isinstance(name, str) or not name:
        problems.append("name: required, e.g. name: skagit-2025-12")
    elif not NAME.fullmatch(name):
        problems.append(
            f"name: {name!r} names the store, so use only letters, digits, '.', '_' and '-'"
        )
    time = (
        check("time", parse_time, doc["time"]) if doc.get("time") is not None else None
    )
    origin = doc.get("origin")
    if origin is not None:
        if not isinstance(origin, dict) or set(origin) != {"latitude", "longitude"}:
            problems.append("origin: must be {latitude, longitude}")
        else:
            check("origin", _check_lonlat, origin["longitude"], origin["latitude"])
    attrs = doc.get("attrs") or {}
    if not isinstance(attrs, dict):
        problems.append("attrs: must be a mapping of attribute names to values")
        attrs = {}
    output = doc.get("output")
    if output is not None and not isinstance(output, str):
        problems.append("output: must be a path or URL")

    aois = doc.get("aois") or {}
    if not isinstance(aois, dict):
        problems.append("aois: must be a mapping of names to AOIs")
        aois = {}
    named = {k: check(f"aois.{k}", parse_aoi, v, base_dir) for k, v in aois.items()}
    default_aoi = doc.get("aoi")
    if isinstance(default_aoi, str) and default_aoi in aois:
        default_aoi = named[default_aoi]
    elif default_aoi is not None:
        default_aoi = check("aoi", parse_aoi, default_aoi, base_dir)
    default_window = doc.get("window")
    if default_window is not None:
        check("window", parse_window, default_window, time)

    datasets = doc.get("datasets")
    if not isinstance(datasets, dict) or not datasets:
        problems.append(
            f"datasets: required, a mapping with any of {', '.join(DATASETS)}"
        )
        datasets = {}
    requests = []
    for ds_name, entry in datasets.items():
        dataset, _, variant = str(ds_name).partition("/")
        if dataset not in DATASETS:
            problems += [
                f"datasets: {p}" for p in _unknown([dataset], list(DATASETS), "dataset")
            ]
            continue
        entry = entry or {}
        if not isinstance(entry, dict):
            problems.append(f"{ds_name}: must be a mapping of options")
            continue
        module = DATASETS[dataset]
        opts = options(module.load)
        # The key names the variant option, which the entry can't repeat
        variant_opt = getattr(module, "VARIANT", None)
        key, variant_value = dataset, None
        if variant_opt is not None:
            hint, variant_value = opts.pop(variant_opt)
            if variant:
                variant_value = check(ds_name, coerce, variant, hint)
            key = f"{dataset}/{plain(variant_value)}"
            if key in [r.name for r in requests]:
                problems.append(f"{ds_name}: {key} is listed twice")
        elif variant:
            problems.append(f"{ds_name}: {dataset} has no variants; use {dataset}")
            continue
        timed = module.TIMED
        keys = [k for k in entry if k != "window"]
        problems += [
            f"{ds_name}: {p}"
            for p in _unknown(keys, ["aoi", "window", *opts], "option")
        ]
        if "window" in entry and not timed:
            problems.append(f"{ds_name}: is static and takes no window")

        aoi = entry.get("aoi")
        if isinstance(aoi, str) and aoi in aois:
            aoi = named[aoi]
        elif aoi is not None:
            aoi = check(f"{ds_name}.aoi", parse_aoi, aoi, base_dir)
        elif default_aoi is not None:
            aoi = default_aoi
        else:
            problems.append(f"{ds_name}: no aoi, and the file has no default aoi")

        windows = []
        if timed:
            spec = entry.get("window", default_window)
            if spec is None:
                problems.append(
                    f"{ds_name}: no window, and the file has no default window"
                )
            else:
                windows = check(f"{ds_name}.window", parse_window, spec, time) or []
            if len(windows) > 1 or any(w.name for w in windows):
                problems.append(f"{ds_name}.window: takes one window, not pre/post")
            if any(w.select for w in windows):
                problems.append(
                    f"{ds_name}.window: select is only for scene-based datasets"
                )
            if any(w.start is None or w.end is None for w in windows):
                problems.append(f"{ds_name}.window: needs both a start and an end")

        values = {variant_opt: variant_value} if variant_opt else {}
        for opt, (hint, default) in opts.items():
            if opt in entry:
                values[opt] = check(f"{ds_name}.{opt}", coerce, entry[opt], hint)
            else:
                values[opt] = default
        if aoi is not None:
            requests.append(Request(key, aoi, windows, values))

    if problems:
        raise StageFileError(path, problems)
    assert isinstance(name, str)  # checked above
    return StageConfig(name, output, time, origin, attrs, requests)


def plain(v):
    """An option value as YAML data."""
    if isinstance(v, Enum):
        return v.value
    if isinstance(v, tuple):
        return list(v)
    return v


def resolved(
    cfg: StageConfig,
    output: str | None,
    aois: dict[str, gpd.GeoDataFrame] | None = None,
    skipped: dict[str, str] | None = None,
) -> dict:
    """The configuration after defaults and AOIs are resolved, as plain YAML data."""
    out: dict[str, Any] = {"name": cfg.name, "output": output}
    if cfg.time is not None:
        out["time"] = cfg.time.isoformat()
    if cfg.origin:
        out["origin"] = cfg.origin
    if cfg.attrs:
        out["attrs"] = cfg.attrs
    out["datasets"] = {}
    for req in cfg.requests:
        entry: dict[str, Any] = {"group": group(req), "aoi": req.aoi}
        if aois is not None:
            entry["aoi_bbox"] = [
                round(float(v), 6) for v in aois[_aoi_key(req.aoi)].total_bounds
            ]
        if req.windows:
            entry["window"] = str(req.windows[0])
            entry["read"] = str(req.windows[0].padded())
        entry.update(_entry_options(req))
        if skipped and req.name in skipped:
            entry["skipped"] = skipped[req.name]
        out["datasets"][req.name] = entry
    return out


def stage_file(
    cfg: StageConfig, output: str, aois: dict[str, gpd.GeoDataFrame]
) -> dict:
    """The resolved stage file written next to the store.

    Every dataset has its own AOI, window and options, so
    ``gaia stage --config`` repeats the run exactly. URLs are kept; local AOI
    files are inlined as ``intersects`` Polygons in EPSG:4326, so the file
    works on other machines.
    """
    out: dict[str, Any] = {"name": cfg.name}
    if cfg.time is not None:
        out["time"] = cfg.time.isoformat()
    if cfg.origin:
        out["origin"] = cfg.origin
    if cfg.attrs:
        out["attrs"] = cfg.attrs
    out["output"] = output
    out["datasets"] = {}
    for req in cfg.requests:
        aoi = req.aoi
        if isinstance(aoi, str) and "://" not in aoi:
            geom = aois[aoi].to_crs("EPSG:4326").geometry.iloc[0]
            aoi = {"intersects": json.loads(shapely.to_geojson(geom))}
        entry: dict[str, Any] = {"aoi": aoi}
        if req.windows:
            entry["window"] = str(req.windows[0])
        entry.update(_entry_options(req))
        out["datasets"][req.name] = entry
    return out


def _entry_options(req: Request) -> dict:
    """*req*'s options as YAML data, without the variant its key names."""
    variant_opt = getattr(DATASETS[req.dataset], "VARIANT", None)
    return {k: plain(v) for k, v in req.options.items() if k != variant_opt}


def write_stage_file(doc: dict, output: str, merge: bool = False):
    """Write *doc* as the output folder's ``stage.yaml``; with *merge*, keep
    the datasets that the folder's existing file lists and *doc* doesn't."""
    path = f"{output.rstrip('/')}/{CONFIG_FILE}"
    if merge and (text := io.read_text(path)) is not None:
        old = yaml.safe_load(text) or {}
        doc = {**doc, "datasets": {**old.get("datasets", {}), **doc["datasets"]}}
    io.write_text(path, yaml.safe_dump(doc, sort_keys=False, default_flow_style=None))
    log.info("Wrote %s", path)


def group(req: Request) -> str:
    """The store group: the key for variant datasets (``dem/cop30``), else
    ``<GROUPS[KIND]>/<name>``."""
    module = DATASETS[req.dataset]
    if hasattr(module, "VARIANT"):
        return req.name
    return f"{GROUPS[module.KIND]}/{req.name}"


def _aoi_key(spec: str | dict) -> str:
    return spec if isinstance(spec, str) else yaml.safe_dump(spec, sort_keys=True)


def load_aois(cfg: StageConfig, path: str) -> dict[str, gpd.GeoDataFrame]:
    """Read every distinct AOI once; raise :class:`StageFileError` on failures."""
    aois, problems = {}, []
    for req in cfg.requests:
        key = _aoi_key(req.aoi)
        if key not in aois:
            try:
                aois[key] = load_aoi(req.aoi)
            except Exception as e:  # noqa: BLE001  any read failure is a problem to report
                problems.append(f"{req.name}.aoi: could not read {req.aoi!r} ({e})")
                aois[key] = None
    if problems:
        raise StageFileError(path, problems)
    return aois


# ---- coverage ------------------------------------------------------------------


def _split_antimeridian(b: list) -> list[list]:
    """*b* as one or two bboxes with west <= east."""
    if b[0] <= b[2]:
        return [b]
    return [[b[0], b[1], 180.0, b[3]], [-180.0, b[1], b[2], b[3]]]


def _intersects(a, b) -> bool:
    return a[0] <= b[2] and a[2] >= b[0] and a[1] <= b[3] and a[3] >= b[1]


def _utc(value: str | None) -> pd.Timestamp | None:
    return None if value is None else parse_time(value)


def outside_extent(req: Request, aoi: gpd.GeoDataFrame) -> str | None:
    """Why *req* misses its dataset's ``EXTENT``, or None if it may have data.

    A coarse check: the AOI's bounds against the extent's bboxes (the
    sub-regions after the first, when there are any), and each window
    against the temporal intervals.
    """
    extent = DATASETS[req.dataset].EXTENT
    bboxes = extent["spatial"]["bbox"]
    aoi_bbox = [round(float(v), 4) for v in aoi.to_crs("EPSG:4326").total_bounds]
    regions = [p for b in bboxes[1:] or bboxes[:1] for p in _split_antimeridian(b)]
    if not any(_intersects(aoi_bbox, r) for r in regions):
        return f"AOI {aoi_bbox} is outside the dataset's extent {bboxes[0]}"
    intervals = [(_utc(s), _utc(e)) for s, e in extent["temporal"]["interval"]]
    for w in req.windows:
        if not any(
            (w.start is None or e is None or w.start <= e)
            and (w.end is None or s is None or w.end > s)
            for s, e in intervals
        ):
            spans = ", ".join(
                f"{'..' if s is None else s.isoformat()}/{'..' if e is None else e.isoformat()}"
                for s, e in intervals
            )
            return f"window {w} is outside the dataset's time extent {spans}"
    return None


# ---- running -----------------------------------------------------------------


def run(
    cfg: StageConfig,
    aois: dict[str, gpd.GeoDataFrame],
    store: str,
    overwrite: bool = False,
    append: bool = False,
    assets: dict | None = None,
) -> dict[str, str]:
    """Load every requested dataset and write each as a group of one store.

    Each group is written, and the store consolidated and its STAC Item
    rewritten, as soon as it loads. A dataset that fails doesn't stop the
    others; returns ``{name: error}`` for the ones that failed. With
    *append*, the groups are added to the existing store.
    """
    failures = {}
    for req in cfg.requests:
        module = DATASETS[req.dataset]
        kwargs = dict(req.options)
        if req.windows:
            w = req.windows[0].padded()
            assert w.start is not None and w.end is not None  # checked in load_config
            # load() takes inclusive bounds; the window's end is exclusive
            kwargs["start_date"] = w.start.isoformat()
            kwargs["end_date"] = (w.end - pd.Timedelta(seconds=1)).isoformat()
        log.info(
            "%s: loading (%s)",
            req.name,
            f"{req.windows[0]}, read as {req.windows[0].padded()}"
            if req.windows
            else "static",
        )
        try:
            obj = module.load(aois[_aoi_key(req.aoi)], **kwargs)
            # Station data is xvec CF-encoded for the store
            node = module.to_dataset(obj) if module.KIND == "stations" else obj
            if method := getattr(module, "REGRID_METHOD", None):
                for var in node.data_vars.values():
                    if {"x", "y"} <= set(var.dims):
                        var.attrs.setdefault(regrid.METHOD_ATTR, method)
            tree = xr.DataTree.from_dict({group(req): node})
            tree.attrs = dict(cfg.attrs)
            io.save_datatree(
                tree, store, overwrite=overwrite, append=append, assets=assets
            )
            append = True  # later groups go into the same store
        except Exception as e:
            log.error("%s: failed: %s", req.name, e)
            log.debug("%s: traceback", req.name, exc_info=True)
            failures[req.name] = f"{type(e).__name__}: {e}"
    # Close cached remote GDAL handles before interpreter shutdown, where their
    # cleanup can crash sys.excepthook on Python 3.14
    FILE_CACHE.clear()
    return failures


def execute(
    cfg: StageConfig,
    path: str,
    output: str | None,
    only: str | None = None,
    dry_run: bool = False,
    overwrite: bool = False,
    skip_existing: bool = False,
):
    """Check *cfg* (from *path*) and stage it: the one code path for every run.

    Exits with an error listing every problem, or every failed dataset.
    """
    try:
        problems = []
        if only is not None:
            names = [n.strip() for n in only.split(",") if n.strip()]
            # `--only dem` selects every dem/<variant>
            listed = [r.name for r in cfg.requests] + [r.dataset for r in cfg.requests]
            problems += [
                f"--only: {p}"
                for p in _unknown(
                    names, list(dict.fromkeys(listed)), "dataset in the file"
                )
            ]
            cfg.requests = [
                r for r in cfg.requests if r.name in names or r.dataset in names
            ]
        if output is None and not dry_run:
            problems.append("output: required, here or with -o")
        elif output is not None and output.rstrip("/").endswith(".zarr"):
            problems.append(
                f"output: -o is a folder for <name>.zarr, {CONFIG_FILE} and the "
                f"STAC Item; drop the .zarr from {output}"
            )
        if overwrite and skip_existing:
            problems.append("--overwrite and --skip-existing can't be combined")
        if problems:
            raise StageFileError(path, problems)
        aois = load_aois(cfg, path)
    except StageFileError as e:
        sys.exit(f"Error: {e}")

    skipped = {
        req.name: reason
        for req in cfg.requests
        if (reason := outside_extent(req, aois[_aoi_key(req.aoi)]))
    }
    if dry_run:
        print(
            yaml.safe_dump(
                resolved(cfg, output, aois, skipped),
                sort_keys=False,
                default_flow_style=None,
            ),
            end="",
        )
        return
    assert output is not None
    for name, reason in skipped.items():
        log.warning("%s: skipped, %s", name, reason)
    requests = [r for r in cfg.requests if r.name not in skipped]
    if not requests:
        reasons = "; ".join(f"{n}: {r}" for n, r in skipped.items())
        sys.exit(f"Error: nothing to stage ({reasons})")

    store = io.store_path(output, cfg.name)
    append = skip_existing and io.output_exists(store)
    if skip_existing:
        io.check_s3_env(output)
    else:
        io.check_output(output, overwrite)
    write_stage_file(
        stage_file(dataclasses.replace(cfg, requests=requests), output, aois),
        output,
        merge=append,
    )
    if append:
        existing = io.groups(store)
        for req in requests:
            if group(req) in existing:
                log.info("%s: skipped, %s already exists", req.name, group(req))
        requests = [r for r in requests if group(r) not in existing]
        if not requests:
            log.info("Done: every group already exists")
            return
    failures = run(
        dataclasses.replace(cfg, requests=requests),
        aois,
        store,
        overwrite=overwrite,
        append=append,
        assets={"config": provenance.config_asset(f"./{CONFIG_FILE}")},
    )
    if failures:
        lines = "\n".join(f"  - {n}: {e}" for n, e in failures.items())
        written = len(requests) - len(failures)
        hint = " Rerun with --skip-existing to retry them." if written else ""
        sys.exit(
            f"Error: {len(failures)} of {len(requests)} dataset(s) failed; "
            f"{written} written to {store}.{hint}\n{lines}"
        )
    log.info("Done")


def stage(
    # Collects any other tokens, so data flags get a clear error, not "Unknown option"
    *flags: Annotated[str, Parameter(allow_leading_hyphen=True, show=False)],
    config: Annotated[str, Parameter(name=["--config", "-c"])],
    vectorPath: Annotated[str | None, Parameter(name=["--input", "-i"])] = None,
    start_date: Annotated[str | None, Parameter(name=["--start", "-s"])] = None,
    end_date: Annotated[str | None, Parameter(name=["--end", "-e"])] = None,
    output_path: Annotated[str | None, Parameter(name=["--output", "-o"])] = None,
    only: Annotated[str | None, Parameter(name=["--only"])] = None,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
):
    """Stage every dataset in a stage file, or a built-in preset, into one store.

    The file gives the AOIs, time windows and options for each dataset (see
    dev/stage-file.md). It replaces ``-i``, ``-s``, ``-e`` and dataset
    options, which can't be combined with it. A preset (``conus-precip``) is
    a built-in file without an AOI or window: give those with ``-i``, ``-s``
    and ``-e``.

    Datasets whose extent the AOI or window misses are skipped with a
    warning. A dataset that fails doesn't stop the others; every failure is
    reported at the end.

    Parameters
    ----------
    config : str
        Path to a ``stage.yaml`` file, or the name of a preset.
    vectorPath : str
        With a preset only: path to a vector file defining the AOI.
    start_date : str
        With a preset only: start (inclusive), a date or UTC time.
    end_date : str
        With a preset only: end; a date includes the whole day.
    output_path : str
        Output folder for ``<name>.zarr``, its STAC Item ``<name>.json`` and
        the resolved ``stage.yaml``. Default: the file's ``output``.
    only : str
        Comma-separated datasets from the file to stage, e.g. ``prism,hrrr``.
    dry_run : bool
        Validate the file, read the AOIs and print the resolved configuration,
        including datasets skipped for coverage, without fetching any data.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.
    skip_existing : bool
        Add to an existing store, skipping datasets whose group it already
        has. Use it to retry the datasets that failed. Default False.

    Examples
    --------
    Check, then stage, the Skagit Basin December 2025 atmospheric river:
        $ gaia stage --config examples/skagit-2025.yaml --dry-run
        $ gaia stage --config examples/skagit-2025.yaml -o skagit

    Station, PRISM and HRRR precipitation for any CONUS AOI:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage --config conus-precip -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit
    """
    preset = PRESETS / f"{config}.yaml"
    is_preset = config in presets() and not Path(config).exists()
    data_flags = [
        f
        for f, v in (("-i", vectorPath), ("-s", start_date), ("-e", end_date))
        if v is not None
    ]
    if flags or (data_flags and not is_preset):
        sys.exit(
            f"Error: --config can't be combined with {' '.join([*data_flags, *flags])}; "
            f"set the AOI, window and options in {config} instead"
        )
    try:
        if is_preset:
            if vectorPath is None:
                raise StageFileError(
                    config, ["the preset has no AOI; give one with -i (and -s, -e)"]
                )
            doc = read_config(str(preset))
            doc["aoi"] = vectorPath
            # Like a flag-only run, the store is named after the output folder
            if output_path is not None:
                doc["name"] = io.folder_name(output_path)
            window = {"start": start_date, "end": end_date}
            if any(window.values()):
                doc["window"] = {k: v for k, v in window.items() if v is not None}
            cfg = parse_config(doc, config, Path.cwd())
        elif not Path(config).exists():
            raise StageFileError(
                config,
                [f"file not found, and not a preset ({', '.join(presets())})"],
            )
        else:
            cfg = load_config(config)
    except StageFileError as e:
        sys.exit(f"Error: {e}")
    execute(
        cfg,
        config,
        output_path or cfg.output,
        only=only,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
    )


def stage_dataset(
    name: str,
    vectorPath: str,
    start_date: str | None,
    end_date: str | None,
    output_path: str,
    dry_run: bool = False,
    overwrite: bool = False,
    skip_existing: bool = False,
    **options,
):
    """Stage one dataset from its CLI flags (``gaia stage prism -i … -s …``).

    Builds the stage file the flags describe and runs it like any other. The
    store is named after the output folder.
    """
    key = name
    if variant_opt := getattr(DATASETS[name], "VARIANT", None):
        # `gaia stage dem -p usgs10m` is the key dem/usgs10m
        key = f"{name}/{plain(options.pop(variant_opt))}"
    doc: dict[str, Any] = {
        "name": io.folder_name(output_path),
        "aoi": vectorPath,
        "datasets": {key: {k: plain(v) for k, v in options.items()}},
    }
    if DATASETS[name].TIMED:
        doc["window"] = {"start": start_date, "end": end_date}
    try:
        cfg = parse_config(doc, f"gaia stage {name}", Path.cwd())
    except StageFileError as e:
        sys.exit(f"Error: {e}")
    execute(
        cfg,
        f"gaia stage {name}",
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
    )
