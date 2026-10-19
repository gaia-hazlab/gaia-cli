"""``gaia stage --config`` – stage several datasets described by one ``stage.yaml``.

See docs/stage-file.md for the file format.
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
from typing import Annotated, Any

import geopandas as gpd
import pandas as pd
import shapely
import xarray as xr
import yaml
from cyclopts import Parameter
from xarray.backends.file_manager import FILE_CACHE

from . import hrrr, io, nlcd, prism, synoptic

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
# name -> (load function, group prefix in the output store)
DATASETS: dict[str, tuple[typing.Callable, str]] = {
    "synoptic": (synoptic.load, "stations"),
    "prism": (prism.load, "rasters"),
    "hrrr": (hrrr.load, "rasters"),
    "nlcd": (nlcd.load, "rasters"),
}
WINDOW_ARGS = ("start_date", "end_date")
WINDOW_KEYS = ("start", "end", "before", "after", "select", "count")
NAMED_WINDOWS = ("pre", "post")
SELECT = ("all", "first", "last", "nearest")
DURATION = re.compile(r"(\d+(?:\.\d+)?)(s|min|h|d|w)")
SECONDS = {"s": 1, "min": 60, "h": 3600, "d": 86400, "w": 604800}
DATE_ONLY = re.compile(r"\d{4}-\d{2}-\d{2}")
AOI_FORMS = ({"bbox"}, {"intersects"}, {"latitude", "longitude", "radius_km"})


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


@dataclasses.dataclass
class Request:
    """One dataset to stage: its AOI, windows and options."""

    name: str
    aoi: str | dict
    windows: list[Window]
    options: dict[str, Any]


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
    """Resolve a window spec (see docs/stage-file.md) into one or more windows."""
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
            shapely.from_geojson(json.dumps(spec["intersects"]))
        except (shapely.errors.ShapelyError, TypeError) as e:
            raise ValueError(f"intersects is not a GeoJSON geometry ({e})") from None
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
        if not (isinstance(items, list) and all(isinstance(v, str) for v in items)):
            raise ValueError(f"{value!r} is not a list of strings")
        return tuple(items)
    if hint in (int, float, str, bool):
        ok = isinstance(value, hint) and not (
            hint is not bool and isinstance(value, bool)
        )
        if hint is float and _is_number(value):
            ok, value = True, float(value)
        if not ok:
            raise ValueError(f"{value!r} is not of type {hint.__name__}")
    return value


# ---- loading a file ----------------------------------------------------------


def load_config(path: str) -> StageConfig:
    """Read and validate a stage file; raise :class:`StageFileError` listing
    every problem found. No network calls are made."""
    try:
        doc = yaml.safe_load(Path(path).read_text())
    except (OSError, yaml.YAMLError) as e:
        raise StageFileError(path, [str(e)]) from None
    if not isinstance(doc, dict):
        raise StageFileError(path, ["the file must be a mapping of keys to values"])
    base_dir = Path(path).resolve().parent
    problems = _unknown(doc, TOP_KEYS, "key")

    def check(where: str, fn, *args):
        try:
            return fn(*args)
        except (ValueError, TypeError) as e:
            problems.append(f"{where}: {e}")

    name = doc.get("name")
    if not isinstance(name, str) or not name:
        problems.append("name: required, e.g. name: skagit-2025-12")
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
        if ds_name not in DATASETS:
            problems += [
                f"datasets: {p}" for p in _unknown([ds_name], list(DATASETS), "dataset")
            ]
            continue
        entry = entry or {}
        if not isinstance(entry, dict):
            problems.append(f"{ds_name}: must be a mapping of options")
            continue
        load, _ = DATASETS[ds_name]
        opts = options(load)
        timed = all(a in inspect.signature(load).parameters for a in WINDOW_ARGS)
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

        values = {}
        for opt, (hint, default) in opts.items():
            if opt in entry:
                values[opt] = check(f"{ds_name}.{opt}", coerce, entry[opt], hint)
            else:
                values[opt] = default
        if aoi is not None:
            requests.append(Request(ds_name, aoi, windows, values))

    if problems:
        raise StageFileError(path, problems)
    return StageConfig(name, output, time, origin, attrs, requests)


def resolved(
    cfg: StageConfig,
    output: str | None,
    aois: dict[str, gpd.GeoDataFrame] | None = None,
) -> dict:
    """The configuration after defaults and AOIs are resolved, as plain YAML data."""

    def plain(v):
        if isinstance(v, Enum):
            return v.value
        if isinstance(v, tuple):
            return list(v)
        return v

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
        entry.update({k: plain(v) for k, v in req.options.items()})
        out["datasets"][req.name] = entry
    return out


def group(req: Request) -> str:
    return f"{DATASETS[req.name][1]}/{req.name}"


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


# ---- running -----------------------------------------------------------------


def _as_dataset(req: Request, obj) -> xr.Dataset:
    if req.name == "synoptic":
        return synoptic.to_dataset(obj)
    if req.name == "nlcd":
        ds = obj.to_dataset()
        ds.attrs = nlcd.source_attrs(req.options["year"])
        return ds
    return obj


def run(
    cfg: StageConfig,
    aois: dict[str, gpd.GeoDataFrame],
    output: str,
    overwrite: bool = False,
):
    """Load every requested dataset and write them as groups of one store."""
    nodes = {}
    for req in cfg.requests:
        load, _ = DATASETS[req.name]
        kwargs = dict(req.options)
        if req.windows:
            w = req.windows[0]
            assert w.start is not None and w.end is not None  # checked in load_config
            # load() takes inclusive bounds; the window's end is exclusive
            kwargs["start_date"] = w.start.isoformat()
            kwargs["end_date"] = (w.end - pd.Timedelta(seconds=1)).isoformat()
        log.info(
            "%s: loading (%s)", req.name, req.windows[0] if req.windows else "static"
        )
        nodes[group(req)] = _as_dataset(req, load(aois[_aoi_key(req.aoi)], **kwargs))

    tree = xr.DataTree.from_dict(nodes)
    tree.attrs = dict(cfg.attrs)
    io.save_datatree(tree, output, overwrite=overwrite)
    # Close cached remote GDAL handles before interpreter shutdown (see nlcd.stage)
    FILE_CACHE.clear()


def stage(
    # Collects any other tokens, so data flags get a clear error, not "Unknown option"
    *flags: Annotated[str, Parameter(allow_leading_hyphen=True, show=False)],
    config: Annotated[str, Parameter(name=["--config", "-c"])],
    output_path: Annotated[str | None, Parameter(name=["--output", "-o"])] = None,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
):
    """Stage every dataset in a stage file into one store.

    The file gives the AOIs, time windows and options for each dataset (see
    docs/stage-file.md). It replaces ``-i``, ``-s``, ``-e`` and dataset
    options, which can't be combined with it.

    Parameters
    ----------
    config : str
        Path to a ``stage.yaml`` file.
    output_path : str
        Path to write the Zarr store. Default: the file's ``output``.
    dry_run : bool
        Validate the file, read the AOIs and print the resolved configuration
        without fetching any data.
    overwrite : bool
        Replace the output store if it already exists. Default False.

    Examples
    --------
    Check, then stage, the Skagit Basin December 2025 atmospheric river:
        $ gaia stage --config examples/skagit-2025.yaml --dry-run
        $ gaia stage --config examples/skagit-2025.yaml -o skagit.zarr
    """
    if flags:
        sys.exit(
            f"Error: --config can't be combined with {' '.join(flags)}; "
            f"set the AOI, window and options in {config} instead"
        )
    try:
        cfg = load_config(config)
        output = output_path or cfg.output
        if output is None and not dry_run:
            raise StageFileError(config, ["output: required, here or with -o"])
        aois = load_aois(cfg, config)
    except StageFileError as e:
        sys.exit(f"Error: {e}")

    if dry_run:
        print(
            yaml.safe_dump(
                resolved(cfg, output, aois), sort_keys=False, default_flow_style=None
            ),
            end="",
        )
        return
    assert output is not None
    io.check_output(output, overwrite)
    run(cfg, aois, output, overwrite=overwrite)
    log.info("Done")
