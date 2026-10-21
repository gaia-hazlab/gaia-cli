"""``stage.yaml``: read, validate and resolve a stage file.

A stage file gives the AOIs, time windows and options for each dataset; see
dev/stage-file.md. :mod:`gaia_cli.runner` runs it.
"""

from __future__ import annotations

import dataclasses
import inspect
import json
import logging
import re
import typing
from enum import Enum
from pathlib import Path
from types import ModuleType
from typing import Any

import geopandas as gpd
import pandas as pd
import requests
import shapely
import yaml

from . import io
from .aoi import aoi_key, check_lonlat, load_aoi, parse_aoi
from .datasets import DATASETS
from .window import Window, parse_time, parse_window, unknown

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
# KIND -> group prefix in the output store
GROUPS = {"raster": "rasters", "stations": "stations"}
# Built-in stage files without AOIs or windows, used as --config NAME -i -s -e
PRESETS = Path(__file__).parent / "presets"
WINDOW_ARGS = ("start_date", "end_date")
# `name` names the store and STAC Item in the output folder
NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
# The resolved configuration, written to the output folder next to the store
CONFIG_FILE = "stage.yaml"


class StageFileError(ValueError):
    """Every problem found in a stage file, reported together."""

    def __init__(self, path: str, problems: list[str]):
        self.problems = problems
        lines = "\n".join(f"  - {p}" for p in problems)
        super().__init__(f"{path} has {len(problems)} problem(s):\n{lines}")


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


# ---- options -----------------------------------------------------------------


def timed(module: ModuleType) -> bool:
    """Whether the dataset takes a window: its ``load()`` takes start and end dates."""
    return all(a in inspect.signature(module.load).parameters for a in WINDOW_ARGS)


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


def is_url(path: str) -> bool:
    """Whether *path* is a stage file to download rather than a local file."""
    return path.startswith(("http://", "https://"))


def read_config(path: str) -> dict:
    """The YAML document at *path* (a local file or an http(s) URL); raise
    :class:`StageFileError` if unreadable."""
    try:
        if is_url(path):
            r = requests.get(path, timeout=60)
            r.raise_for_status()
            text = r.text
        else:
            text = Path(path).read_text()
        doc = yaml.safe_load(text)
    except (OSError, requests.RequestException, yaml.YAMLError) as e:
        raise StageFileError(path, [str(e)]) from None
    if not isinstance(doc, dict):
        raise StageFileError(path, ["the file must be a mapping of keys to values"])
    return doc


def load_config(path: str) -> StageConfig:
    """Read and validate a stage file; raise :class:`StageFileError` listing
    every problem found. No network calls are made, except to download a
    stage file given as a URL, whose AOIs must then be URLs or inline."""
    base_dir = None if is_url(path) else Path(path).resolve().parent
    return parse_config(read_config(path), path, base_dir)


def parse_config(doc: dict, path: str, base_dir: Path | None) -> StageConfig:
    """Validate a stage file's contents; relative AOI paths resolve against
    *base_dir* (``None`` allows no AOI paths), and *path* names the source in
    errors."""
    problems = unknown(doc, TOP_KEYS, "key")

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
            check("origin", check_lonlat, origin["longitude"], origin["latitude"])
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
                f"datasets: {p}" for p in unknown([dataset], list(DATASETS), "dataset")
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
        is_timed = timed(module)
        keys = [k for k in entry if k != "window"]
        problems += [
            f"{ds_name}: {p}" for p in unknown(keys, ["aoi", "window", *opts], "option")
        ]
        if "window" in entry and not is_timed:
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
        if is_timed:
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
    aois: dict[str, gpd.GeoDataFrame],
    skipped: dict[str, str] | None = None,
) -> dict:
    """The configuration after defaults and AOIs are resolved, as YAML data.

    This is the ``stage.yaml`` written next to the store: every dataset has
    its own AOI, window and options, so ``gaia stage --config`` repeats the
    run exactly. URLs are kept; local AOI files are inlined as ``intersects``
    Polygons in EPSG:4326, so the file works on other machines.

    With *skipped* (``--dry-run``), AOIs stay as given, and each dataset also
    lists its group, its AOI's bounds, the interval read and why it is skipped.
    """
    dry_run = skipped is not None
    out: dict[str, Any] = {"name": cfg.name, "output": output}
    if cfg.time is not None:
        out["time"] = cfg.time.isoformat()
    if cfg.origin:
        out["origin"] = cfg.origin
    if cfg.attrs:
        out["attrs"] = cfg.attrs
    out["datasets"] = {}
    for req in cfg.requests:
        entry: dict[str, Any] = {}
        if dry_run:
            bounds = aois[aoi_key(req.aoi)].total_bounds
            entry["group"] = group(req)
            entry["aoi"] = req.aoi
            entry["aoi_bbox"] = [round(float(v), 6) for v in bounds]
        elif isinstance(req.aoi, str) and "://" not in req.aoi:
            geom = aois[req.aoi].to_crs("EPSG:4326").geometry.iloc[0]
            entry["aoi"] = {"intersects": json.loads(shapely.to_geojson(geom))}
        else:
            entry["aoi"] = req.aoi
        if req.windows:
            entry["window"] = str(req.windows[0])
            if dry_run:
                entry["read"] = str(req.windows[0].padded())
        entry.update(_entry_options(req))
        if dry_run and req.name in skipped:
            entry["skipped"] = skipped[req.name]
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


def load_aois(cfg: StageConfig, path: str) -> dict[str, gpd.GeoDataFrame]:
    """Read every distinct AOI once; raise :class:`StageFileError` on failures."""
    aois, problems = {}, []
    for req in cfg.requests:
        key = aoi_key(req.aoi)
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
