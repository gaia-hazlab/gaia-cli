"""``gaia check`` – summarize a store's groups and report what ``gaia align`` can't use.

Works on any Zarr v3 store, not only one ``gaia stage`` wrote: from a fork, another
tool, or by hand. Each group is recognized from its structure, not its path: a
raster has ``x`` and ``y`` dims, a stations group a CF ragged array (``row_size``
along ``obs``), a waveforms group ``sites`` and ``channels`` children.

Errors make ``gaia align`` skip a group or give wrong values (a period without
``cell_methods`` is averaged, not summed). Warnings are missing metadata (units,
the source's title, license, providers, links). See dev/adding-a-dataset.md for
the conventions.
"""

from __future__ import annotations

import logging
import re
import sys
from typing import Annotated, NamedTuple, get_args

import numpy as np
import pandas as pd
import pyproj
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from cyclopts import Parameter

from . import io, provenance, regrid

log = logging.getLogger(__name__)

ERROR, WARNING = "error", "warning"
# Group attributes describing the source (dev/adding-a-dataset.md)
SOURCE_ATTRS = ("title", "license", "providers", "links")
# cell_methods time reductions gaia align recognizes; any other is averaged
TIME_METHODS = ("sum", "mean", "point", "minimum", "maximum")
WAVEFORM_TABLES = ("sites", "channels")


class Finding(NamedTuple):
    level: str
    group: str
    message: str


def kind(node: xr.DataTree) -> str | None:
    """``raster``, ``stations`` or ``waveforms``, from the group's structure."""
    ds = node.dataset
    if {"x", "y"} <= set(ds.dims):
        return "raster"
    if "row_size" in ds and "obs" in ds.dims:
        return "stations"
    if all(t in node.children for t in WAVEFORM_TABLES):
        return "waveforms"
    return None


def _bounds_name(ds: xr.Dataset) -> str | None:
    if "time" not in ds.coords:
        return None
    return ds.time.attrs.get("bounds") or ds.time.encoding.get("bounds")


def _duration(td) -> str:
    """``1D``, ``6h``, ``-12h``, ``20ms``: the largest whole unit."""
    s = float(np.timedelta64(td, "ns").astype("int64")) / 1e9
    if s == 0:
        return "0"
    sign, s = ("-", -s) if s < 0 else ("", s)
    for size, label in ((86400, "D"), (3600, "h"), (60, "min"), (1, "s")):
        if s >= size and s % size == 0:
            return f"{sign}{s // size:g}{label}"
    return f"{sign}{s * 1000:g}ms"


def _time(t) -> str:
    return f"{pd.Timestamp(t):%Y-%m-%dT%H:%M:%S}".removesuffix(":00").removesuffix(
        "T00:00"
    )


def time_findings(path: str, ds: xr.Dataset) -> list[Finding]:
    """The time coordinate is datetime64 and increasing, and its CF bounds, if
    it names any, are a ``(time, 2)`` variable of periods."""
    if "time" not in ds.coords:
        return []
    time = ds.time
    if not np.issubdtype(time.dtype, np.datetime64):
        return [Finding(ERROR, path, f"time is {time.dtype}, not datetime64")]
    out = []
    if time.dims == ("time",) and not time.to_index().is_monotonic_increasing:
        out.append(Finding(ERROR, path, "time is not increasing"))
    name = _bounds_name(ds)
    if name is None:
        return out
    if name not in ds.variables:
        msg = f"time names bounds {name!r}, which is missing: gaia align treats the values as instantaneous"
        return [*out, Finding(ERROR, path, msg)]
    b = ds[name]
    if b.ndim != 2 or b.dims[0] != "time" or b.shape[1] != 2:
        msg = f"{name} has dims {b.dims} and shape {b.shape}, not (time, 2)"
        out.append(Finding(ERROR, path, msg))
    elif (b.isel({b.dims[1]: 0}) > b.isel({b.dims[1]: 1})).any():
        out.append(Finding(ERROR, path, f"{name} has periods ending before they start"))
    return out


def _cell_methods_findings(
    path: str, name: str, da: xr.DataArray, has_bounds: bool
) -> list[Finding]:
    cm = da.attrs.get("cell_methods")
    if cm is None:
        if has_bounds:
            msg = f"{name} has time bounds but no cell_methods: gaia align averages it (time: mean)"
            return [Finding(ERROR, path, msg)]
        return []
    m = re.search(r"\btime:\s*(\w+)", cm)
    if m is None or m.group(1) not in TIME_METHODS:
        msg = f"{name}: cell_methods {cm!r} is not time: {', '.join(TIME_METHODS)}; gaia align averages it"
        return [Finding(WARNING, path, msg)]
    return []


def raster_findings(path: str, ds: xr.Dataset) -> list[Finding]:
    """A CRS and regular grid, and per variable: units, a known regrid method,
    and ``cell_methods`` for values that cover periods."""
    out = []
    if "spatial_ref" not in ds.coords:
        out.append(
            Finding(
                ERROR, path, "no spatial_ref coordinate: gaia align skips this group"
            )
        )
    elif ds.rio.crs is None:
        out.append(Finding(ERROR, path, "spatial_ref holds no CRS"))
    for dim in ("x", "y"):
        step = np.diff(ds[dim].values)
        if step.size and not np.allclose(step, step[0]):
            msg = f"{dim} spacing is irregular; gaia align expects a regular grid"
            out.append(Finding(ERROR, path, msg))
    bounds = _bounds_name(ds)
    has_bounds = bounds is not None and bounds in ds.variables
    for name, da in ds.data_vars.items():
        if not {"x", "y"} <= set(da.dims):
            continue
        if "units" not in da.attrs and "flag_values" not in da.attrs:
            out.append(Finding(WARNING, path, f"{name} has no units"))
        method = da.attrs.get(regrid.METHOD_ATTR)
        if method is not None and method not in get_args(regrid.Method):
            msg = f"{name}: {regrid.METHOD_ATTR} {method!r} is not one of {', '.join(get_args(regrid.Method))}"
            out.append(Finding(ERROR, path, msg))
        if "time" in da.dims:
            out += _cell_methods_findings(path, str(name), da, has_bounds)
    return out


def stations_findings(path: str, ds: xr.Dataset) -> list[Finding]:
    """The ragged array's ``row_size`` accounts for every report, and each
    report variable has units."""
    out = []
    total, reports = int(ds.row_size.sum()), ds.sizes["obs"]
    if total != reports:
        msg = f"row_size adds up to {total}, but obs has {reports} reports"
        out.append(Finding(ERROR, path, msg))
    for name, da in ds.data_vars.items():
        if da.dims == ("obs",) and "units" not in da.attrs:
            out.append(Finding(WARNING, path, f"{name} has no units"))
    return out


def source_findings(path: str, node: xr.DataTree) -> list[Finding]:
    missing = [a for a in SOURCE_ATTRS if a not in node.attrs]
    if not missing:
        return []
    return [Finding(WARNING, path, f"no {', '.join(missing)} (source provenance)")]


def _variables(ds: xr.Dataset, dims: set[str]) -> list[str]:
    rows = []
    for name, da in ds.data_vars.items():
        if not dims <= set(da.dims) or name == _bounds_name(ds):
            continue
        a = da.attrs
        cols = [
            str(name),
            str(da.dtype),
            a.get("units", ""),
            a.get("cell_methods", ""),
            a.get(regrid.METHOD_ATTR, "categorical" if "flag_values" in a else ""),
            a.get(provenance.COMMON_NAME_ATTR, ""),
        ]
        rows.append(cols)
    if not rows:
        return []
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    return [
        "    " + "  ".join(c.ljust(w) for c, w in zip(r, widths, strict=True)).rstrip()
        for r in rows
    ]


def _time_summary(ds: xr.Dataset) -> list[str]:
    if "time" not in ds.coords or ds.time.size == 0:
        return []
    if not np.issubdtype(ds.time.dtype, np.datetime64):
        return []
    t = np.sort(ds.time.values.ravel())
    line = f"  time: {t.size} step(s), {_time(t[0])} to {_time(t[-1])}"
    if ds.time.dims == ("time",) and t.size > 1:
        line += f", every {_duration(np.median(np.diff(t)))}"
    bounds = _bounds_name(ds)
    if bounds in ds.variables and ds[bounds].shape[1:] == (2,):
        b = ds[bounds].values[0] - ds.time.values[0]
        line += (
            f", each covering ({_duration(b[0])}, {_duration(b[1])}] around its time"
        )
    elif ds.time.dims == ("time",):
        line += ", instantaneous (no time bounds)"
    return [line]


def _raster_summary(ds: xr.Dataset) -> str:
    line = f"{ds.sizes['y']} x {ds.sizes['x']} cells"
    crs = ds.rio.crs if "spatial_ref" in ds.coords else None
    if crs is None:
        return line
    rx, ry = ds.rio.resolution()
    p = pyproj.CRS.from_user_input(crs)
    unit = p.axis_info[0].unit_name
    unit = {"metre": "m", "degree": "°"}.get(unit.lower(), f" {unit}")
    epsg = crs.to_epsg()
    name = f"EPSG:{epsg}" if epsg else " ".join(p.name.split())
    if name.lower() in ("unknown", "unnamed"):
        name = f"custom {p.type_name.split()[0].lower()} CRS"
    return f"{line} of {abs(rx):g}{unit} x {abs(ry):g}{unit}, {name}"


def describe(path: str, node: xr.DataTree) -> list[str]:
    """Summary lines for the group at *path*."""
    ds = node.to_dataset()
    k = kind(node)
    if k == "raster":
        head = f"{path}: raster, {_raster_summary(ds)}"
        if ds.attrs.get(provenance.SCENES_ATTR):
            head += ", scenes (gaia align keeps their times)"
        return [head, *_time_summary(ds), *_variables(ds, {"x", "y"})]
    if k == "stations":
        head = f"{path}: stations, {ds.sizes[ds.row_size.dims[0]]} station(s), {ds.sizes['obs']} report(s)"
        return [head, *_time_summary(ds), *_variables(ds, {"obs"})]
    if k == "waveforms":
        sites, channels = (node[t].dataset for t in WAVEFORM_TABLES)
        rates = [
            f"{name} ({child.dataset.sizes.get('channel', 0)} channel(s))"
            for name, child in node.children.items()
            if name not in WAVEFORM_TABLES
        ]
        head = (
            f"{path}: waveforms, {max(sites.sizes.values(), default=0)} site(s), "
            f"{max(channels.sizes.values(), default=0)} channel(s); {', '.join(rates)}"
        )
        return [head]
    dims = ", ".join(f"{d}={n}" for d, n in ds.sizes.items())
    return [f"{path}: {dims} (not a raster: gaia align ignores it)"]


def check(tree: xr.DataTree) -> list[tuple[str, list[str], list[Finding]]]:
    """``(path, summary lines, findings)`` for each group of *tree* holding
    data. The children of a waveforms group are checked but described with it."""
    reports = []
    in_waveforms: set[str] = set()
    for path, node in tree.subtree_with_keys:
        if path == ".":
            continue
        k = kind(node)
        if k == "waveforms":
            in_waveforms.update(f"{path}/{c}" for c in node.children)
        ds = node.to_dataset()
        if not ds.data_vars and k is None:
            continue
        findings = time_findings(path, ds)
        if k == "raster":
            findings += raster_findings(path, ds)
        elif k == "stations":
            findings += stations_findings(path, ds)
        if k is not None:
            findings += source_findings(path, node)
        if path in in_waveforms:
            if findings:
                reports.append((path, [], findings))
            continue
        reports.append((path, describe(path, node), findings))
    return reports


def _header(store: str, tree: xr.DataTree) -> list[str]:
    a = tree.attrs
    lines = [f"{store}" + (f": {a['title']}" if "title" in a else "")]
    if "gaia:cli-version" in a:
        made = f"  written by gaia-cli {a['gaia:cli-version']}"
        if "created" in a:
            made += f" at {a['created']}"
        lines.append(made)
        if "gaia:command" in a:
            lines.append(f"  command: {a['gaia:command']}")
    else:
        lines.append("  no gaia:cli-version: not written by gaia stage or gaia align")
    return lines


def run(
    input_path: Annotated[str, Parameter(name=["--input", "-i"])],
    strict: Annotated[bool, Parameter(name=["--strict"], negative="")] = False,
):
    """Summarize a store's groups and report what gaia align can't use.

    Prints each group's kind, grid, time steps and variables (dtype, units,
    cell_methods, regrid method, common name), followed by its problems.
    Errors make ``gaia align`` skip the group or give wrong values; warnings
    are missing metadata. Works on any Zarr v3 store, not only one written by
    ``gaia stage``. Exits 1 if there are errors.

    Parameters
    ----------
    input_path : str
        Zarr store (local path or ``s3://`` URL).
    strict : bool
        Exit 1 on warnings too. Default False.

    Examples
    --------
        $ gaia check skagit/skagit.zarr
        $ gaia check s3://gaia-scratch/gaia-cli/github-workflows/pioneer/pioneer.zarr --strict
    """
    try:
        tree = io.open_datatree(input_path)
    except Exception as e:  # noqa: BLE001  any read failure is the answer
        sys.exit(f"Error: could not open {input_path} as a Zarr v3 store: {e}")
    lines = _header(input_path, tree)
    errors = warnings = 0
    for _, summary, findings in check(tree):
        lines += summary
        for f in findings:
            prefix = "  " if summary else f"  {f.group}: "
            lines.append(f"{prefix}{f.level.upper()}: {f.message}")
        errors += sum(f.level == ERROR for f in findings)
        warnings += sum(f.level == WARNING for f in findings)
    lines.append(f"{errors} error(s), {warnings} warning(s)")
    print("\n".join(lines))
    if errors or (strict and warnings):
        sys.exit(1)
