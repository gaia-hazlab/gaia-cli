"""``gaia align`` – put staged rasters on a common grid and time step.

Reads a store written by ``gaia stage`` (native grids), aggregates every
raster in time to the coarsest product's intervals, then regrids once onto
the coarsest product's grid, and writes a *new* store. See dev/reprojecting.md.

:func:`common_window` and :func:`select_window` pick the time steps that a
comparison over a requested interval can use.
"""

from __future__ import annotations

import functools
import logging
import sys
from collections.abc import Mapping
from typing import Annotated, cast, get_args

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from cyclopts import Parameter
from shapely.geometry import box

from . import io, provenance, regrid
from .aoi import load_aoi
from .window import parse_time

log = logging.getLogger(__name__)

# degC per km: the ICAO standard atmosphere's, for --elevation
LAPSE_RATE = 6.5


def _is_raster(ds: xr.Dataset) -> bool:
    return {"x", "y"} <= set(ds.dims) and "spatial_ref" in ds.coords


def _raster_vars(ds: xr.Dataset) -> list[str]:
    return [str(v) for v in ds.data_vars if {"x", "y"} <= set(ds[v].dims)]


def _bounds_name(ds: xr.Dataset) -> str | None:
    name = ds.time.attrs.get("bounds") or ds.time.encoding.get("bounds")
    return name if name in ds.variables else None


def _intervals(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray, bool]:
    """(start, end, is_point) of each time step; point data has start == end."""
    name = _bounds_name(ds)
    if name is None:
        t = ds.time.values
        return t, t, True
    b = ds[name].values
    return b[:, 0], b[:, 1], False


def common_window(tree: xr.DataTree, start, end) -> tuple[pd.Timestamp, pd.Timestamp]:
    """The longest interval inside ``[start, end]`` made of whole time steps
    of every group with time bounds.

    *start* and *end* are parsed like ``gaia stage -s``/``-e``: UTC, and a
    date-only *end* includes that day. For example, PRISM days run 12–12 UTC,
    so a request for 2025-12-01/2025-12-16 gives 2025-12-01T12 to
    2025-12-16T12 when the tree holds PRISM, and the request itself when it
    holds only hourly products. Groups without time bounds (station running
    totals, instantaneous values) don't constrain the window.
    """
    lo, hi = parse_time(start), parse_time(end, end=True)
    starts, ends = [], []
    for node in tree.subtree:
        ds = node.dataset
        if "time" not in ds.dims:
            continue
        s, e, is_point = _intervals(ds)
        if is_point:
            continue
        inside = (s >= lo.to_datetime64()) & (e <= hi.to_datetime64())
        starts.append(pd.DatetimeIndex(s[inside]))
        ends.append(pd.DatetimeIndex(e[inside]))
    if not starts:
        return lo, hi
    # Edges that are a step edge in every product
    s = functools.reduce(pd.DatetimeIndex.intersection, starts)
    e = functools.reduce(pd.DatetimeIndex.intersection, ends)
    first, last = cast(pd.Timestamp, s.min()), cast(pd.Timestamp, e.max())
    if s.empty or e.empty or first >= last:
        raise ValueError(f"no whole time steps of every product lie in {start}/{end}")
    return first, last


def select_window(ds: xr.Dataset, t0, t1) -> xr.Dataset:
    """The time steps of *ds* that lie entirely inside ``[t0, t1]``.

    Point data (no time bounds) keeps the times with ``t0 <= time <= t1``.
    """
    s, e, _ = _intervals(ds)
    t0, t1 = pd.Timestamp(t0).to_datetime64(), pd.Timestamp(t1).to_datetime64()
    return ds.isel(time=(s >= t0) & (e <= t1))


def _step(ds: xr.Dataset) -> np.timedelta64:
    start, end, is_point = _intervals(ds)
    if is_point:
        return (
            np.median(np.diff(ds.time.values))
            if ds.sizes["time"] > 1
            else np.timedelta64(0, "ns")
        )
    return np.median(end - start)


def _time_coords(ds: xr.Dataset, steps: np.ndarray, bins: np.ndarray, n: int) -> dict:
    """*ds*'s other coordinates along ``time``, one value per target interval.

    *steps* are the source steps used and *bins* their target intervals. A
    coordinate is kept when the steps in each interval agree on its value, as
    for a Sentinel-1 pass's orbit when each pass is its own interval; intervals
    with no step get ``""`` (text) or NaN. One that an interval would mix, such
    as two tracks' passes in one interval, is dropped.
    """
    out = {}
    for name, coord in ds.coords.items():
        if name == "time" or coord.dims != ("time",):
            continue
        values = pd.Series(coord.values[steps])
        if (values.groupby(bins).nunique() > 1).any():
            log.info("%s: dropped, it differs within a target time step", name)
            continue
        per_bin = values.groupby(bins).first().reindex(range(n))
        if coord.dtype.kind in "OSTU":
            data = np.asarray(per_bin.fillna("").tolist(), dtype=str)
        elif per_bin.isna().any():
            data = per_bin.to_numpy(dtype=float)
        else:
            data = per_bin.to_numpy(dtype=coord.dtype)
        out[str(name)] = ("time", data, coord.attrs)
    return out


def _overlap(
    start: np.ndarray, end: np.ndarray, target_start, target_end
) -> xr.DataArray:
    """Overlap of each source interval with each target interval, in seconds
    (dims ``step``, ``time``)."""
    s = xr.DataArray(start, dims="time")
    e = xr.DataArray(end, dims="time")
    ts = xr.DataArray(np.asarray(target_start), dims="step")
    te = xr.DataArray(np.asarray(target_end), dims="step")
    inside = np.minimum(e, te) - np.maximum(s, ts)
    return (inside / np.timedelta64(1, "s")).clip(min=0).transpose("step", "time")


def overlap(bounds: xr.DataArray, target_start, target_end) -> xr.DataArray:
    """The fraction of each source interval that lies inside each target interval.

    *bounds* are the source's ``time_bounds`` (``time``, ``nv``); the targets are
    ``(start, end]`` intervals. Gives weights (``step``, ``time``): a source sum
    spread over the targets is ``xr.dot(overlap(...), da, dim="time")``, e.g.
    hourly precipitation between satellite passes, or a daily total spread
    evenly over its hours.
    """
    start, end = bounds.isel(nv=0).values, bounds.isel(nv=1).values
    seconds = (end - start) / np.timedelta64(1, "s")
    return _overlap(start, end, target_start, target_end) / xr.DataArray(
        seconds, dims="time"
    )


def aggregate_time(
    ds: xr.Dataset,
    labels: np.ndarray,
    target_start: np.ndarray,
    target_end: np.ndarray,
) -> xr.Dataset:
    """Aggregate *ds* into the target intervals ``(start, end]``; variables
    without time pass through.

    Each source interval is weighted by how much of it lies inside each target
    (:func:`overlap`): a source step that straddles a target edge is split
    between the two by time, and a longer source step is spread evenly over the
    shorter targets it covers. Point values stand for the step ending at them.
    Values combine according to the variable's ``cell_methods``: ``time: sum``
    as a weighted sum, ``time: minimum``/``maximum`` over every overlapping
    step, otherwise a mean weighted by overlap. Missing values within the
    source's time span are skipped, so a sum counts the data present; a target
    that reaches beyond the span (e.g. a day starting before an hourly source
    does), or that has no data, is NaN.
    """
    start, end, is_point = _intervals(ds)
    step = _step(ds)
    if is_point:  # each point stands for the step ending at it
        start = end - step
    seconds = xr.DataArray((end - start) / np.timedelta64(1, "s"), dims="time")
    if step > np.timedelta64(0, "ns"):
        covers = _overlap(start, end, target_start, target_end)
    else:  # one time step: it counts wholly for the target it falls in
        ts = xr.DataArray(target_start, dims="step")
        te = xr.DataArray(target_end, dims="step")
        t = xr.DataArray(end, dims="time")
        covers = ((ts < t) & (t <= te)).astype(float)
        seconds = xr.ones_like(seconds)

    # Targets that reach beyond the source's time span are NaN, not partial
    inside_span = xr.DataArray(
        (target_start >= start.min()) & (target_end <= end.max()), dims="step"
    )

    # A source step's coordinates go to the target holding its end
    j = np.searchsorted(target_end, end, side="left")
    jj = np.minimum(j, len(target_end) - 1)
    held = np.flatnonzero((j < len(target_end)) & (target_start[jj] < end))

    bounds = _bounds_name(ds)
    timed = [str(v) for v in ds.data_vars if "time" in ds[v].dims and v != bounds]
    out = {}
    for name in timed:
        da = ds[name]
        method = da.attrs.get("cell_methods", "time: mean")
        if not np.issubdtype(da.dtype, np.floating):
            da = da.astype(np.float32)
        present = da.notnull()
        has_data = inside_span & (xr.dot(covers, present, dim="time") > 0)
        if "minimum" in method or "maximum" in method:
            # extremes can't be spread: each target takes them over every
            # source step it overlaps
            reduce = "min" if "minimum" in method else "max"
            agg = xr.concat(
                [getattr(da.where(row > 0), reduce)("time") for row in covers],
                dim="step",
            )
        elif "sum" in method:
            agg = xr.dot(covers / seconds, da.fillna(0), dim="time")
        else:
            weight = covers * present
            agg = xr.dot(covers, da.fillna(0), dim="time") / xr.dot(
                weight, xr.ones_like(da), dim="time"
            )
        agg = agg.where(has_data).rename(step="time").transpose(*da.dims)
        agg.attrs = {
            **ds[name].attrs,
            "cell_methods": method.replace("point", "mean"),
            "aggregated_from": f"steps of {np.timedelta64(step, 'm')}, "
            "weighted by their overlap with each target step",
        }
        out[name] = agg.astype(da.dtype, copy=False)

    # Variables without time (e.g. ERA5's orography) pass through unchanged
    for name in ds.data_vars:
        if "time" not in ds[name].dims and name not in out:
            out[str(name)] = ds[name]
    result = xr.Dataset(out, attrs=ds.attrs).assign_coords(time=labels)
    result = result.assign_coords(_time_coords(ds, held, jj[held], len(labels)))
    result["time_bounds"] = (("time", "nv"), np.stack([target_start, target_end], -1))
    result.time.attrs["bounds"] = "time_bounds"
    if "spatial_ref" in ds.coords:
        result = result.assign_coords(spatial_ref=ds.spatial_ref)
    if "start_datetime" in ds.attrs:
        result.attrs.update(provenance.time_extent(result))
    return result


def _template(ds: xr.Dataset) -> xr.DataArray:
    """A 2D (y, x) DataArray carrying *ds*'s grid."""
    da = ds[_raster_vars(ds)[0]]
    return da.isel({d: 0 for d in da.dims if d not in ("y", "x")}, drop=True)


def _extent(rasters: dict[str, xr.Dataset], crs) -> np.ndarray:
    """(minx, miny, maxx, maxy) in *crs* covering every raster."""
    bounds = np.array(
        [_template(ds).rio.transform_bounds(crs) for ds in rasters.values()]
    )
    return np.concatenate([bounds[:, :2].min(axis=0), bounds[:, 2:].max(axis=0)])


def utm_crs(rasters: dict[str, xr.Dataset], aoi=None) -> str:
    """The AOI's own CRS when it is UTM, else the UTM zone of its center (of the
    rasters' extent without an AOI)."""
    if aoi is None:
        aoi = gpd.GeoDataFrame(geometry=[box(*_extent(rasters, 4326))], crs=4326)
    crs = aoi.crs if aoi.crs.utm_zone else aoi.estimate_utm_crs()
    return crs.to_string()


def _grid_from_crs_res(
    rasters: dict[str, xr.Dataset], crs: str, res: float, aoi=None
) -> xr.DataArray:
    """A grid with *res* spacing in *crs* covering *aoi* (else every raster),
    its edges snapped outward to multiples of *res*."""
    if aoi is not None:
        minx, miny, maxx, maxy = aoi.to_crs(crs).total_bounds
    else:
        minx, miny, maxx, maxy = _extent(rasters, crs)
    # 1e-6 of a cell absorbs float noise in bounds already on a multiple of res
    minx, miny = (np.floor(v / res + 1e-6) * res for v in (minx, miny))
    maxx, maxy = (np.ceil(v / res - 1e-6) * res for v in (maxx, maxy))
    x = minx + res * (np.arange(round((maxx - minx) / res)) + 0.5)
    y = maxy - res * (np.arange(round((maxy - miny) / res)) + 0.5)
    da = xr.DataArray(
        np.zeros((len(y), len(x))), coords={"y": y, "x": x}, dims=("y", "x")
    )
    return da.rio.write_crs(crs)


def _method(
    method: regrid.Method | Mapping[str, regrid.Method] | None,
    path: str,
    da: xr.DataArray,
) -> regrid.Method:
    """*da*'s resampling: named for its group or variable, else categorical
    ``mode``, else the default for every variable, else its dataset's own."""
    methods = method if isinstance(method, Mapping) else {"*": method}
    named = methods.get(f"{path}/{da.name}") or methods.get(path)
    if named:
        return named
    if "flag_values" in da.attrs:
        return "mode"
    return methods.get("*") or da.attrs.get(regrid.METHOD_ATTR, "nearest")


def _cell_area(da: xr.DataArray) -> float:
    """Median cell area (m²) of *da*'s grid, from a few cells at its center."""
    cy, cx = da.sizes["y"] // 2, da.sizes["x"] // 2
    center = da.isel(y=slice(max(cy - 1, 0), cy + 1), x=slice(max(cx - 1, 0), cx + 1))
    return float(np.median(regrid.cell_areas(center)))


def _same_grid(a: xr.DataArray, b: xr.DataArray) -> bool:
    return (
        a.rio.crs == b.rio.crs
        and a.sizes["x"] == b.sizes["x"]
        and a.sizes["y"] == b.sizes["y"]
        and np.allclose(a.x, b.x)
        and np.allclose(a.y, b.y)
    )


def align(
    dt: xr.DataTree,
    *,
    like: str | None = None,
    time_like: str | None = None,
    crs: str | None = None,
    res: float | None = None,
    method: regrid.Method | Mapping[str, regrid.Method] | None = None,
    aoi=None,
    align_scenes: bool = False,
    elevation: str | None = None,
    lapse_rate: float = LAPSE_RATE,
) -> xr.DataTree:
    """Aggregate in time, then regrid once, every raster group in *dt*.

    Non-raster groups (e.g. station data) pass through unchanged. Scene groups
    (``gaia:scenes``, e.g. Sentinel-1/2) keep their acquisition times and are
    only regridded, unless *align_scenes*.

    Parameters
    ----------
    dt : xarray.DataTree
        Store written by ``gaia stage`` with rasters on their native grids.
    like : str, optional
        Group whose grid is the target. Default: the coarsest raster.
    time_like : str, optional
        Group whose time intervals are the target, a scene group included (its
        intervals run from one acquisition to the next). Default: the coarsest
        group that isn't scenes.
    crs, res : optional
        Explicit target grid instead of ``like``: *res* spacing in *crs*, over
        *aoi*'s bounds (else every raster's), snapped outward to multiples of
        *res*. ``crs="utm"`` picks the AOI's own CRS if it is UTM, else the UTM
        zone of its center.
    method : {"conservative", "average", "bilinear", "cubic", "nearest", "mode"} or dict, optional
        Spatial resampling for every continuous variable. Default: each
        variable's ``gaia:regrid-method`` attribute (its dataset's default,
        e.g. ``conservative`` for PRISM and HRRR, ``bilinear`` for DEMs), else
        ``nearest``, as in GDAL and rasterio. ``conservative`` becomes
        ``nearest`` onto a finer grid than the source's, where the two differ
        only along source cell edges. Categorical variables (with
        ``flag_values``) use ``mode``. A dict maps groups or group/variable
        paths to methods, categorical or not, with ``"*"`` for all others,
        e.g. ``{"dem/usgs10m": "cubic", "*": "bilinear"}``.
    aoi : geopandas.GeoDataFrame, optional
        If given, cells whose center is outside it in the target CRS are NaN.
        The grid keeps its extent.
    align_scenes : bool
        Aggregate scene groups in time too, like every other group.
    elevation : str, optional
        A DEM group (its ``elevation``) or ``group/variable``. Every air
        temperature in a group with ``orography`` (ERA5's ``t2m``) gets a
        ``<name>_adjusted`` copy, moved from the source's terrain height to the
        DEM's on the target grid with a constant *lapse_rate*. ``orography`` is
        regridded with the temperature's own method for this, so the result is
        the regridded sea-level temperature lapsed down to the DEM.
    lapse_rate : float
        degC per km, for *elevation*. Default 6.5, the standard atmosphere's.
    """
    nodes = {path: node.to_dataset() for path, node in dt.subtree_with_keys}
    rasters = {p: ds for p, ds in nodes.items() if _is_raster(ds)}
    if not rasters:
        raise ValueError("No raster groups (with x, y and spatial_ref) found.")
    if elevation is not None:
        dem_group, dem_name = _elevation(rasters, elevation)
        if not any(_to_adjust(ds) for ds in rasters.values()):
            raise ValueError(
                f"elevation {elevation!r}: no group has an air temperature and "
                "orography (e.g. ERA5's t2m) to adjust"
            )

    # ---- time ----------------------------------------------------------------
    timed = {p: ds for p, ds in rasters.items() if "time" in ds.dims}
    if time_like is not None and time_like not in timed:
        raise ValueError(f"time_like {time_like!r}: no raster group with time")
    scenes = [
        p
        for p, ds in timed.items()
        if ds.attrs.get(provenance.SCENES_ATTR) and not align_scenes
    ]
    for p in scenes:
        log.info("%s: keeping its acquisition times (scenes)", p)
    aggregated = {p: ds for p, ds in timed.items() if p not in scenes}
    if aggregated:
        time_like = time_like or max(aggregated, key=lambda p: _step(aggregated[p]))
        ref = timed[time_like]
        start, end, is_point = _intervals(ref)
        if is_point:  # each point ends the interval from the previous one
            start = np.r_[end[0] - _step(ref), end[:-1]]
        labels = ref.time.values
        log.info("time target %s (%d step(s))", time_like, len(labels))
        for p, ds in aggregated.items():
            rasters[p] = aggregate_time(ds, labels, start, end)

    # ---- space ---------------------------------------------------------------
    if crs is not None and res is not None:
        if crs.lower() == "utm":
            crs = utm_crs(rasters, aoi)
        template = _grid_from_crs_res(rasters, crs, res, aoi)
        target_desc = f"{crs} at {res}"
    else:
        like = like or max(
            rasters,
            key=lambda p: float(np.median(regrid.cell_areas(_template(rasters[p])))),
        )
        template = _template(rasters[like])
        target_desc = like
    log.info(
        "grid target %s (%d x %d cells)",
        target_desc,
        template.sizes["y"],
        template.sizes["x"],
    )

    _check_size(rasters, template, adjust=elevation is not None)

    # TODO: convert heights to a common vertical datum. DEMs keep their geoid
    # (``vertical_datum``: cop30 EGM2008, usgs10m NAVD88), so differences
    # between them include the geoid separation (~1 m). See dev/dem.md.
    target_area = _cell_area(template)
    for p, ds in rasters.items():
        out = ds.drop_vars(_raster_vars(ds))
        out = out.drop_vars([c for c in out.coords if {"x", "y"} & set(out[c].dims)])
        out = out.drop_vars("spatial_ref", errors="ignore")
        for name in _raster_vars(ds):
            da = ds[name]
            if _same_grid(_template(ds), template):
                out[name] = da
                continue
            m = _method(method, p, da)
            if m == "conservative" and target_area < _cell_area(_template(ds)):
                # Refining: each target cell inside one source cell gets that cell's
                # value either way, so only cells on source edges differ, and
                # nearest avoids a polygon per target cell
                m = "nearest"
            log.info("regridding %s/%s (%s)", p, name, m)
            out[name] = regrid.regrid(da, template, m)
            out[name].attrs["regrid_method"] = m
        if aoi is not None:
            # rio.clip on a Dataset fails on non-spatial vars such as time_bounds
            spatial = _raster_vars(out)
            clipped = out[spatial].rio.clip(
                aoi.to_crs(template.rio.crs).geometry, drop=False
            )
            out = clipped.assign(out.drop_dims(["x", "y"]).data_vars)
        rasters[p] = out

    if elevation is not None:
        dem = rasters[dem_group][dem_name]
        for p, ds in rasters.items():
            for name in _to_adjust(ds):
                m = ds[name].attrs.get("regrid_method")
                orography = (
                    ds["orography"]
                    if m is None  # already on the target grid
                    else regrid.regrid(nodes[p]["orography"], template, m)
                )
                log.info("adjusting %s/%s to %s", p, name, elevation)
                ds[f"{name}_adjusted"] = lapse_rate_adjust(
                    ds[name], orography, dem, lapse_rate
                )

    nodes.update(rasters)
    aligned = xr.DataTree.from_dict(nodes)
    aligned.attrs = {
        **dt.attrs,
        "gaia_align": f"grid={target_desc}; time={time_like}; "
        f"method={method or 'per variable'}; "
        f"scenes={'aligned' if align_scenes else 'kept'}"
        + (f"; elevation={elevation} at {lapse_rate:g} degC/km" if elevation else ""),
    }
    # A title that tells it from its input in a catalog
    title = dt.attrs.get("title")
    if title and not title.endswith(" (aligned)"):
        aligned.attrs["title"] = f"{title} (aligned)"
    return aligned


def _to_adjust(ds: xr.Dataset) -> list[str]:
    """Air temperatures in a group that has its source's terrain height
    (``orography``, e.g. ERA5's), which ``--elevation`` adjusts."""
    if "orography" not in ds:
        return []
    return [
        v
        for v in _raster_vars(ds)
        if ds[v].attrs.get("standard_name") == "air_temperature"
    ]


def _elevation(rasters: dict[str, xr.Dataset], spec: str) -> tuple[str, str]:
    """(group, variable) of ``GROUP`` (its ``elevation``) or ``GROUP/VARIABLE``."""
    if spec in rasters:
        group, name = spec, "elevation"
    else:
        group, _, name = spec.rpartition("/")
    if group not in rasters or name not in rasters[group]:
        raise ValueError(
            f"elevation {spec!r}: no raster group with that variable "
            "(e.g. dem/cop30, or GROUP/VARIABLE)"
        )
    return group, name


def lapse_rate_adjust(
    temperature: xr.DataArray,
    orography: xr.DataArray,
    elevation: xr.DataArray,
    lapse_rate: float = LAPSE_RATE,
) -> xr.DataArray:
    """*temperature*, valid at its source's terrain height *orography*, moved to
    *elevation* (both m) with a constant *lapse_rate* (degC per km)."""
    out = temperature - lapse_rate / 1000 * (elevation - orography)
    # not regridded itself, and one variable per store answers to the common name
    drop = {provenance.COMMON_NAME_ATTR, "regrid_method"}
    attrs = {k: v for k, v in temperature.attrs.items() if k not in drop}
    attrs["long_name"] = (
        f"{temperature.attrs.get('long_name', temperature.name)}, adjusted to elevation"
    )
    attrs["comment"] = (
        f"{temperature.name} - {lapse_rate:g} degC/km x ({elevation.name} - "
        f"{orography.name}): moved from the source's terrain height to the "
        "elevation with a constant lapse rate. An estimate: real lapse rates vary, "
        "and inversions reverse them."
    )
    attrs["gaia:lapse-rate"] = f"{lapse_rate:g} degC/km"
    out = out.rename(f"{temperature.name}_adjusted")
    out.attrs = attrs  # arithmetic keeps all of temperature's
    return out


def _check_size(
    rasters: dict[str, xr.Dataset], template: xr.DataArray, adjust: bool = False
) -> None:
    """Stop before regridding if the rasters on *template* would exceed
    ``io.MAX_BYTES``, naming the largest group. With *adjust*, temperatures
    that ``--elevation`` adjusts count twice."""
    cells = template.sizes["y"] * template.sizes["x"]
    per_group = {
        p: sum(
            float(np.prod([n for d, n in ds[v].sizes.items() if d not in ("x", "y")]))
            * cells
            * ds[v].dtype.itemsize
            * (2 if adjust and v in _to_adjust(ds) else 1)
            for v in _raster_vars(ds)
        )
        for p, ds in rasters.items()
    }
    total = sum(per_group.values())
    if total <= io.MAX_BYTES:
        return
    big = max(per_group, key=lambda p: per_group[p])
    steps = rasters[big].sizes.get("time", 1)
    raise ValueError(
        f"the aligned rasters would be {io.size(total)}, over the "
        f"{io.size(io.MAX_BYTES)} limit; {big} alone is {io.size(per_group[big])} "
        f"({steps} time step(s) on {template.sizes['y']} x {template.sizes['x']} "
        "cells). Use a coarser grid (--like, --res) or time step (--time-like, "
        "e.g. a scene dataset for totals between its passes)."
    )


def _parse_methods(specs: list[str]) -> dict[str, regrid.Method] | None:
    """``["bilinear", "dem/usgs10m=cubic"]`` → ``{"*": "bilinear", "dem/usgs10m": "cubic"}``."""
    methods = {}
    for spec in specs:
        target, _, m = spec.rpartition("=")
        if m not in get_args(regrid.Method):
            raise ValueError(
                f"--method {spec}: {m!r} is not one of {', '.join(get_args(regrid.Method))}"
            )
        methods[target or "*"] = cast(regrid.Method, m)
    return methods or None


def run(
    input_path: Annotated[str, Parameter(name=["--input", "-i"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    like: Annotated[str | None, Parameter(name=["--like"])] = None,
    time_like: Annotated[str | None, Parameter(name=["--time-like"])] = None,
    crs: Annotated[str | None, Parameter(name=["--crs"])] = None,
    res: Annotated[float | None, Parameter(name=["--res"])] = None,
    method: Annotated[list[str] | None, Parameter(name=["--method", "-m"])] = None,
    aoi_path: Annotated[str | None, Parameter(name=["--aoi"])] = None,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    align_scenes: Annotated[
        bool, Parameter(name=["--align-scenes"], negative="")
    ] = False,
    elevation: Annotated[str | None, Parameter(name=["--elevation"])] = None,
    lapse_rate: Annotated[float, Parameter(name=["--lapse-rate"])] = LAPSE_RATE,
):
    """Align staged rasters onto one grid and time step, writing a new store.

    Each raster is first aggregated in time to the coarsest product's
    intervals (using its ``time_bounds`` and ``cell_methods``), then regridded
    once from its native grid onto the coarsest product's grid. Scene datasets
    (Sentinel-1/2) keep their acquisition times and are only regridded.

    Parameters
    ----------
    input_path : str
        Store written by ``gaia stage`` (local path or ``s3://`` URL).
    output_path : str
        Output folder for the aligned store (named after the folder) and its
        STAC Item. The input is never modified.
    like : str
        Group whose grid is the target, e.g. ``rasters/hrrr``. Default: coarsest.
    time_like : str
        Group whose time intervals are the target, e.g. ``rasters/prism``, or a
        scene group such as ``rasters/sentinel1-rtc`` for totals between its
        passes. Default: the coarsest that isn't scenes.
    crs : str
        Target CRS for an explicit grid (use with ``--res``), or ``utm`` for
        the UTM zone of the AOI. The grid covers ``--aoi`` (else every raster),
        its edges snapped outward to multiples of ``--res``.
    res : float
        Target resolution in ``--crs`` units.
    method : list[str]
        Spatial resampling: ``conservative``, ``average``, ``bilinear``,
        ``cubic``, ``nearest`` or ``mode``, for all continuous data, or
        ``GROUP=METHOD`` (or ``GROUP/VARIABLE=METHOD``) for one dataset,
        categorical or not. Repeatable. Default: each dataset's own
        (``conservative`` for PRISM and HRRR, ``bilinear`` for DEMs, ``mode``
        for categories), else ``nearest``. ``conservative`` becomes
        ``nearest`` onto a finer grid.
    aoi_path : str
        Optional vector file; cells whose center is outside it are NaN.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.
    align_scenes : bool
        Aggregate scene datasets (Sentinel-1/2) in time too, onto the same
        intervals as the rest. Default False: they keep their acquisition times.
    elevation : str
        Experimental. A DEM group, e.g. ``dem/cop30`` (or ``GROUP/VARIABLE``):
        adds ``t2m_adjusted``, ERA5's 2 m temperature moved from its smoothed
        terrain (``orography``) to the DEM's elevation on the target grid with
        ``--lapse-rate``. Onto a grid coarser than the DEM, use
        ``-m dem/cop30=average``.
    lapse_rate : float
        degC per km, for ``--elevation``. Default 6.5 (standard atmosphere).

    Examples
    --------
    Daily totals on the HRRR 3 km grid:
        $ gaia stage --config conus-precip -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-native
        $ gaia align -i skagit-native/skagit-native.zarr -o skagit-daily-3km --aoi $AOI

    A 10 m UTM grid over the AOI, with a cubic DEM:
        $ gaia align -i pioneer/pioneer.zarr -o pioneer-10m --aoi $AOI --crs utm --res 10 -m dem/usgs10m=cubic

    ERA5 temperature at the DEM's elevation on a 100 m grid:
        $ gaia align -i peters/peters.zarr -o peters-100m --aoi $AOI --crs utm --res 100 --elevation dem/cop30
    """
    if (crs is None) != (res is None):
        raise ValueError("--crs and --res must be given together.")
    methods = _parse_methods(method or [])
    if output_path.rstrip("/").endswith(".zarr"):
        sys.exit(
            "Error: -o is a folder for <name>.zarr and its STAC Item; "
            f"drop the .zarr from {output_path}"
        )
    io.check_output(output_path, overwrite)
    store = io.store_path(output_path, io.folder_name(output_path))
    dt = io.open_datatree(input_path)
    aoi = load_aoi(aoi_path) if aoi_path else None
    try:
        aligned = align(
            dt,
            like=like,
            time_like=time_like,
            crs=crs,
            res=res,
            method=methods,
            aoi=aoi,
            align_scenes=align_scenes,
            elevation=elevation,
            lapse_rate=lapse_rate,
        )
    except ValueError as e:
        sys.exit(f"Error: {e}")
    aligned.attrs = provenance.derived_from(aligned.attrs, input_path)
    io.save_datatree(aligned, store, overwrite=overwrite)
    log.info("Done")
