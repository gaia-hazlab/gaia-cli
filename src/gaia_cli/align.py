"""``gaia align`` – put staged rasters on a common grid and time step.

Reads a store written by ``gaia stage`` (native grids), aggregates every
raster in time to the coarsest product's intervals, then regrids once onto
the coarsest product's grid, and writes a *new* store. See docs/reprojecting.md.

:func:`common_window` and :func:`select_window` pick the time steps that a
comparison over a requested interval can use.
"""

from __future__ import annotations

import functools
import logging
import sys
from typing import Annotated

import numpy as np
import pandas as pd
import xarray as xr
from cyclopts import Parameter

from . import io, provenance, regrid

log = logging.getLogger(__name__)


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
    from .stagefile import parse_time  # stagefile imports the dataset modules

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
    if s.empty or e.empty or s.min() >= e.max():
        raise ValueError(f"no whole time steps of every product lie in {start}/{end}")
    return s.min(), e.max()


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


def aggregate_time(
    ds: xr.Dataset,
    labels: np.ndarray,
    target_start: np.ndarray,
    target_end: np.ndarray,
) -> xr.Dataset:
    """Aggregate *ds* into the target intervals ``(start, end]``.

    A source interval contributes to the target interval that fully contains
    it (point values: ``start < t <= end``). Values combine according to the
    variable's ``cell_methods`` (``time: sum`` → sum, ``time: minimum`` → min,
    ``time: maximum`` → max, otherwise mean). Bins
    missing any expected source step become NaN instead of a silent undercount.
    """
    start, end, is_point = _intervals(ds)
    step = _step(ds)
    j = np.searchsorted(target_end, end, side="left")
    jj = np.minimum(j, len(target_end) - 1)
    inside = (j < len(target_end)) & (
        (target_start[jj] < start) if is_point else (target_start[jj] <= start)
    )
    expected = round(float(np.median(target_end - target_start) / step))

    bounds = _bounds_name(ds)
    timed = [str(v) for v in ds.data_vars if "time" in ds[v].dims and v != bounds]
    sub = ds[timed].isel(time=np.flatnonzero(inside))
    sub = sub.assign_coords(_bin=("time", jj[inside]))

    out = {}
    for name in timed:
        da = sub[name]
        method = da.attrs.get("cell_methods", "time: mean")
        grouped = da.groupby("_bin")
        if "sum" in method:
            agg = grouped.sum()
        elif "minimum" in method:
            agg = grouped.min()
        elif "maximum" in method:
            agg = grouped.max()
        else:
            agg = grouped.mean()
        agg = agg.where(grouped.count() >= expected)
        agg = agg.reindex(_bin=np.arange(len(labels))).rename(_bin="time")
        agg.attrs = {
            **da.attrs,
            "cell_methods": method.replace("point", "mean"),
            "aggregated_from": f"{expected} values of {np.timedelta64(step, 'm')}",
        }
        out[name] = agg

    result = xr.Dataset(out, attrs=ds.attrs).assign_coords(time=labels)
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


def _grid_from_crs_res(
    rasters: dict[str, xr.Dataset], crs: str, res: float
) -> xr.DataArray:
    """A grid with *res* spacing in *crs* covering every raster's extent."""
    bounds = np.array(
        [_template(ds).rio.transform_bounds(crs) for ds in rasters.values()]
    )
    minx, miny = bounds[:, :2].min(axis=0)
    maxx, maxy = bounds[:, 2:].max(axis=0)
    x = np.arange(minx + res / 2, maxx, res)
    y = np.arange(maxy - res / 2, miny, -res)
    da = xr.DataArray(
        np.zeros((len(y), len(x))), coords={"y": y, "x": x}, dims=("y", "x")
    )
    return da.rio.write_crs(crs)


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
    method: regrid.Method | None = None,
    aoi=None,
) -> xr.DataTree:
    """Aggregate in time, then regrid once, every raster group in *dt*.

    Non-raster groups (e.g. station data) pass through unchanged.

    Parameters
    ----------
    dt : xarray.DataTree
        Store written by ``gaia stage`` with rasters on their native grids.
    like : str, optional
        Group whose grid is the target. Default: the coarsest raster.
    time_like : str, optional
        Group whose time intervals are the target. Default: the coarsest.
    crs, res : optional
        Explicit target grid instead of ``like``.
    method : {"conservative", "average", "bilinear", "cubic", "nearest"}, optional
        Spatial resampling for every continuous variable. Default: each
        variable's ``gaia:regrid-method`` attribute (its dataset's default,
        e.g. ``bilinear`` for DEMs), else ``conservative``. ``conservative``
        becomes ``nearest`` onto a finer grid than the source's, where the two
        differ only along source cell edges. Categorical variables (with
        ``flag_values``) always use ``mode``.
    aoi : geopandas.GeoDataFrame, optional
        If given, clip the output to it (cell centers) in the target CRS.
    """
    nodes = {path: node.to_dataset() for path, node in dt.subtree_with_keys}
    rasters = {p: ds for p, ds in nodes.items() if _is_raster(ds)}
    if not rasters:
        raise ValueError("No raster groups (with x, y and spatial_ref) found.")

    # ---- time ----------------------------------------------------------------
    timed = {p: ds for p, ds in rasters.items() if "time" in ds.dims}
    if timed:
        time_like = time_like or max(timed, key=lambda p: _step(timed[p]))
        ref = timed[time_like]
        start, end, is_point = _intervals(ref)
        if is_point:  # label each point as the end of a (t - step, t] interval
            start = end - _step(ref)
        labels = ref.time.values
        log.info("time target %s (%d step(s))", time_like, len(labels))
        for p, ds in timed.items():
            rasters[p] = aggregate_time(ds, labels, start, end)

    # ---- space ---------------------------------------------------------------
    if crs is not None and res is not None:
        template = _grid_from_crs_res(rasters, crs, res)
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

    # TODO: convert heights to a common vertical datum. DEMs keep their geoid
    # (``vertical_datum``: cop30 EGM2008, usgs10m NAVD88), so differences
    # between them include the geoid separation (~1 m). See docs/dem.md.
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
            if "flag_values" in da.attrs:
                m = "mode"
            else:
                m = method or da.attrs.get(regrid.METHOD_ATTR, "conservative")
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
            clipped = out[spatial].rio.clip(aoi.to_crs(template.rio.crs).geometry)
            out = clipped.assign(out.drop_dims(["x", "y"]).data_vars)
        rasters[p] = out

    nodes.update(rasters)
    aligned = xr.DataTree.from_dict(nodes)
    aligned.attrs = {
        **dt.attrs,
        "gaia_align": f"grid={target_desc}; time={time_like}; "
        f"method={method or 'per variable'}",
    }
    return aligned


def run(
    input_path: Annotated[str, Parameter(name=["--input", "-i"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    like: Annotated[str | None, Parameter(name=["--like"])] = None,
    time_like: Annotated[str | None, Parameter(name=["--time-like"])] = None,
    crs: Annotated[str | None, Parameter(name=["--crs"])] = None,
    res: Annotated[float | None, Parameter(name=["--res"])] = None,
    method: Annotated[regrid.Method | None, Parameter(name=["--method", "-m"])] = None,
    aoi_path: Annotated[str | None, Parameter(name=["--aoi"])] = None,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
):
    """Align staged rasters onto one grid and time step, writing a new store.

    Each raster is first aggregated in time to the coarsest product's
    intervals (using its ``time_bounds`` and ``cell_methods``), then regridded
    once from its native grid onto the coarsest product's grid.

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
        Group whose time intervals are the target, e.g. ``rasters/prism``. Default: coarsest.
    crs : str
        Target CRS for an explicit grid (use with ``--res``).
    res : float
        Target resolution in ``--crs`` units.
    method : str
        Spatial resampling for all continuous data: ``conservative``,
        ``average``, ``bilinear``, ``cubic`` or ``nearest``. Default: each
        dataset's own (``bilinear`` for DEMs), else ``conservative``, which
        becomes ``nearest`` onto a finer grid.
    aoi_path : str
        Optional vector file; clip the output to it in the target CRS.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.

    Examples
    --------
    Daily totals on the HRRR 3 km grid:
        $ gaia stage --config conus-precip -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-native
        $ gaia align -i skagit-native/skagit-native.zarr -o skagit-daily-3km --aoi $AOI
    """
    if (crs is None) != (res is None):
        raise ValueError("--crs and --res must be given together.")
    if output_path.rstrip("/").endswith(".zarr"):
        sys.exit(
            "Error: -o is a folder for <name>.zarr and its STAC Item; "
            f"drop the .zarr from {output_path}"
        )
    io.check_output(output_path, overwrite)
    store = io.store_path(output_path, io.folder_name(output_path))
    dt = io.open_datatree(input_path)
    aoi = io.load_aoi(aoi_path) if aoi_path else None
    aligned = align(
        dt, like=like, time_like=time_like, crs=crs, res=res, method=method, aoi=aoi
    )
    aligned.attrs = provenance.derived_from(aligned.attrs, input_path)
    io.save_datatree(aligned, store, overwrite=overwrite)
    log.info("Done")
