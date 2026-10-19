"""Station observations stored as CF DSG contiguous ragged arrays.

A ``stations/*`` group holds every report, not a dense grid. Station
variables (``station``, ``elevation``, ``row_size``, …) lie along the
``geometry`` dimension, and each report's ``time`` and values along ``obs``:
the first ``row_size[0]`` reports belong to the first station, and so on
(CF ``featureType: timeSeries``). See docs/stations.md.

:func:`to_dataframe` gives the long table and :func:`to_grid` a dense
``(geometry, time)`` view. Both accept a group as read from the store
(xvec CF-encoded) or as returned by a dataset's ``load()``.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
import xarray as xr
import xvec  # noqa: F401  registers the .xvec accessor


def ragged(
    stations: xr.Dataset, reports: list[pd.DataFrame], crs: str = "EPSG:4326"
) -> xr.Dataset:
    """A contiguous ragged array from per-station tables of reports.

    *stations* has one row per station along ``station`` (ids), with a
    ``geometry`` coordinate of shapely points. *reports* holds one DataFrame
    per station, in the same order, with a ``time`` index and a column per
    variable.
    """
    obs = xr.Dataset.from_dataframe(pd.concat(reports).reset_index())
    obs = obs.rename(index="obs")
    obs = obs.drop_vars("obs").set_coords("time")
    sizes = np.array([len(r) for r in reports], dtype="int64")
    ds = xr.merge([stations.assign(row_size=("station", sizes)), obs])
    ds.row_size.attrs = {
        "long_name": "number of reports for each station",
        "sample_dimension": "obs",
    }
    ds.station.attrs["cf_role"] = "timeseries_id"
    ds.attrs["featureType"] = "timeSeries"
    ds = ds.swap_dims(station="geometry")
    return ds.xvec.set_geom_indexes("geometry", crs=crs)


def _decoded(ds: xr.Dataset) -> xr.Dataset:
    return ds.xvec.decode_cf() if "geometry_container" in ds else ds


def _report_vars(ds: xr.Dataset) -> list[str]:
    return [str(v) for v in ds.data_vars if ds[v].dims == ("obs",)]


def to_dataframe(ds: xr.Dataset) -> pd.DataFrame:
    """Every report as a row: ``station``, ``time`` and one column per variable."""
    ds = _decoded(ds)
    df = ds[["time", *_report_vars(ds)]].reset_coords().to_dataframe()
    df.insert(0, "station", np.repeat(ds.station.values, ds.row_size.values))
    return df.reset_index(drop=True)


def to_grid(
    ds: xr.Dataset,
    freq: str = "1h",
    how: Literal["last", "sum"] = "last",
    variables: list[str] | None = None,
) -> xr.Dataset:
    """A dense ``(geometry, time)`` view of station reports on a regular grid.

    For each variable, grid time *t* takes the latest report of it at or
    before *t* and less than one *freq* old; times without such a report are
    NaN. With ``how="last"`` the values are those reports' values and
    ``<variable>_report_time`` records when each was reported, so a station
    reporting at :28 lags the clock hour by 28 minutes (see
    docs/precip-timing.md). ``how="sum"`` is for running
    totals such as ``precipitation``: the amount in each interval
    ``(t - freq, t]``, the difference of the running total at both ends.

    Parameters
    ----------
    ds : xarray.Dataset
        A ``stations/*`` group, CF-encoded or decoded.
    freq : str
        A fixed pandas frequency, e.g. ``"1h"`` or ``"1D"``.
    how : {"last", "sum"}
    variables : list of str, optional
        Report variables to grid. Default: all of them.

    Returns
    -------
    xarray.Dataset
        Station variables as coordinates along ``geometry``, which keeps its
        xvec index.
    """
    ds = _decoded(ds)
    names = list(variables or _report_vars(ds))
    step = pd.Timedelta(freq)
    # Just under one step, so a step with no report of its own is NaN rather
    # than a copy of the previous one
    tolerance = step - pd.Timedelta(seconds=1)
    times = pd.DatetimeIndex(ds.time.values)
    if len(times):
        grid = pd.date_range(
            times.min().ceil(freq), times.max().ceil(freq), freq=freq, name="time"
        )
    else:
        grid = pd.DatetimeIndex([], name="time")

    if how not in ("last", "sum"):
        raise ValueError(f"how must be 'last' or 'sum', not {how!r}")

    ends = np.cumsum(ds.row_size.values)
    shape = (ds.sizes["geometry"], len(grid))
    values = {n: np.full(shape, np.nan) for n in names}
    report_time = {n: np.full(shape, np.datetime64("NaT", "ns")) for n in names}
    table = to_dataframe(ds)
    for i, (start, end) in enumerate(zip(ends - ds.row_size.values, ends, strict=True)):
        # Each variable takes its own latest report: stations may report some
        # variables more often than others
        for n in names:
            rows = table.iloc[start:end].dropna(subset=[n]).set_index("time")[n]
            rows = rows[~rows.index.duplicated(keep="last")].sort_index()
            if rows.empty:
                continue
            asof = rows.index.to_series().reindex(
                grid, method="ffill", tolerance=tolerance
            )
            report_time[n][i] = asof.values
            values[n][i] = rows.reindex(grid, method="ffill", tolerance=tolerance)

    out = ds.drop_dims("obs").drop_vars("row_size")
    out = out.set_coords(list(out.data_vars)).assign_coords(time=grid)
    for n in names:
        out[n] = (("geometry", "time"), values[n], ds[n].attrs)
    out.attrs = {k: v for k, v in ds.attrs.items() if k != "featureType"}
    if how == "last":
        for n in names:
            out[n].attrs["ancillary_variables"] = f"{n}_report_time"
            out.coords[f"{n}_report_time"] = (("geometry", "time"), report_time[n])
        return out
    out = out.diff("time", label="upper")
    for n in names:
        out[n].attrs = {**ds[n].attrs, "cell_methods": "time: sum"}
    out["time_bounds"] = (
        ("time", "nv"),
        np.stack([(out.time - step).values, out.time.values], -1),
    )
    out.time.attrs["bounds"] = "time_bounds"
    return out
