import numpy as np
import pandas as pd
import pytest
import xarray as xr
from shapely.geometry import Point

from gaia_cli import io, stations


def _ragged(**reports):
    """A ragged array with one station per keyword: ``{name: {time: {var: value}}}``."""
    ids = list(reports)
    meta = xr.Dataset(
        {"elevation": ("station", np.arange(len(ids), dtype=float))},
        coords={
            "station": ids,
            "geometry": ("station", [Point(-121 - i, 48) for i in range(len(ids))]),
        },
    )
    tables = [
        pd.DataFrame.from_dict(r, orient="index").rename_axis("time")
        for r in reports.values()
    ]
    for t in tables:
        t.index = pd.DatetimeIndex(t.index)
    return stations.ragged(meta, tables)


def _hours(*hours, minute=0, **values):
    return {
        f"2025-12-01T{h:02d}:{minute:02d}": {k: v[i] for k, v in values.items()}
        for i, h in enumerate(hours)
    }


def test_ragged_layout_is_cf_contiguous():
    ds = _ragged(
        A=_hours(0, 1, precipitation=[0.0, 1.0]),
        B=_hours(0, 1, 2, precipitation=[0.0, 2.0, 2.0]),
    )
    assert ds.attrs["featureType"] == "timeSeries"
    assert ds.row_size.values.tolist() == [2, 3]
    assert ds.row_size.attrs["sample_dimension"] == "obs"
    assert ds.station.attrs["cf_role"] == "timeseries_id"
    assert ds.precipitation.dims == ("obs",)
    assert ds.xindexes["geometry"].crs.to_epsg() == 4326


def test_to_dataframe_repeats_station_ids():
    ds = _ragged(
        A=_hours(0, 1, precipitation=[0.0, 1.0]),
        B=_hours(0, precipitation=[5.0]),
    )
    df = stations.to_dataframe(ds)
    assert df.columns.tolist() == ["station", "time", "precipitation"]
    assert df.station.tolist() == ["A", "A", "B"]
    assert df.precipitation.tolist() == [0.0, 1.0, 5.0]


def test_to_grid_hourly_and_quarter_hourly():
    quarter = {
        str(t): {"precipitation": 0.5 * i}
        for i, t in enumerate(pd.date_range("2025-12-01", periods=13, freq="15min"))
    }
    ds = _ragged(
        HOURLY=_hours(0, 1, 2, 3, precipitation=[0.0, 1.0, 1.0, 3.0]),
        QUARTER=quarter,
    )
    grid = stations.to_grid(ds).swap_dims(geometry="station")

    assert dict(grid.sizes) == {"station": 2, "time": 4}
    assert grid.precipitation.sel(station="HOURLY").values.tolist() == [
        0.0,
        1.0,
        1.0,
        3.0,
    ]
    assert grid.precipitation.sel(station="QUARTER").values.tolist() == [
        0.0,
        2.0,
        4.0,
        6.0,
    ]
    assert grid.elevation.values.tolist() == [0.0, 1.0]


def test_to_grid_off_hour_reports_never_include_later_precip():
    ds = _ragged(OFFHOUR=_hours(0, 1, 2, 3, minute=28, precipitation=[0, 1, 3, 6]))
    grid = stations.to_grid(ds).isel(geometry=0)

    # 01:00 holds the 00:28 report, not the 01:28 one
    assert grid.time.dt.strftime("%H:%M").values.tolist() == [
        "01:00",
        "02:00",
        "03:00",
        "04:00",
    ]
    assert grid.precipitation.values.tolist() == [0.0, 1.0, 3.0, 6.0]
    report_time = grid.precipitation_report_time.dt.strftime("%H:%M")
    assert report_time.values.tolist() == ["00:28", "01:28", "02:28", "03:28"]


def test_to_grid_leaves_hours_without_a_report_missing():
    ds = _ragged(GAP=_hours(0, 1, 4, precipitation=[0.0, 1.0, 9.0]))
    grid = stations.to_grid(ds).isel(geometry=0)

    # The 8 mm reported at 04:00 is not assigned to 02:00 or 03:00
    assert grid.sizes["time"] == 5
    assert grid.precipitation.isnull().values.tolist() == [
        False,
        False,
        True,
        True,
        False,
    ]


def test_to_grid_takes_each_variables_own_latest_report():
    # Temperature every 15 minutes, precipitation only on the hour
    reports = {
        str(t): {
            "air_temp": float(i),
            "precipitation": i / 4 if t.minute == 0 else None,
        }
        for i, t in enumerate(pd.date_range("2025-12-01", periods=9, freq="15min"))
    }
    grid = stations.to_grid(_ragged(METAR=reports)).isel(geometry=0)

    assert grid.precipitation.values.tolist() == [0.0, 1.0, 2.0]
    assert grid.air_temp.values.tolist() == [0.0, 4.0, 8.0]


def test_to_grid_sum_gives_interval_totals():
    ds = _ragged(A=_hours(0, 1, 2, 3, precipitation=[0.0, 1.0, 1.0, 3.5]))
    grid = stations.to_grid(ds, how="sum").isel(geometry=0)

    assert grid.precipitation.values.tolist() == [1.0, 0.0, 2.5]
    assert grid.precipitation.attrs["cell_methods"] == "time: sum"
    first = pd.DatetimeIndex(grid.time_bounds.values[0])
    assert first.strftime("%H:%M").tolist() == ["00:00", "01:00"]


def test_to_grid_rejects_unknown_how():
    ds = _ragged(A=_hours(0, precipitation=[0.0]))
    with pytest.raises(ValueError, match="how must be"):
        stations.to_grid(ds, how="mean")  # ty: ignore[invalid-argument-type]


def test_helpers_read_the_stored_group(tmp_path):
    ds = _ragged(
        A=_hours(0, 1, precipitation=[0.0, 1.0]),
        B=_hours(0, 1, precipitation=[0.0, 2.0]),
    )
    ds["name"] = ("geometry", ["first", ""])
    store = str(tmp_path / "s.zarr")
    io.save_datatree(xr.DataTree.from_dict({"stations/x": ds.xvec.encode_cf()}), store)
    stored = io.open_datatree(store)["stations/x"].to_dataset()

    assert stations.to_dataframe(stored).station.tolist() == ["A", "A", "B", "B"]
    grid = stations.to_grid(stored)
    assert grid.name.values.tolist() == ["first", ""]
    assert grid.precipitation.values.tolist() == [[0.0, 1.0], [0.0, 2.0]]
    assert grid.xindexes["geometry"].crs.to_epsg() == 4326
