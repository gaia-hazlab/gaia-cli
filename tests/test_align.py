import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import align, io


def _raster(times, res, nx, ny, values):
    x = 560_000 + res * (np.arange(nx) + 0.5)
    y = 5_365_000 - res * (np.arange(ny) + 0.5)
    da = xr.DataArray(
        values,
        coords={"time": times, "y": y, "x": x},
        dims=("time", "y", "x"),
        name="precipitation",
    )
    return da.rio.write_crs("EPSG:32610")


def _hourly(values_per_hour=1.0):
    # hour-ending values from 2025-12-01T13 to 2025-12-03T12 → two full PRISM days
    times = pd.date_range("2025-12-01T13:00", "2025-12-03T12:00", freq="h")
    da = _raster(times, 3_000, 2, 2, np.full((len(times), 2, 2), values_per_hour))
    return io.with_time_bounds(da, -1, 0, "time: sum")


def _daily():
    times = pd.date_range("2025-12-01", "2025-12-03", freq="D")
    da = _raster(times, 1_000, 6, 6, np.full((3, 6, 6), 10.0))
    return io.with_time_bounds(da, -12, 12, "time: sum")


def _day_bins(daily):
    b = daily.time_bounds.values
    return daily.time.values, b[:, 0], b[:, 1]


def test_hours_land_in_prism_day_ending_12z():
    hourly = _hourly(0.0)
    t = pd.DatetimeIndex(hourly.time.values)
    hourly["precipitation"][t == pd.Timestamp("2025-12-02T12:00")] = (
        1.0  # last hour of Dec 2
    )
    hourly["precipitation"][t == pd.Timestamp("2025-12-02T13:00")] = (
        5.0  # first hour of Dec 3
    )

    out = align.aggregate_time(hourly, *_day_bins(_daily()))

    totals = out.precipitation.isel(x=0, y=0).values
    # Dec 1 needs hours from Nov 30 13Z → incomplete → NaN, not an undercount
    assert np.isnan(totals[0])
    np.testing.assert_allclose(totals[1:], [1.0, 5.0])
    assert out.precipitation.attrs["cell_methods"] == "time: sum"


def test_point_values_are_averaged():
    times = pd.date_range("2025-12-01T13:00", "2025-12-02T12:00", freq="h")
    temp = _raster(times, 3_000, 1, 1, np.arange(24.0).reshape(24, 1, 1)).rename("t2m")
    temp.attrs["cell_methods"] = "time: point"

    out = align.aggregate_time(temp.to_dataset(), *_day_bins(_daily()))

    np.testing.assert_allclose(out.t2m.sel(time="2025-12-02").squeeze(), 11.5)


def test_align_daily_on_coarsest_grid_roundtrip(tmp_path):
    dt = xr.DataTree.from_dict({"rasters/prism": _daily(), "rasters/hrrr": _hourly()})
    path = str(tmp_path / "native.zarr")
    io.save_datatree(dt, path)

    aligned = align.align(io.open_datatree(path))
    io.save_datatree(aligned, str(tmp_path / "aligned.zarr"))
    reopened = io.open_datatree(str(tmp_path / "aligned.zarr"))

    prism = reopened["rasters/prism"].to_dataset()
    hrrr = reopened["rasters/hrrr"].to_dataset()
    # HRRR (3 km) is the coarsest grid, PRISM (daily) the coarsest time step
    np.testing.assert_array_equal(prism.x, hrrr.x)
    np.testing.assert_array_equal(prism.time, hrrr.time)
    np.testing.assert_allclose(prism.precipitation, 10.0)
    np.testing.assert_allclose(hrrr.precipitation.sel(time="2025-12-02"), 24.0)
    for ds in (prism, hrrr):
        assert "spatial_ref" in ds.coords
        assert ds.precipitation.rio.crs.to_epsg() == 32610
    assert prism.precipitation.attrs["regrid_method"] == "conservative"


def test_clip_to_aoi_keeps_native_crs_and_touched_cells():
    da = _raster(
        pd.date_range("2025-12-01", periods=1), 1_000, 10, 10, np.ones((1, 10, 10))
    )
    # AOI in EPSG:4326 covering a small part of the UTM grid
    aoi_utm = gpd.GeoDataFrame(
        geometry=[box(562_500, 5_360_500, 565_500, 5_363_500)], crs="EPSG:32610"
    )

    clipped = io.clip_to_aoi(da, aoi_utm.to_crs("EPSG:4326"))

    assert clipped.rio.crs.to_epsg() == 32610
    # edges at .5 km offsets only touch cells: all_touched keeps 4 x 4 instead of 3 x 3
    assert clipped.notnull().sum().item() >= 16
