import json

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from rasterio.features import geometry_mask
from shapely.geometry import box

from gaia_cli import align, io, regrid


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
    da.attrs[regrid.METHOD_ATTR] = "conservative"  # as staged
    return io.with_time_bounds(da, -1, 0, "time: sum")


def _daily():
    times = pd.date_range("2025-12-01", "2025-12-03", freq="D")
    da = _raster(times, 1_000, 6, 6, np.full((3, 6, 6), 10.0))
    da.attrs[regrid.METHOD_ATTR] = "conservative"  # as staged
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


def test_minimum_and_maximum_keep_their_extremes():
    times = pd.date_range("2025-12-01T13:00", "2025-12-02T12:00", freq="h")
    da = _raster(times, 3_000, 1, 1, np.arange(24.0).reshape(24, 1, 1)).rename("tmin")
    ds = io.with_time_bounds(da, -1, 0, "time: minimum")
    ds["tmax"] = ds.tmin.assign_attrs(cell_methods="time: maximum")

    out = align.aggregate_time(ds, *_day_bins(_daily()))

    day = {"time": "2025-12-02"}
    np.testing.assert_allclose(out.tmin.sel(day).squeeze(), 0.0)
    np.testing.assert_allclose(out.tmax.sel(day).squeeze(), 23.0)


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


def _elevation(res, n):
    da = _raster(pd.date_range("2025-12-01", periods=1), res, n, n, np.ones((1, n, n)))
    return da.isel(time=0, drop=True).rename("elevation").to_dataset()


def test_each_variable_defaults_to_its_datasets_method():
    dem = _elevation(1_000, 6)
    dem.elevation.attrs[regrid.METHOD_ATTR] = "bilinear"
    dt = xr.DataTree.from_dict({"dem/x": dem, "rasters/other": _elevation(1_000, 6)})
    template = {"crs": "EPSG:32610", "res": 2_000.0}

    out = align.align(dt, **template)
    assert out["dem/x"].elevation.attrs["regrid_method"] == "bilinear"
    assert out["rasters/other"].elevation.attrs["regrid_method"] == "nearest"
    # --method overrides every dataset's default
    out = align.align(dt, method="cubic", **template)
    assert out["dem/x"].elevation.attrs["regrid_method"] == "cubic"
    assert out["rasters/other"].elevation.attrs["regrid_method"] == "cubic"


def test_conservative_becomes_nearest_onto_a_finer_grid():
    coarse = _elevation(1_000, 4)
    coarse["elevation"][:] = np.arange(16.0).reshape(4, 4)
    coarse.elevation.attrs[regrid.METHOD_ATTR] = "conservative"
    dt = xr.DataTree.from_dict({"rasters/coarse": coarse})

    out = align.align(dt, crs="EPSG:32610", res=250.0)["rasters/coarse"].elevation
    assert out.attrs["regrid_method"] == "nearest"
    # Each 250 m cell inside one 1 km cell keeps its value
    np.testing.assert_array_equal(
        np.unique(out.values[~np.isnan(out.values)]), np.arange(16.0)
    )
    out = align.align(dt, crs="EPSG:32610", res=2_000.0)["rasters/coarse"].elevation
    assert out.attrs["regrid_method"] == "conservative"


def test_run_writes_a_folder_named_store(tmp_path):
    dt = xr.DataTree.from_dict({"rasters/prism": _daily(), "rasters/hrrr": _hourly()})
    dt.attrs["title"] = "Skagit"
    native = str(tmp_path / "native.zarr")
    io.save_datatree(dt, native)
    out = tmp_path / "daily"

    align.run(native, str(out))

    assert sorted(p.name for p in out.iterdir()) == ["daily.json", "daily.zarr"]
    # The Item tells an aligned store from its input
    props = json.loads((out / "daily.json").read_text())["properties"]
    assert props["title"] == "Skagit (aligned)"
    assert props["description"].startswith("gaia aligned store with groups: ")
    native_props = json.loads((tmp_path / "native.json").read_text())["properties"]
    assert native_props["description"].startswith("gaia store with groups: ")
    with pytest.raises(SystemExit, match="drop the .zarr"):
        align.run(native, str(tmp_path / "daily.zarr"))


def test_common_window_snaps_to_every_products_steps():
    tree = xr.DataTree.from_dict({"rasters/prism": _daily(), "rasters/hrrr": _hourly()})
    # PRISM days end at 12 UTC: Nov 30 12Z → Dec 3 12Z; HRRR hours Dec 1 12Z → Dec 3 12Z
    t0, t1 = align.common_window(tree, "2025-12-01", "2025-12-03")
    assert (t0, t1) == (pd.Timestamp("2025-12-01T12"), pd.Timestamp("2025-12-03T12"))
    hrrr = xr.DataTree.from_dict({"rasters/hrrr": _hourly()})
    assert align.common_window(hrrr, "2025-12-02", "2025-12-02T18:00") == (
        pd.Timestamp("2025-12-02"),
        pd.Timestamp("2025-12-02T18:00"),
    )
    with pytest.raises(ValueError, match="no whole time steps"):
        align.common_window(tree, "2025-12-02T00", "2025-12-02T06")


def test_select_window_keeps_whole_steps():
    daily = align.select_window(_daily(), "2025-12-01T12", "2025-12-03T12")
    assert list(daily.time.dt.day.values) == [2, 3]
    hourly = align.select_window(_hourly(), "2025-12-01T12", "2025-12-01T15")
    assert hourly.sizes["time"] == 3


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


def _geographic(values):
    # 0.001° cells near 48°N, 120.5°W (UTM zone 10)
    ny, nx = values.shape
    da = xr.DataArray(
        values,
        coords={
            "y": 48.1 - 0.001 * (np.arange(ny) + 0.5),
            "x": -120.6 + 0.001 * (np.arange(nx) + 0.5),
        },
        dims=("y", "x"),
        name="elevation",
    )
    return da.rio.write_crs("EPSG:4326").to_dataset()


def test_utm_grid_covers_the_aoi_snapped_outward():
    dt = xr.DataTree.from_dict({"dem/x": _geographic(np.ones((60, 60)))})
    aoi = gpd.GeoDataFrame(
        geometry=[box(-120.59, 48.05, -120.555, 48.09)], crs="EPSG:4326"
    )

    out = align.align(dt, crs="utm", res=10.0, aoi=aoi)["dem/x"].elevation

    assert out.rio.crs.to_epsg() == 32610
    left, bottom, right, top = out.rio.bounds()
    for edge in (left, bottom, right, top):
        assert edge % 10 == pytest.approx(0, abs=1e-6)
    minx, miny, maxx, maxy = aoi.to_crs(32610).total_bounds
    assert left <= minx < left + 10 and right - 10 < maxx <= right
    assert bottom <= miny < bottom + 10 and top - 10 < maxy <= top
    # The grid keeps its extent; cells whose center is outside the AOI are NaN
    inside = geometry_mask(
        aoi.to_crs(32610).geometry,
        out_shape=out.shape,
        transform=out.rio.transform(),
        invert=True,
    )
    np.testing.assert_array_equal(out.notnull().values, inside)


def test_utm_keeps_an_aois_own_utm_crs():
    dt = xr.DataTree.from_dict({"dem/x": _geographic(np.ones((10, 10)))})
    aoi = gpd.GeoDataFrame(
        geometry=[box(-120.598, 48.092, -120.592, 48.098)], crs="EPSG:4326"
    ).to_crs("EPSG:26910")  # NAD83 / UTM 10N, not WGS 84's 32610
    assert align.utm_crs({"dem/x": dt["dem/x"].to_dataset()}, aoi) == "EPSG:26910"
    assert align.utm_crs({"dem/x": dt["dem/x"].to_dataset()}) == "EPSG:32610"


def test_methods_per_group_override_defaults_and_categories():
    dem = _elevation(1_000, 6)
    dem.elevation.attrs[regrid.METHOD_ATTR] = "bilinear"
    classes = _elevation(1_000, 6).rename(elevation="landcover")
    classes.landcover.attrs["flag_values"] = [1, 2]
    dt = xr.DataTree.from_dict(
        {"dem/x": dem, "rasters/nlcd": classes, "rasters/other": _elevation(1_000, 6)}
    )
    methods = {"dem/x": "cubic", "rasters/nlcd/landcover": "nearest", "*": "average"}

    out = align.align(dt, crs="EPSG:32610", res=2_000.0, method=methods)

    assert out["dem/x"].elevation.attrs["regrid_method"] == "cubic"
    assert out["rasters/nlcd"].landcover.attrs["regrid_method"] == "nearest"
    assert out["rasters/other"].elevation.attrs["regrid_method"] == "average"
    # "*" leaves categories on mode
    out = align.align(dt, crs="EPSG:32610", res=2_000.0, method={"*": "average"})
    assert out["rasters/nlcd"].landcover.attrs["regrid_method"] == "mode"


def test_parse_methods():
    assert align._parse_methods(["bilinear", "dem/usgs10m=cubic"]) == {
        "*": "bilinear",
        "dem/usgs10m": "cubic",
    }
    assert align._parse_methods([]) is None
    with pytest.raises(ValueError, match="not one of"):
        align._parse_methods(["dem/usgs10m=bicubic"])
