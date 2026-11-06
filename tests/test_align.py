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


def _hourly_in(year):
    # the same hourly series, but in another year, bounds included
    times = pd.date_range(f"{year}-12-01T13:00", f"{year}-12-03T12:00", freq="h")
    da = _raster(times, 3_000, 2, 2, np.ones((len(times), 2, 2)))
    da.attrs[regrid.METHOD_ATTR] = "conservative"
    return io.with_time_bounds(da, -1, 0, "time: sum")


def test_no_overlap_gives_all_nan_not_an_error():
    # e.g. multitemporal DEM strips from 2012-2020 next to a 2025 event window
    old = _hourly_in(2020)
    old["precipitation"] = old.precipitation.astype("float32")

    out = align.aggregate_time(old, *_day_bins(_daily()))

    assert out.precipitation.sizes == {"time": 3, "y": 2, "x": 2}
    assert out.precipitation.dims == old.precipitation.dims
    assert np.isnan(out.precipitation).all()
    assert out.precipitation.dtype == old.precipitation.dtype
    assert out.precipitation.attrs["cell_methods"] == "time: sum"
    assert "aggregated_from" in out.precipitation.attrs


def test_align_keeps_a_group_outside_the_time_window():
    old = _hourly_in(2020)
    dt = xr.DataTree.from_dict(
        {"rasters/prism": _daily(), "rasters/hrrr": _hourly(), "rasters/strips": old}
    )

    aligned = align.align(dt, time_like="rasters/prism")

    strips = aligned["rasters/strips"].to_dataset()
    np.testing.assert_array_equal(strips.time, aligned["rasters/prism"].time)
    assert np.isnan(strips.precipitation).all()
    np.testing.assert_allclose(aligned["rasters/prism"].precipitation, 10.0)


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
    with pytest.raises(SystemExit, match="Error: time_like 'rasters/x'"):
        align.run(native, str(tmp_path / "x"), time_like="rasters/x")


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


def _passes():
    """Four instantaneous Sentinel-1-like passes on two tracks, with scene coordinates."""
    times = pd.to_datetime(
        ["2023-09-07T16:45", "2023-09-09T16:29", "2023-09-19T16:45", "2023-09-21T16:29"]
    )
    da = _raster(times, 10, 4, 4, np.arange(64.0).reshape(4, 4, 4) + 1).rename("vv")
    da.attrs[regrid.METHOD_ATTR] = "average"
    return da.to_dataset().assign_coords(
        platform=("time", ["SENTINEL-1A"] * 4),
        relative_orbit=("time", np.array([102, 131, 102, 131])),
        scene_ids=("time", ["a", "b", "c", "d e"]),
    )


def test_a_scene_groups_own_steps_keep_its_scene_coordinates(tmp_path):
    s1 = _passes()
    dt = xr.DataTree.from_dict(
        {"rasters/sentinel1-rtc": s1, "dem/cop30": _elevation(30, 2)}
    )
    path = str(tmp_path / "native.zarr")
    io.save_datatree(dt, path)

    aligned = align.align(io.open_datatree(path), crs="EPSG:32610", res=20)
    io.save_datatree(aligned, str(tmp_path / "aligned.zarr"))
    out = io.open_datatree(str(tmp_path / "aligned.zarr"))[
        "rasters/sentinel1-rtc"
    ].to_dataset()

    np.testing.assert_array_equal(out.time, s1.time)
    assert out.relative_orbit.values.tolist() == [102, 131, 102, 131]
    assert out.relative_orbit.dtype.kind == "i"
    assert out.platform.values.tolist() == ["SENTINEL-1A"] * 4
    assert out.scene_ids.values.tolist() == ["a", "b", "c", "d e"]
    assert out.vv.sizes["time"] == 4 and out.vv.notnull().any()


def test_scene_coordinates_that_a_coarser_step_would_mix_are_dropped(caplog):
    s1 = _passes()
    labels = pd.to_datetime(["2023-09-10", "2023-09-22", "2023-10-04"]).values
    start = labels - np.timedelta64(12, "D")

    with caplog.at_level("INFO", logger="gaia_cli.align"):
        out = align.aggregate_time(s1, labels, start, labels)

    # Each 12-day step mixes tracks 102 and 131, so their passes average together
    assert "relative_orbit" not in out.coords
    assert "relative_orbit: dropped" in caplog.text
    # The same platform in each step is kept; the empty last step gets ""
    assert out.platform.values.tolist() == ["SENTINEL-1A", "SENTINEL-1A", ""]
    assert "scene_ids" not in out.coords


def test_numeric_scene_coordinates_with_empty_steps_become_float():
    s1 = _passes().isel(time=[1, 3])  # track 131 only
    labels = pd.to_datetime(["2023-09-10", "2023-09-16", "2023-09-22"]).values
    out = align.aggregate_time(s1, labels, labels - np.timedelta64(6, "D"), labels)
    np.testing.assert_array_equal(out.relative_orbit.values, [131.0, np.nan, 131.0])


def _hours(start, end):
    # hour-ending hourly sums of 1 mm, e.g. ERA5 precipitation
    times = pd.date_range(start, end, freq="h")
    da = _raster(times, 25_000, 2, 2, np.ones((len(times), 2, 2)))
    return io.with_time_bounds(da, -1, 0, "time: sum")


# Sentinel-1-like bins: edges at acquisition times, off the hour; the last one spans a
# missed pass (24 days)
EDGES = pd.to_datetime(
    [
        "2023-08-04T16:29:38",
        "2023-08-16T16:29:38",
        "2023-08-28T16:29:38",
        "2023-09-21T16:29:38",
    ]
).values


def test_bins_off_the_hour_count_the_whole_hours_inside():
    hourly = _hours("2023-08-01T01:00", "2023-09-30T00:00")
    out = align.aggregate_time(hourly, EDGES[1:], EDGES[:-1], EDGES[1:])

    # 12 days hold 287 whole hours between edges at hh:29:38; 24 days hold 575
    assert out.precipitation.isel(y=0, x=0).values.tolist() == [287, 287, 575]
    assert (
        out.precipitation.attrs["aggregated_from"] == "287 to 575 values of 60 minutes"
    )


def test_a_longer_bin_missing_hours_is_nan_not_an_undercount():
    hourly = _hours("2023-08-01T01:00", "2023-09-30T00:00")
    gap = (hourly.time > np.datetime64("2023-09-05")) & (
        hourly.time < np.datetime64("2023-09-07")
    )
    out = align.aggregate_time(
        hourly.isel(time=~gap.values), EDGES[1:], EDGES[:-1], EDGES[1:]
    )

    values = out.precipitation.isel(y=0, x=0).values
    assert values[:2].tolist() == [287, 287]
    assert np.isnan(values[2])


def test_source_ending_inside_a_bin_leaves_it_nan():
    hourly = _hours("2023-08-01T01:00", "2023-08-20T00:00")
    out = align.aggregate_time(hourly, EDGES[1:], EDGES[:-1], EDGES[1:])
    assert out.precipitation.isel(y=0, x=0).values[0] == 287
    assert np.isnan(out.precipitation.isel(y=0, x=0).values[1:]).all()


def test_expected_counts_for_points_and_intervals():
    hourly = _hours("2023-08-01T01:00", "2023-08-03T00:00")
    days = pd.to_datetime(
        ["2023-08-01T00:00", "2023-08-02T00:00", "2023-08-02T12:30"]
    ).values
    starts, ends = days[:-1], days[1:]
    assert align._expected_counts(hourly, starts, ends).tolist() == [24, 12]
    points = hourly.drop_vars("time_bounds")
    del points.time.attrs["bounds"]
    assert align._expected_counts(points, starts, ends).tolist() == [24, 12]


def test_aligned_groups_keep_names_and_select_by_common_name():
    import cf_xarray as cfxr

    from gaia_cli.datasets import cf_criteria

    def staged(ds, name):
        ds = ds.rename(precipitation=name)
        ds[name].attrs["gaia:common-name"] = "precipitation"
        return ds

    tree = xr.DataTree.from_dict(
        {
            "rasters/era5": staged(_hourly(), "tp"),
            "rasters/prism": staged(_daily(), "ppt"),
        }
    )
    aligned = align.align(tree)
    with cfxr.set_options(custom_criteria=cf_criteria()):
        for group, name in (("rasters/era5", "tp"), ("rasters/prism", "ppt")):
            assert aligned[group].to_dataset().cf["precipitation"].name == name


def _scenes():
    s1 = _passes()
    s1.attrs["gaia:scenes"] = True  # as gaia stage writes it for Sentinel-1/2
    return s1


def test_scene_groups_keep_their_acquisition_times():
    s1 = _scenes()
    tree = xr.DataTree.from_dict(
        {
            "rasters/sentinel1-rtc": s1,
            "rasters/era5": _hours("2023-09-01", "2023-09-30"),
        }
    )
    aligned = align.align(tree, crs="EPSG:32610", res=25_000)
    out = aligned["rasters/sentinel1-rtc"].to_dataset()
    np.testing.assert_array_equal(out.time, s1.time)
    assert out.relative_orbit.values.tolist() == [102, 131, 102, 131]
    # ERA5 is the time target, so it keeps its hours
    assert aligned["rasters/era5"].sizes["time"] == 30 * 24 - 23
    assert "scenes=kept" in aligned.attrs["gaia_align"]

    # As before: every timed group on the coarsest, here Sentinel-1's passes
    together = align.align(tree, crs="EPSG:32610", res=25_000, align_scenes=True)
    assert together["rasters/era5"].sizes["time"] == 4
    assert "scenes=aligned" in together.attrs["gaia_align"]


def test_weather_sums_between_passes_when_the_scenes_are_the_target():
    s1 = _scenes()
    hourly = _hours("2023-09-01", "2023-09-30")
    tree = xr.DataTree.from_dict({"rasters/sentinel1-rtc": s1, "rasters/era5": hourly})
    aligned = align.align(
        tree, crs="EPSG:32610", res=25_000, time_like="rasters/sentinel1-rtc"
    )
    np.testing.assert_array_equal(aligned["rasters/sentinel1-rtc"].time, s1.time)
    era5 = aligned["rasters/era5"].to_dataset()
    np.testing.assert_array_equal(era5.time, s1.time)
    # Each interval runs from the previous pass to this one, with its whole hours
    passes = s1.time.values
    np.testing.assert_array_equal(era5.time_bounds.values[1:, 0], passes[:-1])
    ends = hourly.time.values
    whole = [
        int(((ends - np.timedelta64(1, "h") >= lo) & (ends <= hi)).sum())
        for lo, hi in era5.time_bounds.values
    ]
    assert era5.precipitation.isel(y=0, x=0).values.tolist() == whole
    assert whole[1:] == [47, 239, 47]


def test_variables_without_time_pass_through_a_timed_group():
    era5 = _hours("2023-09-01", "2023-09-03")
    era5["orography"] = era5.precipitation.isel(time=0, drop=True) * 0 + 1500.0
    era5.orography.attrs = {"units": "m", "gaia:regrid-method": "bilinear"}
    tree = xr.DataTree.from_dict({"rasters/era5": era5, "rasters/prism": _daily()})
    aligned = align.align(tree, like="rasters/prism")
    out = aligned["rasters/era5"].to_dataset()
    assert out.orography.dims == ("y", "x")
    assert out.orography.attrs["regrid_method"] == "bilinear"
    assert float(out.orography.max()) == 1500.0


def test_an_aligned_store_over_the_limit_stops_before_regridding(monkeypatch):
    # 696 hours on a 10 m grid over 50 km: far over the limit once regridded
    tree = xr.DataTree.from_dict(
        {"rasters/era5": _hours("2023-09-01", "2023-09-30"), "dem/x": _elevation(10, 2)}
    )
    monkeypatch.setattr(io, "MAX_BYTES", 10 * 1024**2)
    with pytest.raises(
        ValueError, match=r"over the 10\.0 MB limit; rasters/era5 alone"
    ):
        align.align(tree, crs="EPSG:32610", res=10)
    # A coarser grid fits
    align.align(tree, crs="EPSG:32610", res=5_000)
