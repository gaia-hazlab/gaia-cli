import json

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import io


def test_load_aoi_requires_a_single_polygon(tmp_path):
    path = tmp_path / "two.geojson"
    gpd.GeoDataFrame(
        geometry=[box(0, 0, 1, 1), box(2, 2, 3, 3)], crs="EPSG:4326"
    ).to_file(path)
    with pytest.raises(ValueError, match="single Polygon, not Polygon, Polygon"):
        io.load_aoi(str(path))


def test_check_s3_env_raises_without_aws_profile(monkeypatch):
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    with pytest.raises(ValueError):
        io.check_s3_env("s3://some-bucket/output.zarr")


def test_check_s3_env_local_path_ok():
    io.check_s3_env("/tmp/output.zarr")


def test_clip_to_aoi():
    da = xr.DataArray(
        np.ones((10, 10)),
        dims=("y", "x"),
        coords={"y": np.linspace(1, 0, 10), "x": np.linspace(0, 1, 10)},
    )
    da = da.rio.write_crs("EPSG:4326")
    aoi = gpd.GeoDataFrame(geometry=[box(0.2, 0.2, 0.6, 0.6)], crs="EPSG:4326")

    clipped = io.clip_to_aoi(da, aoi)

    assert clipped.rio.crs.to_epsg() == 4326
    assert clipped.sizes["x"] < da.sizes["x"]


def test_save_datatree_root_roundtrip(tmp_path):
    ds = xr.Dataset({"precip": (("time",), np.arange(5.0))})
    output_path = tmp_path / "test.zarr"

    io.save_datatree(xr.DataTree(ds), str(output_path))
    result = xr.open_zarr(output_path, consolidated=False)

    np.testing.assert_array_equal(result["precip"].values, ds["precip"].values)


def test_io_registers_accessors():
    da = xr.DataArray([1])
    assert hasattr(da, "rio")
    assert hasattr(da, "xvec")


def test_check_output_exits_if_exists(tmp_path):
    output_path = tmp_path / "exists.zarr"
    output_path.mkdir()
    with pytest.raises(SystemExit, match="--overwrite"):
        io.check_output(str(output_path))


def test_check_output_overwrite_ok(tmp_path):
    output_path = tmp_path / "exists.zarr"
    output_path.mkdir()
    io.check_output(str(output_path), overwrite=True)
    io.check_output(str(tmp_path / "new.zarr"))


def test_save_datatree_root_overwrite(tmp_path):
    output_path = str(tmp_path / "test.zarr")
    io.save_datatree(
        xr.DataTree(xr.Dataset({"precip": (("time",), np.arange(5.0))})), output_path
    )
    io.save_datatree(
        xr.DataTree(xr.Dataset({"precip": (("time",), np.zeros(3))})),
        output_path,
        overwrite=True,
    )
    result = xr.open_zarr(output_path, consolidated=False)
    np.testing.assert_array_equal(result["precip"].values, np.zeros(3))


def test_save_datatree_writes_fill_value_chunks(tmp_path):
    # Kopah answers anonymous GETs for missing chunks with 403, so chunks equal
    # to the fill value (HRRR's spatial_ref is 0) must still be written
    da = xr.DataArray(
        np.ones((2, 2)),
        coords={"y": [1.0, 0.0], "x": [0.0, 1.0]},
        dims=("y", "x"),
        name="v",
    ).rio.write_crs("EPSG:32610")
    da = da.assign_coords(spatial_ref=0)
    path = tmp_path / "t.zarr"

    io.save_datatree(xr.DataTree.from_dict({"a": da.to_dataset()}), str(path))

    assert (path / "a" / "spatial_ref" / "c").exists()


def test_save_datatree_append_keeps_groups_and_describes_the_whole_store(tmp_path):
    path = str(tmp_path / "out.zarr")
    a = xr.Dataset({"v": ("time", [1.0])}, attrs={"license": "CC-BY-4.0"})
    b = xr.Dataset({"w": ("time", [2.0])}, attrs={"license": "other"})
    io.save_datatree(xr.DataTree.from_dict({"rasters/a": a}), path)
    io.save_datatree(xr.DataTree.from_dict({"rasters/b": b}), path, append=True)

    tree = xr.open_datatree(path, engine="zarr", consolidated=True)
    assert set(tree["rasters"].children) == {"a", "b"}
    assert float(tree["rasters/a/v"][0]) == 1.0
    # The note covers both groups' licenses, not just the appended one
    assert "gaia:license-note" in tree.attrs
    item = json.loads((tmp_path / "out.json").read_text())
    assert "rasters/a" in item["properties"]["description"]
    assert "rasters/b" in item["properties"]["description"]
