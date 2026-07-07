import numpy as np
import geopandas as gpd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import io


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


def test_save_zarr_roundtrip(tmp_path):
    ds = xr.Dataset({"precip": (("time",), np.arange(5.0))})
    output_path = tmp_path / "test.zarr"

    io.save_zarr(ds, str(output_path))
    result = xr.open_zarr(output_path, consolidated=False)

    np.testing.assert_array_equal(result["precip"].values, ds["precip"].values)
