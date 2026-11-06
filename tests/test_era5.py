import os
from types import SimpleNamespace

import arraylake
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from arraylake.exceptions import AuthException
from shapely.geometry import box

from gaia_cli.datasets import era5

STEP_TYPES = {
    "tp": "accum",
    "sf": "accum",
    "t2m": "instant",
    "swvl1": "instant",
    "swvl2": "instant",
    "sd": "instant",
}
UNITS = {
    "tp": "m",
    "sf": "m of water equivalent",
    "t2m": "K",
    "swvl1": "m**3 m**-3",
    "swvl2": "m**3 m**-3",
    "sd": "m of water equivalent",
}


def fake_store(hours=48):
    """A small ERA5-like group as stored: valid_time, latitude 90→-90, longitude 0→360."""
    time = pd.date_range("2023-09-13", periods=hours, freq="h")
    lat = np.arange(64.0, 62.0, -0.25)
    lon = np.arange(207.0, 210.0, 0.25)
    shape = (hours, lat.size, lon.size)
    data = {
        v: xr.DataArray(
            np.full(shape, 0.001 if v != "t2m" else 273.15, "float32"),
            dims=("valid_time", "latitude", "longitude"),
            attrs={"units": UNITS[v], "GRIB_stepType": STEP_TYPES[v], "long_name": v},
        )
        for v in STEP_TYPES
    }
    lsm = xr.DataArray(np.ones(shape[1:], "float32"), dims=("latitude", "longitude"))
    z_sfc = lsm * 1500 * 9.80665  # geopotential of 1500 m terrain
    return xr.Dataset(
        data,
        coords={
            "valid_time": time,
            "latitude": lat,
            "longitude": lon,
            "lsm": lsm,
            "z_sfc": z_sfc,
        },
        attrs={"summary": "ERA5", "acknowledgement": "Copernicus"},
    )


@pytest.fixture
def fake_repo(monkeypatch):
    session = SimpleNamespace(store="store", snapshot_id="SNAP")
    repo = SimpleNamespace(readonly_session=lambda branch: session)
    endpoints = []

    class Client:
        def get_repo(self, name):
            endpoints.append(os.environ.get("AWS_ENDPOINT_URL_S3"))
            return repo

    monkeypatch.setattr(arraylake, "Client", Client)
    monkeypatch.setattr(era5.xr, "open_zarr", lambda store, **kw: fake_store())
    return endpoints


@pytest.fixture
def aoi():
    # Peters Dome, AK
    return gpd.GeoDataFrame(geometry=[box(-151.3, 63.1, -151.1, 63.2)], crs=4326)


def test_open_era5_regularizes_grid(fake_repo):
    with era5._aws_endpoint():
        ds = era5.open_era5()
    assert {"time", "x", "y"} <= set(ds.dims)
    assert float(ds.x.min()) == -153.0 and bool((ds.x.diff("x") > 0).all())
    assert ds.rio.crs.to_epsg() == 4326
    assert ds.attrs["gaia:icechunk-snapshot"] == "SNAP"
    assert fake_repo == [era5.S3_ENDPOINT]


def test_open_era5_explains_missing_token(monkeypatch):
    class Client:
        def get_repo(self, name):
            raise AuthException("not logged in")

    monkeypatch.setattr(arraylake, "Client", Client)
    with pytest.raises(RuntimeError, match="ARRAYLAKE_TOKEN"):
        era5.open_era5()


def test_aws_endpoint_is_restored(monkeypatch):
    monkeypatch.setenv("AWS_ENDPOINT_URL_S3", "https://elsewhere.example")
    with era5._aws_endpoint():
        assert os.environ["AWS_ENDPOINT_URL_S3"] == era5.S3_ENDPOINT
    assert os.environ["AWS_ENDPOINT_URL_S3"] == "https://elsewhere.example"
    monkeypatch.delenv("AWS_ENDPOINT_URL_S3")
    with era5._aws_endpoint():
        pass
    assert "AWS_ENDPOINT_URL_S3" not in os.environ


def test_load_units_and_time(fake_repo, aoi):
    ds = era5.load(aoi, "2023-09-13T00:00", "2023-09-13T23:59")
    assert ds.sizes["time"] == 24
    assert set(ds.data_vars) == {
        "tp",
        "sf",
        "t2m",
        "swvl1",
        "swvl2",
        "sd",
        "time_bounds",
        "orography",
    }
    assert "lsm" not in ds.coords and "z_sfc" not in ds.coords
    # Static terrain height from the surface geopotential
    assert ds.orography.dims == ("y", "x")
    assert float(ds.orography.max()) == pytest.approx(1500, rel=1e-5)
    assert ds.orography.attrs["units"] == "m"
    assert ds.orography.attrs["gaia:regrid-method"] == "bilinear"
    # m → mm, K → degC; names unchanged
    assert ds.tp.attrs["units"] == "mm" and float(ds.tp.max()) == pytest.approx(1.0)
    assert ds.sd.attrs["units"] == "mm"
    assert ds.sf.attrs["cell_methods"] == "time: sum"
    assert ds.t2m.attrs["units"] == "degC" and float(ds.t2m.max()) == pytest.approx(0)
    assert ds.swvl1.attrs["units"] == "m**3 m**-3"
    assert not any(k.startswith("GRIB_") for k in ds.tp.attrs)
    # Hour-ending: (T - 1h, T]
    assert ds.tp.attrs["cell_methods"] == "time: sum"
    assert ds.t2m.attrs["cell_methods"] == "time: point"
    assert ds.time.attrs["bounds"] == "time_bounds"
    assert (ds.time_bounds[:, 1] - ds.time_bounds[:, 0] == np.timedelta64(1, "h")).all()
    assert ds.attrs["start_datetime"] == "2023-09-12T23:00:00Z"
    assert ds.attrs["license"] == "CC-BY-4.0"
    assert ds.attrs["gaia:icechunk-snapshot"] == "SNAP"


def test_load_after_the_store_ends(fake_repo, aoi):
    with pytest.raises(ValueError, match="ERA5 .* ends at 2023-09-14T23:00Z"):
        era5.load(aoi, "2023-10-01", "2023-10-02")
