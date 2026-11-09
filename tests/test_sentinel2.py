import datetime as dt

import geopandas as gpd
import numpy as np
import pandas as pd
import pystac
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import zarr
from shapely.geometry import box, mapping

from gaia_cli import io, stac
from gaia_cli.aoi import load_aoi
from gaia_cli.cli import stage_app
from gaia_cli.datasets import sentinel2
from gaia_cli.window import parse_window

# 10 x 10 cells of 10 m in UTM 5N; each tile is half of it
X0, Y0 = 590000, 7010000
NORTH = (X0, Y0 - 50, X0 + 100, Y0)
SOUTH = (X0, Y0 - 100, X0 + 100, Y0 - 50)
COLLECTION = {"title": "Sentinel-2 Level-2A", "license": "proprietary", "links": []}


def tif(tmp_path, name, bounds, value, res, dtype="uint16"):
    """A GeoTIFF of *value* at *res* over *bounds*, with one 0 (nodata) cell."""
    left, bottom, right, top = bounds
    ny, nx = round((top - bottom) / res), round((right - left) / res)
    values = np.full((ny, nx), value, dtype=dtype)
    values[0, 0] = 0
    da = xr.DataArray(
        values,
        dims=("y", "x"),
        coords={
            "y": np.arange(top - res / 2, bottom, -res),
            "x": np.arange(left + res / 2, right, res),
        },
    )
    path = tmp_path / f"{name}.tif"
    da.rio.write_crs("EPSG:32605").rio.write_nodata(0).rio.to_raster(path)
    return str(path)


def item(
    tmp_path, name, bounds, when, datatake, tile, baseline="05.10", dn=2000, cloud=10.0
):
    footprint = gpd.GeoSeries([box(*bounds)], crs="EPSG:32605").to_crs(4326)[0]
    it = pystac.Item(
        id=name,
        geometry=mapping(footprint),
        bbox=list(footprint.bounds),
        datetime=dt.datetime.fromisoformat(when),
        properties={
            "platform": "Sentinel-2A",
            "sat:relative_orbit": 129,
            "eo:cloud_cover": cloud,
            "s2:mgrs_tile": tile,
            "s2:datatake_id": f"{datatake}_N{baseline}",
            "s2:processing_baseline": baseline,
            "s2:generation_time": f"{when}Z",
            "proj:epsg": 32605,
        },
    )
    assets = {"B02": (dn, 10), "B11": (dn + 1000, 20), "SCL": (4, 20)}
    for key, (value, res) in assets.items():
        dtype = "uint8" if key == "SCL" else "uint16"
        href = tif(tmp_path, f"{name}-{key}", bounds, value, res, dtype)
        it.add_asset(key, pystac.Asset(href, media_type=pystac.MediaType.COG))
    return it


@pytest.fixture
def tiles(tmp_path):
    """One pass on two tiles (one also in an older processing), and a 2021 pass."""
    return [
        item(
            tmp_path,
            "a-north",
            NORTH,
            "2023-09-12T21:45:29",
            "GS2B_20230912_034046",
            "05VNL",
            cloud=40.0,
        ),
        item(
            tmp_path,
            "a-south",
            SOUTH,
            "2023-09-12T21:45:29",
            "GS2B_20230912_034046",
            "05VNK",
            cloud=20.0,
        ),
        item(
            tmp_path,
            "a-south-old",
            SOUTH,
            "2023-09-12T21:45:29",
            "GS2B_20230912_034046",
            "05VNK",
            baseline="05.09",
            dn=9000,
        ),
        item(
            tmp_path,
            "b",
            NORTH,
            "2023-09-19T21:35:29",
            "GS2B_20230919_034146",
            "05VNL",
            baseline="03.01",
        ),
    ]


@pytest.fixture
def pc(monkeypatch, tiles):
    """Planetary Computer stand-ins; returns the searches made."""
    searches = []

    class Client:
        def search(self, **kwargs):
            searches.append(kwargs)
            return type("Search", (), {"items": lambda self: iter(tiles)})()

    monkeypatch.setattr(stac.pystac_client.Client, "open", lambda url: Client())
    monkeypatch.setattr(stac, "token", lambda collection: "sig=x")

    def get(url, timeout):
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(sentinel2.provenance.requests, "get", get)
    return searches


AOI = gpd.GeoDataFrame(
    geometry=[box(X0, Y0 - 100, X0 + 100, Y0)], crs="EPSG:32605"
).to_crs(4326)
BANDS = ("blue", "swir16", "scl")


def test_load_mosaics_tiles_keeps_the_newest_processing_and_scales_by_baseline(pc):
    ds = sentinel2.load(AOI, "2023-09-01", "2023-09-30", BANDS)

    assert list(ds.data_vars) == ["blue", "swir16", "scl"]
    assert ds.rio.crs.to_epsg() == 32605
    assert abs(ds.rio.resolution()[0]) == 10  # 20 m bands upsampled to the finest
    # One step per pass; the older processing of a-south is dropped
    assert ds.time.dt.day.values.tolist() == [12, 19]
    assert ds.scene_ids.values.tolist() == ["a-north a-south", "b"]
    assert ds.mgrs_tiles.values.tolist() == ["05VNL 05VNK", "05VNL"]
    assert ds.cloud_cover.values.tolist() == [30.0, 10.0]
    assert ds.relative_orbit.values.tolist() == [129, 129]

    a, b = ds.isel(time=0).load(), ds.isel(time=1).load()
    # Baseline 05.10: (2000 - 1000) / 10000; 03.01 has no offset: 2000 / 10000
    np.testing.assert_allclose(np.nanmedian(a.blue), 0.1, atol=1e-6)
    np.testing.assert_allclose(np.nanmedian(a.swir16), 0.2, atol=1e-6)
    np.testing.assert_allclose(np.nanmedian(b.blue.values[:5]), 0.2, atol=1e-6)
    assert np.isnan(a.blue.values[0, 0])  # DN 0 is NaN
    assert a.blue.dtype == "float32" and a.blue.attrs["units"] == "1"
    # SCL: uint8 CF flags; cells outside any tile are no_data
    assert a.scl.dtype == "uint8" and a.scl.values[3, 3] == 4
    assert (b.scl.values[5:] == 0).all()
    assert a.scl.attrs["flag_meanings"].split()[11] == "snow_or_ice"
    assert ds.blue.encoding == {
        "dtype": "uint16",
        "scale_factor": 1e-4,
        "add_offset": -0.1,
        "_FillValue": 0,
    }
    hrefs = [
        lk["href"].rsplit("/", 1)[-1]
        for lk in ds.attrs["links"]
        if "/items/" in lk["href"]
    ]
    assert hrefs == ["a-north", "a-south", "b"]
    assert ds.attrs["license"] == "proprietary"


def test_cloud_cover_filters_the_search(pc):
    sentinel2.load(AOI, "2023-09-01", "2023-09-30", BANDS, max_cloud_cover=50)
    sentinel2.load(AOI, "2023-09-01", "2023-09-30", BANDS)
    assert [s["query"] for s in pc] == [{"eo:cloud_cover": {"lte": 50}}, None]
    with pytest.raises(ValueError, match="is a percentage"):
        sentinel2.load(AOI, "2023-09-01", "2023-09-30", BANDS, max_cloud_cover=150)


def test_windows_choose_passes(pc):
    windows = parse_window(
        {
            "pre": {"before": "30d", "select": "last"},
            "post": {"after": "30d", "select": "first"},
        },
        pd.Timestamp("2023-09-14T02:06:49"),
    )
    ds = sentinel2.load(AOI, "2023-08-15", "2023-10-14", BANDS, windows=windows)
    assert ds.window.values.tolist() == ["pre", "post"]
    with pytest.raises(ValueError, match="no Sentinel-2 L2A passes match"):
        sentinel2.load(
            AOI,
            "2023-08-15",
            "2023-10-14",
            BANDS,
            windows=parse_window("2023-10-01/2023-10-02", None),
        )


def test_no_scenes_is_an_error(pc, tiles):
    tiles.clear()
    with pytest.raises(ValueError, match="with at most 5% cloud"):
        sentinel2.load(AOI, "2023-09-01", "2023-09-30", BANDS, max_cloud_cover=5)


def test_pass_key_ignores_the_processing_baseline():
    a = pystac.Item(
        "a",
        None,
        None,
        dt.datetime(2023, 9, 12),
        {"s2:datatake_id": "GS2B_20230912T214529_034046_N05.09"},
    )
    b = pystac.Item(
        "b",
        None,
        None,
        dt.datetime(2023, 9, 12),
        {"s2:datatake_id": "GS2B_20230912T214529_034046_N05.10"},
    )
    assert (
        sentinel2.pass_key(a) == sentinel2.pass_key(b) == "GS2B_20230912T214529_034046"
    )


def test_reflectance_is_stored_packed_and_reads_back(pc, tmp_path):
    ds = sentinel2.load(AOI, "2023-09-01", "2023-09-30", BANDS)
    path = str(tmp_path / "s2.zarr")
    io.save_datatree(xr.DataTree.from_dict({"rasters/sentinel2": ds}), path)

    group = zarr.open_group(path, mode="r")["rasters/sentinel2"]
    assert group["blue"].dtype == "uint16" and group["scl"].dtype == "uint8"
    back = io.open_datatree(path)["rasters/sentinel2"].to_dataset()
    np.testing.assert_allclose(back.blue, ds.blue, atol=1e-5)
    assert back.scl.dtype == "uint8"


def test_cli_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"bbox": [-151.3, 63.1, -151.1, 63.2]}).to_file(aoi)
    stage_app(
        [
            "sentinel2",
            "-i",
            str(aoi),
            "-s",
            "2023-09-01",
            "-e",
            "2023-09-30",
            "-o",
            str(tmp_path / "o"),
        ]
        + ["-b", "red", "green", "blue", "-m", "50", "--dry-run"],
        result_action="return_value",
    )
    out = capsys.readouterr().out
    assert "bands: [red, green, blue]" in out and "max_cloud_cover: 50" in out
