import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import dem, stagefile
from gaia_cli.dem import Product


@pytest.fixture
def local_usgs10m(tmp_path, monkeypatch):
    """A small NAD83 GeoTIFF with the usgs10m nodata, in place of the VRT."""
    values = np.arange(100, dtype="float32").reshape(10, 10)
    values[0, 0] = -999999
    da = xr.DataArray(
        values,
        dims=("y", "x"),
        coords={
            "y": np.linspace(47.95, 47.05, 10),
            "x": np.linspace(-121.95, -121.05, 10),
        },
    )
    da = da.rio.write_crs("EPSG:4269").rio.write_nodata(-999999)
    path = tmp_path / "dem.tif"
    da.rio.to_raster(path)
    monkeypatch.setitem(dem.PRODUCTS[Product.usgs10m], "url", str(path))


def test_load_masks_nodata_and_records_the_vertical_datum(local_usgs10m):
    aoi = gpd.GeoDataFrame(geometry=[box(-122.0, 47.5, -121.5, 48.0)], crs="EPSG:4326")
    ds = dem.load(aoi, Product.usgs10m)

    da = ds.elevation
    assert da.rio.crs.to_epsg() == 4269  # native grid, not reprojected
    assert (da.sizes["y"], da.sizes["x"]) == (5, 5)
    assert np.isnan(da.values[0, 0])  # -999999 is masked
    assert da.attrs["units"] == "m"
    assert da.attrs["vertical_datum"] == "EPSG:5703"
    assert ds.attrs["sci:doi"] == "10.5069/G98K778D"
    rels = {lk["rel"] for lk in ds.attrs["links"]}
    assert {"via", "license", "cite-as"} <= rels


@pytest.mark.parametrize("product", list(Product))
def test_every_product_has_provenance(product):
    attrs = dem.source_attrs(product)
    assert attrs["title"] and attrs["license"] and attrs["providers"]
    assert dem.PRODUCTS[product]["vertical_datum"].startswith("EPSG:")


def test_stage_file_product_option(tmp_path):
    path = tmp_path / "stage.yaml"
    path.write_text(
        "name: x\naoi: {bbox: [-121.6, 47.6, -121.0, 47.9]}\n"
        "datasets: {dem: {product: usgs10m}}\n"
    )
    (req,) = stagefile.load_config(str(path)).requests
    assert req.options == {"product": Product.usgs10m}
    assert req.windows == []

    path.write_text(path.read_text().replace("usgs10m", "srtm"))
    with pytest.raises(stagefile.StageFileError, match="not one of cop30, usgs10m"):
        stagefile.load_config(str(path))
