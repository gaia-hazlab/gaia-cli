import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import dem, runner, stagefile
from gaia_cli.aoi import load_aoi
from gaia_cli.cli import stage_app
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


AOI = "aoi: {bbox: [-121.6, 47.6, -121.0, 47.9]}\n"


def config(tmp_path, datasets):
    path = tmp_path / "stage.yaml"
    path.write_text(f"name: x\n{AOI}datasets:\n{datasets}")
    return str(path)


def test_stage_file_keys_name_the_product_and_group(tmp_path):
    cfg = stagefile.load_config(config(tmp_path, "  dem: {}\n  dem/usgs10m: {}\n"))
    names = [r.name for r in cfg.requests]
    assert names == ["dem/cop30", "dem/usgs10m"]  # a bare dem is cop30
    assert [r.options for r in cfg.requests] == [
        {"product": Product.cop30},
        {"product": Product.usgs10m},
    ]
    assert [stagefile.group(r) for r in cfg.requests] == names
    assert all(r.windows == [] for r in cfg.requests)


@pytest.mark.parametrize(
    "datasets, message",
    [
        ("  dem/srtm: {}\n", "dem/srtm: 'srtm' is not one of cop30, usgs10m"),
        ("  dem: {}\n  dem/cop30: {}\n", "dem/cop30: dem/cop30 is listed twice"),
        ("  dem: {product: usgs10m}\n", "dem: unknown option 'product'"),
        ("  nlcd/2021: {}\n", "nlcd/2021: nlcd has no variants; use nlcd"),
    ],
)
def test_bad_variant_keys(tmp_path, datasets, message):
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(config(tmp_path, datasets))
    assert message in e.value.problems


def fake_load(calls):
    def load(aoi, product: Product = Product.cop30):
        calls.append(product)
        return xr.Dataset(
            {"elevation": (("y", "x"), [[1.0]])}, attrs={"title": product.value}
        )

    return load


def test_both_products_go_in_one_store(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(dem, "load", fake_load(calls))
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"bbox": [-121.6, 47.6, -121.0, 47.9]}).to_file(aoi)
    out = tmp_path / "out"
    stage_app(["dem", "-i", str(aoi), "-o", str(out)], result_action="return_value")
    stage_app(
        ["dem", "-i", str(aoi), "-o", str(out), "-p", "usgs10m", "--skip-existing"],
        result_action="return_value",
    )

    assert calls == [Product.cop30, Product.usgs10m]
    tree = xr.open_datatree(out / "out.zarr", engine="zarr")
    assert set(tree["dem"].children) == {"cop30", "usgs10m"}
    assert tree["dem/usgs10m"].attrs["title"] == "usgs10m"
    # gaia align's default for the DEM, stamped by the runner
    assert tree["dem/usgs10m"].elevation.attrs["gaia:regrid-method"] == "bilinear"
    # The stage file repeats both, keyed by product, without a product option
    doc = stagefile.read_config(f"{out}/stage.yaml")
    assert list(doc["datasets"]) == ["dem/cop30", "dem/usgs10m"]
    assert "product" not in doc["datasets"]["dem/usgs10m"]


def test_only_dem_selects_every_product(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(dem, "load", fake_load(calls))
    path = config(tmp_path, "  dem/cop30: {}\n  dem/usgs10m: {}\n")
    runner.stage(config=path, output_path=str(tmp_path / "a"), only="dem")
    runner.stage(config=path, output_path=str(tmp_path / "b"), only="dem/usgs10m")
    assert calls == [Product.cop30, Product.usgs10m, Product.usgs10m]
