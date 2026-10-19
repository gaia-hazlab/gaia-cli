import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import io, soil, stagefile
from gaia_cli.soil import Product

# Near Bolt Creek, WA, in the SOLUS grid's CRS (CONUS Albers)
X0, Y0 = -1890000.0, 3000000.0


@pytest.fixture
def local_solus(tmp_path, monkeypatch):
    """Small 100 m GeoTIFFs with SOLUS's dtypes and nodata, in place of the COGs.

    Every layer holds 0..99 times its scalar, with nodata in the first cell.
    """
    paths = {}
    for layer, spec in soil.SOLUS_LAYERS.items():
        dtype, nodata = ("uint8", 255) if spec["scalar"] == 1 else ("uint16", 65535)
        values = (np.arange(100).reshape(10, 10) * spec["scalar"]).astype(dtype)
        values[0, 0] = nodata
        da = xr.DataArray(
            values,
            dims=("y", "x"),
            coords={
                "y": Y0 - 50 - 100 * np.arange(10),
                "x": X0 + 50 + 100 * np.arange(10),
            },
        )
        da = da.rio.write_crs("EPSG:5070").rio.write_nodata(nodata)
        paths[layer] = tmp_path / f"{layer}.tif"
        da.rio.to_raster(paths[layer])
    monkeypatch.setattr(soil, "_url", lambda layer: str(paths[layer]))


def whole_fixture():
    return gpd.GeoDataFrame(
        geometry=[box(X0, Y0 - 1000, X0 + 1000, Y0)], crs="EPSG:5070"
    )


def test_load_decodes_scaled_integers(local_solus):
    ds = soil.load(whole_fixture())

    assert set(ds.data_vars) == set(soil.SOLUS_LAYERS)
    assert ds.rio.crs.to_epsg() == 5070  # native grid, not reprojected
    assert (ds.sizes["y"], ds.sizes["x"]) == (10, 10)
    for name, da in ds.data_vars.items():
        assert np.isnan(da.values[0, 0])  # nodata is masked
        assert da.values[9, 9] == 99  # stored 99 * scalar
        assert da.attrs["units"] == soil.SOLUS_LAYERS[name]["units"]
    assert ds.ph1to1h2o.attrs["depth"] == "0 cm"
    assert ds.attrs["license"] == "CC-BY-4.0"
    rels = {lk["rel"] for lk in ds.attrs["links"]}
    assert {"via", "license", "cite-as", "describedby"} <= rels


def test_store_keeps_the_source_integers(local_solus, tmp_path):
    ds = soil.load(whole_fixture())
    path = tmp_path / "t.zarr"
    io.save_datatree(xr.DataTree.from_dict({"soil/solus100": ds}), str(path))

    raw = xr.open_zarr(path, group="soil/solus100", mask_and_scale=False)
    assert raw.claytotal.dtype == "uint8"
    assert raw.ph1to1h2o.dtype == "uint16"
    assert int(raw.ph1to1h2o.values[9, 9]) == 9900
    assert int(raw.ph1to1h2o.values[0, 0]) == 65535
    # Readers decode to the same values load() returned
    back = xr.open_zarr(path, group="soil/solus100").load()
    for name in soil.SOLUS_LAYERS:
        np.testing.assert_array_equal(back[name].values, ds[name].values)


AOI = "aoi: {bbox: [-121.6, 47.6, -121.0, 47.9]}\n"


def config(tmp_path, datasets):
    path = tmp_path / "stage.yaml"
    path.write_text(f"name: x\n{AOI}datasets:\n{datasets}")
    return str(path)


def test_stage_file_key_names_the_product_and_group(tmp_path):
    cfg = stagefile.load_config(config(tmp_path, "  soil: {}\n"))
    (req,) = cfg.requests
    assert req.name == "soil/solus100"  # a bare soil is solus100
    assert req.options == {"product": Product.solus100}
    assert stagefile.group(req) == "soil/solus100"
    assert req.windows == []


def test_unknown_product(tmp_path):
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(config(tmp_path, "  soil/soilgrids: {}\n"))
    assert "soil/soilgrids: 'soilgrids' is not one of solus100" in e.value.problems


def test_extent_is_conus_only():
    def outside(bbox):
        req = stagefile.Request("soil/solus100", {"bbox": bbox}, [], {})
        return stagefile.outside_extent(req, stagefile.load_aoi({"bbox": bbox}))

    assert outside([-121.6, 47.6, -121.0, 47.9]) is None  # Bolt Creek
    assert "outside the dataset's extent" in outside([-152.2, 62.75, -150.2, 63.55])
