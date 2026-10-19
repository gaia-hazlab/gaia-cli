import json

import geopandas as gpd
import numpy as np
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box

from gaia_cli import burn_severity, io, stagefile

# A 10 x 10 fire raster at 20 m in UTM 10N, like the service's per-fire rasters
X0, Y0, RES = 680000.0, 5350000.0, 20.0
EXTENT = {
    "xmin": X0,
    "ymin": Y0 - 10 * RES,
    "xmax": X0 + 10 * RES,
    "ymax": Y0,
    "spatialReference": {"wkid": 26910, "latestWkid": 26910},
}
FIRES = [
    {"objectid": 7, "name": "bolt_creek_sbs", "year": 2022},
    {"objectid": 8, "name": "irving_peak_wa4789712102720220811_sbs", "year": 2022},
    {"objectid": 9, "name": "Bolt_Creek_sbs2", "year": 2017},
]


def footprint(xmin, ymin, xmax, ymax, crs="EPSG:26910"):
    return gpd.GeoSeries([box(xmin, ymin, xmax, ymax)], crs=crs).to_crs(4326).iloc[0]


class Response:
    def __init__(self, data=None, content=b"", content_type="application/json"):
        self.data, self.content = data, content
        self.headers = {"content-type": content_type}

    def raise_for_status(self):
        pass

    def json(self):
        return self.data


@pytest.fixture
def service(tmp_path, monkeypatch):
    """The image service, answering from one small fire raster.

    Classes 1-4 by row, the service's 15 in the last row and a source 0 in the
    first cell. Returns the list of requests made.
    """
    values = np.repeat(np.array([1, 2, 3, 4, 1, 2, 3, 4, 1, 15]), 10).reshape(10, 10)
    values[0, 0] = 0
    da = xr.DataArray(
        values.astype("uint8"),
        dims=("y", "x"),
        coords={
            "y": Y0 - RES / 2 - RES * np.arange(10),
            "x": X0 + RES / 2 + RES * np.arange(10),
        },
    ).rio.write_crs("EPSG:26910")
    tif = tmp_path / "fire.tif"
    da.rio.to_raster(tif)
    far = footprint(X0 + 50000, Y0, X0 + 51000, Y0 + 1000)  # in the bbox, not the AOI
    features = [
        (
            7,
            "bolt_creek_sbs",
            footprint(*[EXTENT[k] for k in ("xmin", "ymin", "xmax", "ymax")]),
        ),
        (8, "irving_peak_sbs", far),
    ]
    calls = []

    def get(url, params, timeout):
        calls.append((url, params))
        if url.endswith("/query"):
            return Response(
                {
                    "features": [
                        {
                            "attributes": {"objectid": i, "name": n, "beginyear": 2022},
                            "geometry": {"rings": [list(g.exterior.coords)]},
                        }
                        for i, n, g in features
                    ]
                }
            )
        if url.endswith("/info"):
            return Response({"extent": EXTENT, "pixelSizeX": RES, "pixelSizeY": RES})
        assert url.endswith("/exportImage")
        return Response(content=tif.read_bytes(), content_type="image/tiff")

    monkeypatch.setattr(burn_severity.requests, "get", get)
    return calls


def whole_fire():
    return gpd.GeoDataFrame(
        geometry=[box(X0, Y0 - 10 * RES, X0 + 10 * RES, Y0)], crs="EPSG:26910"
    )


def test_select_fire_by_name_and_year():
    select = burn_severity.select_fire
    assert select(FIRES, "irving")["objectid"] == 8
    assert select(FIRES, "Bolt Creek", 2022)["objectid"] == 7  # case, spaces
    assert select(FIRES, year=2017)["objectid"] == 9
    assert select(FIRES[:1])["objectid"] == 7  # one fire needs no choice


def test_select_fire_errors_list_the_fires():
    with pytest.raises(ValueError, match="No BAER soil burn severity map"):
        burn_severity.select_fire([])
    with pytest.raises(ValueError, match=r"3 fires .* Choose one"):
        burn_severity.select_fire(FIRES)
    with pytest.raises(ValueError, match=r"2 fires .*bolt_creek_sbs \(2022\)"):
        burn_severity.select_fire(FIRES, "bolt")
    with pytest.raises(ValueError, match=r"it has bolt_creek_sbs \(2022\)"):
        burn_severity.select_fire(FIRES, "dixie")


def test_find_fires_keeps_footprints_meeting_the_aoi(service):
    fires = burn_severity.find_fires(whole_fire())
    assert fires == [{"objectid": 7, "name": "bolt_creek_sbs", "year": 2022}]
    _, params = service[0]
    assert params["where"] == "category = 1"


def test_window_snaps_to_the_raster_grid():
    info = {"extent": EXTENT, "pixelSizeX": RES, "pixelSizeY": RES}
    aoi = gpd.GeoDataFrame(
        geometry=[box(X0 + 25, Y0 - 75, X0 + 55, Y0 - 45)], crs="EPSG:26910"
    )
    bbox, width, height = burn_severity._window(info, aoi)
    assert bbox == [X0 + 20, Y0 - 80, X0 + 60, Y0 - 40]
    assert (width, height) == (2, 2)
    # Beyond the raster: cut at its edges, or nothing at all
    big = gpd.GeoDataFrame(
        geometry=[box(X0 - 1e4, Y0 - 1e4, X0 + 1e4, Y0 + 1e4)], crs="EPSG:26910"
    )
    assert burn_severity._window(info, big)[1:] == (10, 10)
    off = gpd.GeoDataFrame(
        geometry=[box(X0 + 1e4, Y0, X0 + 2e4, Y0 + 1e4)], crs="EPSG:26910"
    )
    assert 0 in burn_severity._window(info, off)[1:]


def test_load_keeps_classes_and_masks_the_rest(service):
    ds = burn_severity.load(whole_fire())

    da = ds.burn_severity
    assert da.rio.crs.to_epsg() == 26910  # native grid, not reprojected
    assert da.shape == (10, 10)
    assert np.isnan(da.values[0, 0])  # a source 0
    assert np.isnan(da.values[9]).all()  # the service's 15
    assert da.values[1:9, 0].tolist() == [2, 3, 4, 1, 2, 3, 4, 1]
    assert da.attrs["flag_values"] == [1, 2, 3, 4]
    assert da.attrs["flag_meanings"] == "unburned_to_very_low low moderate high"
    assert ds.attrs["fire_name"] == "bolt_creek_sbs"
    assert ds.attrs["start_datetime"].startswith("2022-01-01")

    url, params = service[-1]
    assert url.endswith("/exportImage")
    assert params["size"] == "10,10"
    assert params["bboxSR"] == params["imageSR"] == 26910
    assert json.loads(params["mosaicRule"])["lockRasterIds"] == [7]


def test_store_keeps_the_source_bytes(service, tmp_path):
    ds = burn_severity.load(whole_fire())
    path = tmp_path / "t.zarr"
    io.save_datatree(xr.DataTree.from_dict({"rasters/burn_severity": ds}), str(path))

    raw = xr.open_zarr(path, group="rasters/burn_severity", mask_and_scale=False)
    assert raw.burn_severity.dtype == "uint8"
    assert int(raw.burn_severity.values[9, 0]) == burn_severity.FILL
    back = xr.open_zarr(path, group="rasters/burn_severity").load()
    np.testing.assert_array_equal(back.burn_severity.values, ds.burn_severity.values)


AOI = "aoi: {bbox: [-121.6, 47.6, -121.0, 47.9]}\n"


def config(tmp_path, datasets):
    path = tmp_path / "stage.yaml"
    path.write_text(f"name: x\n{AOI}datasets:\n{datasets}")
    return str(path)


def test_stage_file_options_and_group(tmp_path):
    cfg = stagefile.load_config(
        config(tmp_path, "  burn_severity: {fire: bolt creek, year: 2022}\n")
    )
    (req,) = cfg.requests
    assert req.options == {"fire": "bolt creek", "year": 2022}
    assert stagefile.group(req) == "rasters/burn_severity"
    assert req.windows == []


def test_stage_file_rejects_a_window(tmp_path):
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(
            config(tmp_path, "  burn_severity: {window: 2022-01-01/2022-12-31}\n")
        )
    assert "burn_severity: is static and takes no window" in e.value.problems
