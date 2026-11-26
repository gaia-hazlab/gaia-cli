import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box, mapping

from gaia_cli import earthdata, fetch, io, regrid
from gaia_cli.cli import stage_app
from gaia_cli.datasets import nisar_sme2
from gaia_cli.datasets.sentinel1_rtc import OrbitState
from gaia_cli.window import parse_window

# 10 x 10 cells of the 200 m EASE-Grid 2.0 near Odessa, WA: whole cells from
# the grid's origin, as the files' are. Frames 25 and 26 of a pass overlap in
# rows 4-5
CELL = 200.17900466918945
X0, Y0 = -58635 * CELL, 28080 * CELL
BOUNDS = (X0, Y0 - 10 * CELL, X0 + 10 * CELL, Y0)
SOUTH = (X0, Y0 - 10 * CELL, X0 + 10 * CELL, Y0 - 4 * CELL)
NORTH = (X0, Y0 - 6 * CELL, X0 + 10 * CELL, Y0)
COLLECTION = {
    "title": "NISAR Provisional Soil Moisture (Version 1)",
    "license": "CC0-1.0",
    "links": [{"rel": "license", "href": "https://www.earthdata.nasa.gov/"}],
}
AOI = gpd.GeoDataFrame(geometry=[box(*BOUNDS)], crs="EPSG:6933").to_crs(4326)


def tif(path, bounds, value, dtype, nodata):
    """A GeoTIFF of *value* over *bounds* on the 200 m grid, its first cell
    *nodata*."""
    left, bottom, right, top = bounds
    ny, nx = round((top - bottom) / CELL), round((right - left) / CELL)
    data = np.full((ny, nx), value, dtype=dtype)
    data[0, 0] = nodata
    da = xr.DataArray(
        data,
        dims=("y", "x"),
        coords={
            "y": top - CELL * (np.arange(ny) + 0.5),
            "x": left + CELL * (np.arange(nx) + 0.5),
        },
    )
    da.rio.write_crs("EPSG:6933").rio.write_nodata(nodata).rio.to_raster(path)


def feature(root, frame, bounds, value, when, orbit, track=77, collection=nisar_sme2.COLLECTION, processed="2026-10-19T00:00:00Z"):  # fmt: skip
    """An asf_search GeoJSON feature for one SME2 frame; its grids are
    GeoTIFFs under *root*, named as :func:`local` finds them."""
    t = pd.Timestamp(when)
    name = (
        f"NISAR_L3_PR_SME2_033_{track:03d}_A_{frame:03d}_4005_DHDH_A_"
        f"{t:%Y%m%dT%H%M%S}_{t:%Y%m%dT%H%M%S}_P05023_N_F_J_"
        f"{pd.Timestamp(processed):%j}"
    )
    for v, band in nisar_sme2.BANDS.items():
        tif(root / f"{name}_{v}.tif", bounds, value, band["data_type"], band["nodata"])
    footprint = gpd.GeoSeries([box(*bounds)], crs="EPSG:6933").to_crs(4326)[0]
    return {
        "type": "Feature",
        "geometry": mapping(footprint),
        "properties": {
            "sceneName": name,
            "startTime": when,
            "processingDate": processed,
            "platform": "NISAR",
            "flightDirection": "ASCENDING",
            "pathNumber": track,
            "frameNumber": frame,
            "orbit": orbit,
            "collectionName": collection,
            "url": str(root / f"{name}.h5"),
        },
    }


def local(href, name):
    """The GeoTIFF standing in for *name* of the HDF5 file at *href*."""
    return href.removesuffix(".h5") + f"_{name}.tif"


@pytest.fixture
def frames(tmp_path):
    """Two overlapping frames of one pass, a reprocessing of the second, a
    later pass on another track, and a GCOV granule."""
    return [
        feature(tmp_path, 25, SOUTH, 1, "2026-10-17T12:46:48Z", 6404),
        feature(tmp_path, 26, NORTH, 2, "2026-10-17T12:47:22Z", 6404),
        feature(
            tmp_path,
            26,
            NORTH,
            3,
            "2026-10-17T12:47:22Z",
            6404,
            processed="2026-11-10T00:00:00Z",
        ),
        feature(tmp_path, 92, BOUNDS, 4, "2026-10-21T20:37:40Z", 6466, track=139),
        feature(
            tmp_path,
            26,
            NORTH,
            9,
            "2026-10-17T12:47:22Z",
            6404,
            collection="NISAR_L2_GCOV_PROVISIONAL_V1",
        ),
    ]


@pytest.fixture
def asf(monkeypatch, frames):
    """ASF stand-ins, the files local; returns the searches made."""
    searches = []

    def search(aoi, start=None, end=None, **filters):
        searches.append({"start": start, "end": end, **filters})
        return frames

    monkeypatch.setattr(earthdata, "search", search)
    monkeypatch.setattr(earthdata, "sign", lambda items, assets: None)
    monkeypatch.setattr(nisar_sme2, "subdataset", local)

    def get(url, **kwargs):
        assert url.startswith(nisar_sme2.STAC)
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(fetch, "get", get)
    return searches


def test_load_keeps_each_frame_apart(asf, caplog):
    caplog.set_level(logging.INFO)
    ds = nisar_sme2.load(AOI, "2026-10-01T00:00:00", "2026-10-31T23:59:59")
    # The older processing and the GCOV granule are skipped
    assert "skipped 2 older or other processing(s)" in caplog.text
    assert [str(t)[:19] for t in ds.time.values] == [
        "2026-10-17T12:46:48",
        "2026-10-17T12:47:22",
        "2026-10-21T20:37:40",
    ]
    assert ds.frame.values.tolist() == [25, 26, 92]
    assert ds.relative_orbit.values.tolist() == [77, 77, 139]
    assert ds.orbit_state.values.tolist() == ["ascending"] * 3
    assert set(ds.data_vars) == set(nisar_sme2.VARIABLES)
    assert ds.rio.crs.to_epsg() == 6933
    assert ds.rio.resolution() == pytest.approx((CELL, -CELL))
    sm = ds.soilMoisture
    # Overlapping rows 4-5: each frame its own value, not a mosaic
    row = Y0 - 4.5 * CELL
    assert sm.isel(time=0).sel(y=row, x=X0 + 5.5 * CELL, method="nearest") == 1
    assert sm.isel(time=1).sel(y=row, x=X0 + 5.5 * CELL, method="nearest") == 3
    # Outside its frame, and its -9999 cell: NaN
    assert np.isnan(sm.isel(time=0).sel(y=Y0 - 0.5 * CELL, x=X0 + 5.5 * CELL, method="nearest"))  # fmt: skip
    assert int(sm.isel(time=1).isnull().sum()) == 4 * 10 + 1
    assert sm.attrs["units"] == "m3 m-3"
    assert sm.attrs["standard_name"] == "volume_fraction_of_condensed_water_in_soil"
    # The flags keep their dtype and fill value
    flag = ds.retrievalQualityFlag
    assert flag.dtype == np.int16
    assert flag.encoding["_FillValue"] == -9999
    assert flag.attrs["flag_masks"] == [1, 2, 4, 8]
    assert flag.attrs[regrid.METHOD_ATTR] == "nearest"
    assert int(flag.isel(time=1).sel(y=Y0 - 0.5 * CELL, x=X0 + 0.5 * CELL, method="nearest")) == -9999  # fmt: skip
    assert ds.attrs["version"] == nisar_sme2.COLLECTION
    items = [x["href"] for x in ds.attrs["links"] if "/items/" in x["href"]]
    assert len(items) == 3
    assert all(f"{nisar_sme2.COLLECTION}_1/items/" in h for h in items)


def test_options_filter_the_search(asf):
    nisar_sme2.load(AOI, "2026-10-01", "2026-10-31", ("soilMoisture",), OrbitState.ascending, 77)  # fmt: skip
    assert asf[0]["flightDirection"] == "ASCENDING"
    assert asf[0]["relativeOrbit"] == 77
    assert asf[0]["dataset"] == "NISAR" and asf[0]["processingLevel"] == "SME2"


def test_ancillary_variables(asf):
    variables = ("soilMoisture", "landCover", "waterbodyFraction")
    ds = nisar_sme2.load(AOI, "2026-10-01", "2026-10-31", variables)
    assert set(ds.data_vars) == set(variables)
    assert ds.landCover.dtype == np.int8
    assert ds.landCover.encoding["_FillValue"] == -127
    assert ds.waterbodyFraction.attrs["units"] == "1"


def test_no_frames_is_an_error(asf, frames):
    del frames[:4]  # only the GCOV granule
    with pytest.raises(ValueError, match="no NISAR SME2 frames"):
        nisar_sme2.load(AOI, "2026-10-01", "2026-10-31")


def test_one_header_per_frame(asf, monkeypatch):
    keys = []
    set_proj = earthdata.set_proj

    def spy(items, asset, key):
        keys.append({key(it) for it in items})
        set_proj(items, asset, key)

    monkeypatch.setattr(earthdata, "set_proj", spy)
    nisar_sme2.load(AOI, "2026-10-01", "2026-10-31")
    assert keys == [{(77, 25), (77, 26), (139, 92)}]


def test_load_keeps_the_frames_the_windows_choose(asf):
    windows = parse_window(
        {"pre": {"before": "3d", "select": "last"}, "post": {"after": "3d"}},
        pd.Timestamp("2026-10-20"),
    )
    ds = nisar_sme2.load(AOI, "2026-10-17", "2026-10-23", windows=windows)
    assert ds.frame.values.tolist() == [26, 92]
    assert ds.window.values.tolist() == ["pre", "post"]


def test_subdataset_paths():
    assert nisar_sme2.subdataset("https://x/f.h5?sig=1", "soilMoisture") == (
        'HDF5:"/vsicurl/https://x/f.h5?sig=1"://science/LSAR/SME2/grids/soilMoisture'
    )
    assert nisar_sme2.subdataset("f.h5", "landCover").endswith(
        "//science/LSAR/SME2/grids/ancillaryData/landCover"
    )


def test_cli_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["nisar-sme2", "-i", str(aoi), "-s", "2026-10-01", "-e", "2026-10-31"]
        + ["-o", str(tmp_path / "o"), "-v", "soilMoisture", "--dry-run"],
        result_action="return_value",
    )
    assert "variables: [soilMoisture]" in capsys.readouterr().out


def test_the_runner_reads_it_in_processes(asf, tmp_path, monkeypatch):
    calls = []

    def compute(tree, env):
        calls.append(env)
        return tree

    monkeypatch.setattr(io, "compute_in_processes", compute)
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["nisar-sme2", "-i", str(aoi), "-s", "2026-10-01", "-e", "2026-10-31"]
        + ["-o", str(tmp_path / "o")],
        result_action="return_value",
    )
    assert calls == [nisar_sme2.PROCESS_ENV]
    assert calls[0]["CPL_VSIL_CURL_CHUNK_SIZE"] == "2097152"
