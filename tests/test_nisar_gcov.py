import logging

import cf_xarray
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box, mapping

from gaia_cli import earthdata, fetch
from gaia_cli.cli import stage_app
from gaia_cli.datasets import cf_criteria, nisar_gcov, opera_rtc, sentinel1_rtc
from gaia_cli.datasets.nisar_gcov import Maturity
from gaia_cli.datasets.sentinel1_rtc import OrbitState
from gaia_cli.window import parse_window

# 10 x 10 cells of 10 m in UTM 10N, near Paradise; each frame is half of it
X0, Y0 = 595000, 5182000
NORTH = (X0, Y0 - 50, X0 + 100, Y0)
SOUTH = (X0, Y0 - 100, X0 + 100, Y0 - 50)
VALIDATED, PROVISIONAL = nisar_gcov.COLLECTIONS.values()
COLLECTION = {
    "title": "NISAR GCOV",
    "license": "CC0-1.0",
    "links": [{"rel": "license", "href": "https://www.earthdata.nasa.gov/"}],
}
AOI = gpd.GeoDataFrame(
    geometry=[box(X0, Y0 - 100, X0 + 100, Y0)], crs="EPSG:32610"
).to_crs(4326)


def tif(path, bounds, value, dtype="float32", nodata=np.nan):
    """A GeoTIFF of *value* over *bounds* (5 x 10 cells of 10 m)."""
    left, bottom, right, top = bounds
    da = xr.DataArray(
        np.full((5, 10), value, dtype=dtype),
        dims=("y", "x"),
        coords={
            "y": np.arange(top - 5, bottom, -10.0),
            "x": np.arange(left + 5, right, 10.0),
        },
    )
    da.rio.write_crs("EPSG:32610").rio.write_nodata(nodata).rio.to_raster(path)


def feature(root, frame, bounds, value, when, orbit, track=172, collection=VALIDATED, processed="2026-11-01T00:00:00Z", pols=("HH", "HV"), mask=1):  # fmt: skip
    """An asf_search GeoJSON feature for one GCOV frame; its grids are GeoTIFFs
    under *root*, named as :func:`local` finds them."""
    t = pd.Timestamp(when)
    crid = "P05036" if collection == VALIDATED else "P05023"
    name = (
        f"NISAR_L2_PR_GCOV_030_{track:03d}_D_{frame:03d}_4005_DHDH_A_"
        f"{t:%Y%m%dT%H%M%S}_{t:%Y%m%dT%H%M%S}_{crid}_N_F_J_"
        f"{pd.Timestamp(processed):%j}"
    )
    for pol in pols:
        tif(root / f"{name}_{pol}{pol}.tif", bounds, value)
    tif(root / f"{name}_mask.tif", bounds, mask, "uint8", 255)
    for layer in ("numberOfLooks", "rtcGammaToSigmaFactor"):
        tif(root / f"{name}_{layer}.tif", bounds, 0.5)
    footprint = gpd.GeoSeries([box(*bounds)], crs="EPSG:32610").to_crs(4326)[0]
    return {
        "type": "Feature",
        "geometry": mapping(footprint),
        "properties": {
            "sceneName": name,
            "startTime": when,
            "processingDate": processed,
            "platform": "NISAR",
            "flightDirection": "DESCENDING",
            "pathNumber": track,
            "frameNumber": frame,
            "orbit": orbit,
            "rangeBandwidth": ["40+5"],
            "mainBandPolarization": list(pols),
            "collectionName": collection,
            "url": str(root / f"{name}.h5"),
        },
    }


def local(href, name):
    """The GeoTIFF standing in for *name* of the HDF5 file at *href*."""
    return href.removesuffix(".h5") + f"_{name}.tif"


@pytest.fixture
def frames(tmp_path):
    """Two frames of one pass, a reprocessed frame, an HH-only pass, a later
    pass, and a provisional processing of the first pass."""
    return [
        feature(tmp_path, 64, NORTH, 1.0, "2026-10-24T03:18:47Z", 6500, mask=0),
        feature(tmp_path, 65, SOUTH, 2.0, "2026-10-24T03:19:20Z", 6500),
        feature(
            tmp_path,
            65,
            SOUTH,
            9.0,
            "2026-10-24T03:19:20Z",
            6500,
            processed="2026-11-05T00:00:00Z",
        ),
        feature(
            tmp_path, 64, NORTH, 5.0, "2026-10-29T12:55:07Z", 6560, 5, pols=("HH",)
        ),
        feature(tmp_path, 64, NORTH, 3.0, "2026-11-05T03:18:47Z", 6670),
        feature(
            tmp_path,
            64,
            NORTH,
            7.0,
            "2026-10-24T03:18:47Z",
            6500,
            collection=PROVISIONAL,
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
    monkeypatch.setattr(nisar_gcov, "subdataset", local)

    def get(url, **kwargs):
        assert url.startswith(nisar_gcov.STAC)
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(fetch, "get", get)
    return searches


def test_load_mosaics_the_frames_of_each_pass(asf, caplog):
    caplog.set_level(logging.INFO)
    ds = nisar_gcov.load(AOI, "2026-10-21T00:00:00", "2026-11-08T23:59:59")
    # The HH-only pass is skipped; the older and provisional processings too
    assert [str(t)[:10] for t in ds.time.values] == ["2026-10-24", "2026-11-05"]
    assert "skipped 2 older or other processing(s)" in caplog.text
    assert "skipped 1 frame(s) without HHHH+HVHV" in caplog.text
    assert set(ds.data_vars) == {"HHHH", "HVHV", "mask"}
    assert ds.rio.resolution() == (10.0, -10.0)
    first = ds.HHHH.isel(time=0)
    # North frame 1.0, south frame's newest processing 9.0
    assert float(first.sel(y=Y0 - 5, x=X0 + 5)) == 1.0
    assert float(first.sel(y=Y0 - 95, x=X0 + 5)) == 9.0
    assert ds.HHHH.attrs["units"] == "1"
    assert ds.frames.values.tolist() == ["64 65", "64"]
    assert ds.relative_orbit.values.tolist() == [172, 172]
    assert ds.orbit_state.values.tolist() == ["descending", "descending"]
    assert ds.range_bandwidth.values.tolist() == ["40+5", "40+5"]
    # The mask keeps its 0 and its flags, not a _FillValue
    mask = ds["mask"].isel(time=0)
    assert mask.dtype == np.uint8
    assert set(np.unique(mask.values).tolist()) == {0, 1}
    assert ds["mask"].attrs["flag_values"] == [0, 255]
    assert ds.attrs["version"] == VALIDATED
    items = [x["href"] for x in ds.attrs["links"] if "/items/" in x["href"]]
    assert len(items) == 3 and all(f"{VALIDATED}_1/items/" in h for h in items)


def test_maturity_is_validated_unless_asked_and_never_mixed(asf):
    nisar_gcov.load(AOI, "2026-10-21", "2026-11-08")
    ds = nisar_gcov.load(AOI, "2026-10-21", "2026-11-08", maturity=Maturity.provisional)
    assert [s["dataMaturity"] for s in asf] == ["VALIDATED", "PROVISIONAL"]
    # Only the provisional processing, though the search returned both
    assert [str(t)[:10] for t in ds.time.values] == ["2026-10-24"]
    assert float(ds.HHHH.isel(time=0).sel(y=Y0 - 5, x=X0 + 5)) == 7.0
    assert ds.attrs["version"] == PROVISIONAL


def test_orbit_options_filter_the_search(asf):
    nisar_gcov.load(AOI, "2026-10-21", "2026-11-08", ("HHHH",), Maturity.validated, OrbitState.descending, 172)  # fmt: skip
    assert asf[0]["flightDirection"] == "DESCENDING"
    assert asf[0]["relativeOrbit"] == 172
    assert asf[0]["dataset"] == "NISAR" and asf[0]["processingLevel"] == "GCOV"


def test_layers_can_be_added_or_left_out(asf):
    layers = ("numberOfLooks", "rtcGammaToSigmaFactor")
    ds = nisar_gcov.load(AOI, "2026-10-21", "2026-11-08", ("HHHH",), layers=layers)
    assert set(ds.data_vars) == {"HHHH", *layers}
    assert float(ds.rtcGammaToSigmaFactor.max()) == 0.5
    # HH-only pass kept now that HV isn't asked for
    assert ds.sizes["time"] == 3


def test_no_frames_is_an_error(asf):
    with pytest.raises(ValueError, match="no validated NISAR GCOV frames with VVVV"):
        nisar_gcov.load(AOI, "2026-10-21", "2026-11-08", ("VVVV",))


def test_one_header_per_frame_grid(asf, monkeypatch):
    keys = []
    set_proj = earthdata.set_proj

    def spy(items, asset, key):
        keys.append({key(it) for it in items})
        set_proj(items, asset, key)

    monkeypatch.setattr(earthdata, "set_proj", spy)
    nisar_gcov.load(AOI, "2026-10-21", "2026-11-08")
    assert keys == [{(172, 64, "40+5"), (172, 65, "40+5")}]


def test_load_keeps_the_passes_the_windows_choose(asf):
    windows = parse_window(
        {"pre": {"before": "10d", "select": "last"}, "post": {"after": "10d"}},
        pd.Timestamp("2026-10-30"),
    )
    ds = nisar_gcov.load(AOI, "2026-10-20", "2026-11-09", windows=windows)
    assert [str(t)[:10] for t in ds.time.values] == ["2026-10-24", "2026-11-05"]
    assert ds.window.values.tolist() == ["pre", "post"]


def test_subdataset_is_frequency_a_of_the_remote_file():
    assert nisar_gcov.subdataset("https://x/f.h5?sig=1", "HHHH") == (
        'HDF5:"/vsicurl/https://x/f.h5?sig=1"://science/LSAR/GCOV/grids/frequencyA/HHHH'
    )


def test_to_item_has_an_asset_per_term_and_layer(frames):
    item = nisar_gcov.to_item(frames[3])  # HH only
    assert set(item.assets) == {
        "HHHH",
        "mask",
        "numberOfLooks",
        "rtcGammaToSigmaFactor",
    }
    assert {a.href for a in item.assets.values()} == {frames[3]["properties"]["url"]}
    assert item.properties["sat:relative_orbit"] == 5
    assert item.properties["nisar:frame"] == 64


def test_sign_requests_a_shared_file_once(monkeypatch, frames):
    calls = []

    class Response:
        status_code = 303

        def __init__(self):
            self.headers = {"location": "https://cdn.example/f.h5?sig=x"}

    monkeypatch.setattr(
        fetch,
        "get",
        lambda href, headers, allow_redirects: calls.append(href) or Response(),
    )
    monkeypatch.setenv("EARTHDATA_TOKEN", "secret")
    item = nisar_gcov.to_item(frames[0])
    earthdata.sign([item], ["HHHH", "HVHV", "mask"])
    assert len(calls) == 1
    assert item.assets["mask"].href == "https://cdn.example/f.h5?sig=x"


def test_the_sar_datasets_share_common_names():
    terms = {"HHHH": "gamma0_hh", "HVHV": "gamma0_hv"}
    assert {nisar_gcov.COMMON_NAMES[t] for t in terms} == set(terms.values())
    for module in (sentinel1_rtc, opera_rtc):
        assert module.COMMON_NAMES["vh"] == "gamma0_vh"
    ds = xr.Dataset(
        {
            "HVHV": ("x", [1.0], {"gaia:common-name": "gamma0_hv"}),
            "vh": ("x", [2.0], {"gaia:common-name": "gamma0_vh"}),
        }
    )
    with cf_xarray.set_options(custom_criteria=cf_criteria()):
        assert ds.cf["gamma0_hv"].name == "HVHV"
        assert ds.cf["gamma0_vh"].name == "vh"


def test_cli_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["nisar-gcov", "-i", str(aoi), "-s", "2026-07-01", "-e", "2026-10-20"]
        + ["-o", str(tmp_path / "o"), "-m", "provisional", "-t", "HVHV", "--dry-run"],
        result_action="return_value",
    )
    out = capsys.readouterr().out
    assert "maturity: provisional" in out and "terms: [HVHV]" in out
