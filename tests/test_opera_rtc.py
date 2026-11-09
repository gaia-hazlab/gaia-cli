import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import box, mapping

from gaia_cli import earthdata
from gaia_cli.cli import stage_app
from gaia_cli.datasets import opera_rtc
from gaia_cli.datasets.sentinel1_rtc import OrbitState
from gaia_cli.stac import Match
from gaia_cli.window import parse_window

# 10 x 10 cells of 30 m in UTM 5N, near Peters Dome; each burst is half of it
X0, Y0 = 590010, 7010010
NORTH = (X0, Y0 - 150, X0 + 300, Y0)
SOUTH = (X0, Y0 - 300, X0 + 300, Y0 - 150)
BUCKET = "asf-cumulus-prod-opera-products"
COLLECTION = {
    "title": "OPERA RTC-S1",
    "license": "CC0-1.0",
    "links": [{"rel": "license", "href": "https://www.earthdata.nasa.gov/"}],
}
AOI = gpd.GeoDataFrame(
    geometry=[box(X0, Y0 - 300, X0 + 300, Y0)], crs="EPSG:32605"
).to_crs(4326)


def cog(path, bounds, value, dtype="float32", nodata=np.nan):
    """A GeoTIFF of *value* over *bounds* (5 x 10 cells of 30 m)."""
    left, bottom, right, top = bounds
    values = np.full((5, 10), value, dtype=dtype)
    da = xr.DataArray(
        values,
        dims=("y", "x"),
        coords={
            "y": np.arange(top - 15, bottom, -30.0),
            "x": np.arange(left + 15, right, 30.0),
        },
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    da.rio.write_crs("EPSG:32605").rio.write_nodata(nodata).rio.to_raster(path)


def feature(root, burst, bounds, value, when, orbit, track, processed="2025-02-14T15:39:35Z", pols=("VV", "VH"), mask=0):  # fmt: skip
    """An asf_search GeoJSON feature for one burst, with its COGs under *root*."""
    t = pd.Timestamp(when)
    name = (
        f"OPERA_L2_RTC-S1_{burst.replace('_', '-')}_{t:%Y%m%dT%H%M%SZ}_"
        f"{pd.Timestamp(processed):%Y%m%dT%H%M%SZ}_S1A_30_v1.0"
    )
    keys = [f"OPERA_L2_RTC-S1/{name}/{name}_{layer}.tif" for layer in (*pols, "mask")]
    for key, layer in zip(keys, (*pols, "mask"), strict=True):
        if layer == "mask":
            cog(root / "OPERA" / key, bounds, mask, "uint8", 255)
        else:
            cog(root / "OPERA" / key, bounds, value)
    footprint = gpd.GeoSeries([box(*bounds)], crs="EPSG:32605").to_crs(4326)[0]
    return {
        "type": "Feature",
        "geometry": mapping(footprint),
        "properties": {
            "sceneName": name,
            "startTime": when,
            "processingDate": processed,
            "platform": "Sentinel-1A",
            "flightDirection": "DESCENDING",
            "pathNumber": track,
            "orbit": orbit,
            "operaBurstID": burst,
            "polarization": list(pols),
            "s3Urls": [f"s3://{BUCKET}/{k}" for k in keys]
            + [f"s3://asf-cumulus-prod-opera-browse/OPERA_L2_RTC-S1/{name}/x.png"],
        },
    }


def static(root, burst, bounds, track, mask, angle, processed="2023-12-02T03:34:45Z"):
    """An asf_search feature for one burst ID's static layers, the COGs under *root*."""
    name = f"OPERA_L2_RTC-S1-STATIC_{burst.replace('_', '-')}_20140403_S1A_30_v1.0"
    layers = ("mask", "local_incidence_angle", "incidence_angle")
    keys = {
        layer: f"OPERA_L2_RTC-S1_STATIC/{name}/{name}_{layer}.tif" for layer in layers
    }
    for layer, key in keys.items():
        if layer == "mask":
            cog(root / "OPERA" / key, bounds, mask, "uint8", 255)
        else:
            cog(root / "OPERA" / key, bounds, angle)
    footprint = gpd.GeoSeries([box(*bounds)], crs="EPSG:32605").to_crs(4326)[0]
    return {
        "type": "Feature",
        "geometry": mapping(footprint),
        "properties": {
            "sceneName": name,
            "startTime": None,
            "processingDate": processed,
            "platform": "Sentinel-1A",
            "flightDirection": "DESCENDING" if track != 65 else "ASCENDING",
            "pathNumber": track,
            "orbit": 11403,
            "operaBurstID": burst,
            "s3Urls": [f"s3://{BUCKET}/{k}" for k in keys.values()],
        },
    }


@pytest.fixture
def statics(tmp_path):
    """Static layers of the bursts below, and of a track with no passes (65)."""
    return [
        static(tmp_path, "T131_279928_IW2", NORTH, 131, mask=2, angle=40.0),
        static(tmp_path, "T131_279929_IW2", SOUTH, 131, mask=0, angle=41.0),
        static(tmp_path, "T102_217000_IW1", NORTH, 102, mask=1, angle=30.0),
        static(tmp_path, "T065_137842_IW1", NORTH, 65, mask=3, angle=20.0),
    ]


@pytest.fixture
def bursts(tmp_path):
    """Two bursts of one pass, a reprocessed burst, a VV-only pass and a later pass."""
    return [
        feature(
            tmp_path,
            "T131_279928_IW2",
            NORTH,
            1.0,
            "2023-09-09T16:29:31Z",
            50253,
            131,
            mask=2,
        ),
        feature(
            tmp_path, "T131_279929_IW2", SOUTH, 2.0, "2023-09-09T16:29:34Z", 50253, 131
        ),
        feature(
            tmp_path,
            "T131_279929_IW2",
            SOUTH,
            9.0,
            "2023-09-09T16:29:34Z",
            50253,
            131,
            processed="2026-01-01T00:00:00Z",
        ),
        feature(
            tmp_path,
            "T102_217000_IW1",
            NORTH,
            5.0,
            "2023-09-19T16:45:46Z",
            50399,
            102,
            pols=("VV",),
        ),
        feature(
            tmp_path, "T131_279928_IW2", NORTH, 3.0, "2023-09-21T16:29:31Z", 50428, 131
        ),
    ]


@pytest.fixture
def asf(monkeypatch, tmp_path, bursts, statics):
    """ASF stand-ins, the files local; returns the searches made."""
    searches = []

    def search(aoi, start=None, end=None, **filters):
        searches.append({"start": start, "end": end, **filters})
        if filters["processingLevel"] == "RTC-STATIC":
            tracks = filters["relativeOrbit"]
            return [f for f in statics if f["properties"]["pathNumber"] in tracks]
        return bursts

    monkeypatch.setattr(earthdata, "search", search)
    monkeypatch.setattr(earthdata, "HOST", str(tmp_path))
    monkeypatch.setattr(earthdata, "sign", lambda items, assets: None)

    def get(url, timeout):
        assert url in (opera_rtc.COLLECTION_URL, opera_rtc.STATIC_URL)
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(opera_rtc.provenance.requests, "get", get)
    return searches


def test_load_mosaics_the_bursts_of_each_pass(asf, caplog):
    caplog.set_level(logging.INFO)
    ds = opera_rtc.load(AOI, "2023-09-01T00:00:00", "2023-09-30T23:59:59")
    # The VV-only pass is skipped; the older processing of a burst too
    assert [str(t)[:10] for t in ds.time.values] == ["2023-09-09", "2023-09-21"]
    assert "skipped 1 older processing(s)" in caplog.text
    assert "skipped 1 burst(s) without vv+vh" in caplog.text
    assert ds.rio.crs.to_epsg() == 32605 and ds.rio.resolution() == (30.0, -30.0)
    first = ds.isel(time=0)
    assert float(first.vv.sel(y=Y0 - 15, x=X0 + 15)) == 1.0  # north burst
    assert float(first.vv.sel(y=Y0 - 285, x=X0 + 15)) == 9.0  # newest south burst
    assert ds.vv.attrs["units"] == "1" and ds.vv.dtype == "float32"
    # The mask, by default, is the static one of the pass's track
    assert ds.mask.dtype == "uint8"
    assert int(first.mask.sel(y=Y0 - 15, x=X0 + 15)) == 2  # layover
    assert int(first.mask.sel(y=Y0 - 285, x=X0 + 15)) == 0  # valid
    assert asf[1]["processingLevel"] == "RTC-STATIC" and asf[1]["relativeOrbit"] == [
        131
    ]
    assert ds.mask.attrs["flag_values"] == [0, 1, 2, 3, 255]
    assert ds.mask.attrs["flag_meanings"].split()[-1] == "no_data"
    assert ds.burst_ids.values.tolist() == [
        "T131_279928_IW2 T131_279929_IW2",
        "T131_279928_IW2",
    ]
    assert ds.orbit_state.values.tolist() == ["descending", "descending"]
    assert ds.relative_orbit.values.tolist() == [131, 131]
    assert ds.absolute_orbit.values.tolist() == [50253, 50428]
    assert ds.attrs["license"] == "CC0-1.0"
    derived = [ln["href"] for ln in ds.attrs["links"] if ln["rel"] == "derived_from"]
    assert len(derived) == 1 + 3 + 2  # the collection, each burst and static read


def test_orbit_options_filter_the_search(asf):
    opera_rtc.load(AOI, "2023-09-01", "2023-09-30", ("vv",), OrbitState.descending, 131)
    search = asf[0]
    assert search["flightDirection"] == "DESCENDING"
    assert search["relativeOrbit"] == 131
    assert search["dataset"] == "OPERA-S1" and search["processingLevel"] == "RTC"


def test_each_pass_gets_its_tracks_static_layers(asf):
    static = ("mask", "local_incidence_angle")
    ds = opera_rtc.load(AOI, "2023-09-01", "2023-09-30", ("vv",), static=static).load()
    assert ds.relative_orbit.values.tolist() == [131, 102, 131]
    assert asf[1]["relativeOrbit"] == [102, 131]  # track 65 has no passes
    north = ds.sel(y=Y0 - 15, x=X0 + 15)
    assert north["mask"].values.tolist() == [2, 1, 2]
    assert north.local_incidence_angle.values.tolist() == [40.0, 30.0, 40.0]
    assert ds.local_incidence_angle.attrs["units"] == "degree"
    assert "_FillValue" not in ds["mask"].attrs  # 255 stays a flag value
    # Track 102's static layers cover only the north half
    south = ds["mask"].sel(time=ds.relative_orbit == 102, y=Y0 - 285, x=X0 + 15)
    assert south.values.tolist() == [255]


def test_static_can_be_left_out(asf):
    ds = opera_rtc.load(AOI, "2023-09-01", "2023-09-30", static=())
    assert "mask" not in ds and len(asf) == 1


def test_a_track_without_static_layers_is_an_error(asf, statics):
    statics[:] = [f for f in statics if f["properties"]["pathNumber"] != 102]
    with pytest.raises(ValueError, match="for relative orbit\\(s\\) 102"):
        opera_rtc.load(AOI, "2023-09-01", "2023-09-30", ("vv",))


def test_no_bursts_is_an_error(asf):
    with pytest.raises(ValueError, match="no OPERA RTC bursts with hh"):
        opera_rtc.load(AOI, "2023-09-01", "2023-09-30", ("hh",))


def test_load_keeps_the_passes_the_windows_choose(asf):
    windows = parse_window(
        {"pre": {"before": "30d", "select": "last"}, "post": {"after": "30d"}},
        pd.Timestamp("2023-09-14T02:06:49"),
    )
    ds = opera_rtc.load(
        AOI,
        "2023-08-15",
        "2023-10-14",
        ("vv",),
        windows=windows,
        match=Match.relative_orbit,
    )
    # Track 102 has no pass before the event, so only track 131 is kept
    assert ds.relative_orbit.values.tolist() == [131, 131]
    assert ds.window.values.tolist() == ["pre", "post"]


def test_cli_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["opera-rtc", "-i", str(aoi), "-s", "2023-09-01", "-e", "2023-09-30"]
        + ["-o", str(tmp_path / "o"), "-p", "vv", "-r", "131", "--dry-run"],
        result_action="return_value",
    )
    out = capsys.readouterr().out
    assert "polarizations: [vv]" in out and "relative_orbit: 131" in out


# ---- earthdata.py --------------------------------------------------------------


def test_to_item_names_assets_by_layer_on_the_cumulus_host(tmp_path, bursts):
    item = earthdata.to_item(bursts[0], opera_rtc.BANDS)
    assert sorted(item.assets) == ["mask", "vh", "vv"]  # browse PNGs are left out
    name = bursts[0]["properties"]["sceneName"]
    assert item.assets["vv"].href == (
        f"{earthdata.HOST}/OPERA/OPERA_L2_RTC-S1/{name}/{name}_VV.tif"
    )
    assert item.assets["vv"].extra_fields["raster:bands"] == [
        {"data_type": "float32", "nodata": "nan"}
    ]
    assert "raster:bands" not in item.assets["mask"].extra_fields  # not loaded
    assert item.properties["sat:orbit_state"] == "descending"
    assert item.properties["sat:relative_orbit"] == 131
    assert item.properties["opera:burst_id"] == "T131_279928_IW2"
    assert str(item.datetime) == "2023-09-09 16:29:31+00:00"


def test_static_layers_take_their_processing_time(bursts):
    static = {**bursts[0], "properties": {**bursts[0]["properties"], "startTime": None}}
    assert str(earthdata.to_item(static).datetime) == "2025-02-14 15:39:35+00:00"


def test_newest_keeps_the_latest_processing_of_each_burst(bursts):
    items = [earthdata.to_item(f) for f in bursts]
    kept = earthdata.newest(items)
    assert len(kept) == 4
    assert items[1] not in kept and items[2] in kept


class Response:
    def __init__(self, status, location=""):
        self.status_code, self.headers = status, {"location": location}


def test_sign_points_assets_at_the_signed_urls(monkeypatch, bursts):
    calls = []

    def get(href, headers, allow_redirects, timeout):
        calls.append((href, headers))
        return Response(303, f"https://cdn.example/{href.rsplit('/', 1)[-1]}?sig=x")

    monkeypatch.setattr(earthdata.requests, "get", get)
    monkeypatch.setenv("EARTHDATA_TOKEN", "secret")
    item = earthdata.to_item(bursts[0])
    earthdata.sign([item], ["vv", "mask"])
    assert item.assets["vv"].href.startswith("https://cdn.example/")
    assert item.assets["vv"].href.endswith("_VV.tif?sig=x")
    assert item.assets["vh"].href.startswith(earthdata.HOST)  # not asked for
    assert len(calls) == 2
    assert {h["Authorization"] for _, h in calls} == {"Bearer secret"}


@pytest.mark.parametrize(
    "response",
    [Response(302, "https://urs.earthdata.nasa.gov/oauth/authorize"), Response(404)],
)
def test_sign_errors_name_the_token(monkeypatch, bursts, response):
    monkeypatch.setattr(earthdata.requests, "get", lambda *a, **k: response)
    monkeypatch.setenv("EARTHDATA_TOKEN", "expired")
    with pytest.raises(RuntimeError, match="is EARTHDATA_TOKEN current"):
        earthdata.sign([earthdata.to_item(bursts[0])], ["vv"])


def test_a_missing_token_is_an_error(monkeypatch):
    monkeypatch.delenv("EARTHDATA_TOKEN", raising=False)
    with pytest.raises(RuntimeError, match="EARTHDATA_TOKEN is not set"):
        earthdata.token()


def test_set_proj_reads_one_header_per_burst(tmp_path, monkeypatch, bursts):
    monkeypatch.setattr(earthdata, "HOST", str(tmp_path))
    opened = []
    open_ = earthdata.rasterio.open
    monkeypatch.setattr(
        earthdata.rasterio, "open", lambda h: opened.append(h) or open_(h)
    )
    items = [earthdata.to_item(f) for f in bursts]
    earthdata.set_proj(items, "vv")
    assert len(opened) == 3  # three burst IDs
    assert {it.properties["proj:code"] for it in items} == {"EPSG:32605"}
    assert items[0].properties["proj:transform"][:3] == [30.0, 0.0, X0]


def test_an_unreadable_burst_is_an_error(tmp_path, monkeypatch, bursts):
    monkeypatch.setattr(earthdata, "HOST", str(tmp_path / "missing"))
    with pytest.raises(RuntimeError, match="could not read"):
        earthdata.set_proj([earthdata.to_item(bursts[0])], "vv")


def test_search_passes_utc_times_and_the_aoi(monkeypatch):
    seen = {}

    class Results:
        def geojson(self):
            return {"features": ["f"]}

    def search(**kwargs):
        seen.update(kwargs)
        return Results()

    monkeypatch.setattr(earthdata.asf_search, "search", search)
    out = earthdata.search(AOI, "2023-09-01", "2023-09-30T23:59:59", dataset="OPERA-S1")
    assert out == ["f"]
    assert seen["start"] == "2023-09-01T00:00:00+00:00"
    assert seen["end"] == "2023-09-30T23:59:59+00:00"
    # The AOI's bounds, counterclockwise
    assert seen["intersectsWith"] == box(*AOI.total_bounds).wkt
    assert seen["dataset"] == "OPERA-S1"
