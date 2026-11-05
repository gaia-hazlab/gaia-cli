import datetime as dt
import logging

import dask.array
import geopandas as gpd
import numpy as np
import pandas as pd
import pystac
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import yaml
from shapely.geometry import box, mapping

from gaia_cli import runner, stac, stagefile
from gaia_cli.aoi import load_aoi
from gaia_cli.cli import stage_app
from gaia_cli.datasets import sentinel1_rtc
from gaia_cli.datasets.sentinel1_rtc import Match, OrbitState
from gaia_cli.window import parse_window

# 10 x 10 cells of 10 m in UTM 5N, near Peters Dome; each frame is half of it
X0, Y0 = 590000, 7010000
NORTH = (X0, Y0 - 50, X0 + 100, Y0)
SOUTH = (X0, Y0 - 100, X0 + 100, Y0 - 50)
COLLECTION = {
    "title": "Sentinel 1 RTC",
    "license": "CC-BY-4.0",
    "links": [
        {"rel": "license", "href": "https://creativecommons.org/licenses/by/4.0/"}
    ],
}


def frame(tmp_path, name, bounds, value):
    """A GeoTIFF of *value* over *bounds*, with one -32768 and one 0 cell."""
    left, bottom, right, top = bounds
    values = np.full((5, 10), value, dtype="float32")
    values[0, 0], values[0, 1] = -32768, 0
    da = xr.DataArray(
        values,
        dims=("y", "x"),
        coords={
            "y": np.arange(top - 5, bottom, -10.0),
            "x": np.arange(left + 5, right, 10.0),
        },
    )
    path = tmp_path / f"{name}.tif"
    da.rio.write_crs("EPSG:32605").rio.write_nodata(-32768).rio.to_raster(path)
    return str(path)


def item(tmp_path, name, bounds, value, when, orbit, track, pols=("vv", "vh")):
    footprint = gpd.GeoSeries([box(*bounds)], crs="EPSG:32605").to_crs(4326)[0]
    left, top = bounds[0], bounds[3]
    it = pystac.Item(
        id=name,
        geometry=mapping(footprint),
        bbox=list(footprint.bounds),
        datetime=dt.datetime.fromisoformat(when),
        properties={
            "platform": "SENTINEL-1A",
            "sat:orbit_state": "descending",
            "sat:relative_orbit": track,
            "sat:absolute_orbit": orbit,
            "proj:epsg": 32605,
            "proj:shape": [5, 10],
            "proj:transform": [10.0, 0.0, left, 0.0, -10.0, top, 0.0, 0.0, 1.0],
        },
    )
    for pol in pols:
        href = frame(tmp_path, f"{name}-{pol}", bounds, value)
        it.add_asset(pol, pystac.Asset(href, media_type=pystac.MediaType.COG))
    return it


@pytest.fixture
def scenes(tmp_path):
    """Two frames of one pass, a VV-only scene, and a later pass."""
    return [
        item(tmp_path, "a-north", NORTH, 1.0, "2023-09-09T16:29:15", 50253, 131),
        item(tmp_path, "a-south", SOUTH, 2.0, "2023-09-09T16:29:40", 50253, 131),
        item(tmp_path, "b-vv", NORTH, 5.0, "2023-09-19T16:45:36", 50399, 102, ("vv",)),
        item(tmp_path, "c", NORTH, 3.0, "2023-09-21T16:29:15", 50428, 131),
    ]


@pytest.fixture
def pc(monkeypatch, scenes):
    """Planetary Computer stand-ins; returns the searches made."""
    searches = []

    class Client:
        def search(self, **kwargs):
            searches.append(kwargs)
            return type("Search", (), {"items": lambda self: iter(scenes)})()

    monkeypatch.setattr(stac.pystac_client.Client, "open", lambda url: Client())
    monkeypatch.setattr(stac, "token", lambda collection: "sig=x")

    def get(url, timeout):
        assert url == sentinel1_rtc.COLLECTION_URL
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(sentinel1_rtc.provenance.requests, "get", get)
    return searches


AOI = gpd.GeoDataFrame(
    geometry=[box(X0, Y0 - 100, X0 + 100, Y0)], crs="EPSG:32605"
).to_crs(4326)


def test_load_mosaics_the_frames_of_each_pass(pc, caplog):
    caplog.set_level(logging.INFO)
    ds = sentinel1_rtc.load(AOI, "2023-09-01T00:00:00", "2023-09-30T23:59:59")

    assert isinstance(ds.vv.data, dask.array.Array)  # lazy: the size cap comes first
    assert list(ds.data_vars) == ["vv", "vh"]
    assert ds.vv.dtype == "float32"
    assert ds.rio.crs.to_epsg() == 32605
    # A time step per pass, at its first frame; the VV-only scene is skipped
    assert list(ds.time.values) == list(
        pd.to_datetime(["2023-09-09T16:29:15", "2023-09-21T16:29:15"])
    )
    assert "skipped 1 scene(s) without vv+vh" in caplog.text
    assert ds.relative_orbit.values.tolist() == [131, 131]
    assert ds.absolute_orbit.values.tolist() == [50253, 50428]
    assert ds.orbit_state.values.tolist() == ["descending", "descending"]
    assert ds.platform.values.tolist() == ["SENTINEL-1A"] * 2
    assert ds.scene_ids.values.tolist() == ["a-north a-south", "c"]

    vv = ds.vv.sel(time="2023-09-09").squeeze().values
    north = vv[: vv.shape[0] // 2]
    south = vv[vv.shape[0] // 2 :]
    # nodata and zero are NaN; the two frames fill one grid
    assert np.isnan(north[0, :2]).all() and np.isnan(south[0, :2]).all()
    assert np.nanmax(north) == 1.0 and np.nanmin(south) == 2.0
    assert np.isnan(ds.vv.sel(time="2023-09-21").squeeze().values[-1]).all()
    assert ds.vv.attrs["units"] == "1" and "nodata" not in ds.vv.attrs

    assert ds.attrs["license"] == "CC-BY-4.0"
    assert ds.attrs["description"].startswith("Sentinel-1 IW GRD backscatter")
    assert ds.attrs["start_datetime"] == "2023-09-09T16:29:15Z"
    hrefs = [lk["href"] for lk in ds.attrs["links"] if lk["rel"] == "derived_from"]
    assert hrefs[0] == sentinel1_rtc.COLLECTION_URL
    assert hrefs[1:] == [
        f"{sentinel1_rtc.COLLECTION_URL}/items/{i}" for i in ("a-north", "a-south", "c")
    ]


def test_orbit_options_filter_the_search(pc):
    sentinel1_rtc.load(
        AOI, "2023-09-01", "2023-09-30", ("vv",), OrbitState.descending, 131
    )
    (search,) = pc
    assert search["query"] == {
        "sat:orbit_state": {"eq": "descending"},
        "sat:relative_orbit": {"eq": 131},
    }
    assert search["collections"] == ["sentinel-1-rtc"]
    assert search["datetime"] == "2023-09-01/2023-09-30"


def test_vv_alone_keeps_the_vv_only_scene(pc):
    ds = sentinel1_rtc.load(AOI, "2023-09-01", "2023-09-30", ("vv",))
    assert ds.sizes["time"] == 3
    assert pc[0]["query"] is None


def test_no_scenes_is_an_error(pc):
    with pytest.raises(ValueError, match="no Sentinel-1 RTC scenes with hh"):
        sentinel1_rtc.load(AOI, "2023-09-01", "2023-09-30", ("hh",))


def test_cli_flags_and_stage_file_options(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"bbox": [-151.3, 63.1, -151.1, 63.2]}).to_file(aoi)
    stage_app(
        ["sentinel1-rtc", "-i", str(aoi), "-s", "2023-09-01", "-e", "2023-09-30"]
        + ["-o", str(tmp_path / "o"), "-p", "vv", "--orbit-state", "ascending"]
        + ["-r", "131", "--dry-run"],
        result_action="return_value",
    )
    out = capsys.readouterr().out
    assert "polarizations: [vv]" in out
    assert "orbit_state: ascending" in out and "relative_orbit: 131" in out

    path = tmp_path / "stage.yaml"
    path.write_text(
        "name: x\naoi: {bbox: [-151.3, 63.1, -151.1, 63.2]}\nwindow: 2023-09-01/2023-09-30\n"
        "datasets:\n  sentinel1-rtc: {polarizations: [hh, hv], orbit_state: up}\n"
    )
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(str(path))
    assert any(
        "'up' is not one of any, ascending, descending" in p for p in e.value.problems
    )


# ---- stac.py -----------------------------------------------------------------


def test_token_is_cached_until_it_nears_expiry(monkeypatch):
    calls = []
    expiry = [pd.Timestamp.now(tz="UTC") + pd.Timedelta(hours=24)]

    def get(url, timeout):
        calls.append(url)
        doc = {"token": f"sig={len(calls)}", "msft:expiry": expiry[0].isoformat()}
        return type(
            "R", (), {"raise_for_status": lambda self: None, "json": lambda self: doc}
        )()

    monkeypatch.setattr(stac.requests, "get", get)
    monkeypatch.setattr(stac, "_tokens", {})
    assert stac.token("sentinel-1-rtc") == stac.token("sentinel-1-rtc") == "sig=1"
    assert calls == [f"{stac.PC_SAS}/sentinel-1-rtc"]

    expiry[0] = pd.Timestamp.now(tz="UTC") + pd.Timedelta(minutes=30)
    stac._tokens.clear()
    stac.token("sentinel-1-rtc")  # sig=2, expires within TOKEN_MARGIN
    assert stac.token("sentinel-1-rtc") == "sig=3"


def test_sign_adds_the_token_to_blob_hrefs_only(monkeypatch):
    monkeypatch.setattr(stac, "token", lambda collection: "sig=x")
    it = pystac.Item("a", None, None, dt.datetime(2023, 9, 9), {})
    blob = "https://acct.blob.core.windows.net/c/a.tif"
    it.add_asset("vv", pystac.Asset(blob))
    it.add_asset("signed", pystac.Asset(f"{blob}?sig=old"))
    it.add_asset("thumbnail", pystac.Asset("https://example.com/a.png"))
    (signed,) = stac.sign([it], "sentinel-1-rtc")
    assert signed.assets["vv"].href == f"{blob}?sig=x"
    assert signed.assets["signed"].href == f"{blob}?sig=old"
    assert signed.assets["thumbnail"].href == "https://example.com/a.png"
    assert it.assets["vv"].href == blob  # the input is unchanged


def test_native_crs_is_the_most_common(scenes):
    other = scenes[0].clone()
    other.properties.pop("proj:epsg")
    other.properties["proj:code"] = "EPSG:32606"
    assert stac.native_crs([other, *scenes]) == "EPSG:32605"
    assert stac.item_crs(other) == "EPSG:32606"


# ---- choosing scenes (stage-file windows) --------------------------------------

EVENT = pd.Timestamp("2023-09-14T02:06:49")
PRE_POST = {
    "pre": {"before": "30d", "select": "last"},
    "post": {"after": "30d", "select": "first"},
}


def scene(when, track, name=None):
    it = pystac.Item(
        name or f"{when}-{track}",
        None,
        None,
        dt.datetime.fromisoformat(when),
        {"sat:relative_orbit": track},
    )
    return [it]


# Peters Dome, Aug-Oct 2023: tracks 102 and 131, with one 102 pass only before
PETERS_PASSES = [
    scene("2023-08-28T16:29", 131),
    scene("2023-09-07T16:45", 102),
    scene("2023-09-09T16:29", 131),
    scene("2023-09-19T16:45", 102),
    scene("2023-09-21T16:29", 131),
    scene("2023-10-03T16:29", 131),
    scene("2023-08-14T03:00", 160),  # a third track seen only before
]


def chosen(picked):
    return [
        (name, s[0].datetime.strftime("%m-%d"), s[0].properties["sat:relative_orbit"])
        for name, s in picked
    ]


def test_select_last_before_and_first_after():
    windows = parse_window(PRE_POST, EVENT)
    assert chosen(stac.select(PETERS_PASSES, windows)) == [
        ("pre", "09-09", 131),
        ("post", "09-19", 102),  # a different track: what match prevents
    ]


def test_match_pairs_scenes_on_each_track_that_has_both():
    windows = parse_window(PRE_POST, EVENT)
    assert chosen(stac.select(PETERS_PASSES, windows, "sat:relative_orbit")) == [
        ("pre", "09-07", 102),
        ("pre", "09-09", 131),
        ("post", "09-19", 102),
        ("post", "09-21", 131),
    ]  # track 160 has no post scene, so it is dropped


@pytest.mark.parametrize(
    "spec, expected",
    [
        ({"before": "30d", "after": "30d"}, 7 - 1),  # all, minus 08-14 (31 days before)
        (
            {"before": "30d", "after": "30d", "select": "nearest", "count": 2},
            ["09-09", "09-19"],
        ),
        ({"before": "30d", "select": "last", "count": 2}, ["09-07", "09-09"]),
        ({"after": "30d", "select": "first", "count": 3}, ["09-19", "09-21", "10-03"]),
    ],
)
def test_select_all_nearest_and_count(spec, expected):
    picked = stac.select(PETERS_PASSES, parse_window(spec, EVENT))
    if isinstance(expected, int):
        assert len(picked) == expected and all(name is None for name, _ in picked)
    else:
        assert [s[0].datetime.strftime("%m-%d") for _, s in picked] == expected


def test_overlapping_windows_count_a_scene_once():
    windows = parse_window(
        {
            "pre": {"before": "10d", "after": "10d"},
            "post": {"before": "10d", "after": "10d"},
        },
        EVENT,
    )
    picked = stac.select(PETERS_PASSES, windows)
    ids = [s[0].id for _, s in picked]
    assert len(ids) == len(set(ids)) == 4
    assert {name for name, _ in picked} == {"pre"}


def test_load_keeps_the_passes_the_windows_choose(pc):
    windows = parse_window(PRE_POST, EVENT)
    ds = sentinel1_rtc.load(AOI, "2023-08-15", "2023-10-14", windows=windows)
    assert ds.window.values.tolist() == ["pre", "post"]
    assert ds.scene_ids.values.tolist() == ["a-north a-south", "c"]
    assert [
        lk["href"].rsplit("/", 1)[-1]
        for lk in ds.attrs["links"]
        if "/items/" in lk["href"]
    ] == ["a-north", "a-south", "c"]


def test_load_errors_when_no_track_has_both_windows(pc):
    windows = parse_window({"pre": {"before": "1d"}, "post": {"after": "30d"}}, EVENT)
    with pytest.raises(ValueError, match="on any track with passes in every window"):
        sentinel1_rtc.load(
            AOI, "2023-09-13", "2023-10-14", match=Match.relative_orbit, windows=windows
        )


STAGE = """name: peters
time: 2023-09-14T02:06:49Z
aoi: {bbox: [-151.3, 63.1, -151.1, 63.2]}
datasets:
  sentinel1-rtc:
    window:
      pre: {before: 30d, select: last}
      post: {after: 30d, select: first}
    match: relative_orbit
"""


def test_stage_file_takes_pre_post_select_and_match(tmp_path):
    path = tmp_path / "stage.yaml"
    path.write_text(STAGE)
    (req,) = stagefile.load_config(str(path)).requests
    assert [w.name for w in req.windows] == ["pre", "post"]
    assert req.options["match"] is Match.relative_orbit
    assert stagefile.read_window(req).start == EVENT - pd.Timedelta(days=30)


def test_resolved_stage_file_repeats_the_windows(tmp_path, capsys):
    path = tmp_path / "stage.yaml"
    path.write_text(STAGE)
    stage_app(["--config", str(path), "--dry-run"], result_action="return_value")
    out = capsys.readouterr().out
    assert "read: 2023-08-15T02:06:49/2023-10-14T02:06:49" in out  # exact, not padded
    doc = yaml.safe_load(out)
    window = doc["datasets"]["sentinel1-rtc"]["window"]
    assert window["pre"] == {
        "start": "2023-08-15T02:06:49",
        "end": "2023-09-14T02:06:49",
        "select": "last",
    }
    # The resolved windows parse back to the same windows
    again = parse_window(window, pd.Timestamp(doc["time"]))
    assert [(w.name, w.start, w.end, w.select) for w in again] == [
        (w.name, w.start, w.end, w.select)
        for w in stagefile.load_config(str(path)).requests[0].windows
    ]


def test_runner_passes_the_windows_and_their_exact_span(tmp_path, monkeypatch):
    calls = []

    def load(aoi, start_date, end_date, match: Match = Match.none, windows=None):
        calls.append((start_date, end_date, match, [w.name for w in windows]))
        return xr.Dataset({"vv": ("time", [1.0])}, attrs={"title": "s1"})

    monkeypatch.setattr(sentinel1_rtc, "load", load)
    path = tmp_path / "stage.yaml"
    path.write_text(STAGE)
    runner.stage(config=str(path), output_path=str(tmp_path / "out"))
    assert calls == [
        (
            "2023-08-15T02:06:49",
            "2023-10-14T02:06:48",
            Match.relative_orbit,
            ["pre", "post"],
        )
    ]
