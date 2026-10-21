import inspect
import typing

import numpy as np
import pandas as pd
import pytest
import requests
import xarray as xr

from gaia_cli import regrid, runner, stagefile
from gaia_cli.aoi import load_aoi, parse_aoi
from gaia_cli.cli import stage_app
from gaia_cli.datasets import hrrr, prism, synoptic
from gaia_cli.datasets.hrrr import HRRRVariable
from gaia_cli.window import parse_duration, parse_time, parse_window

T = pd.Timestamp("2023-09-14T02:06:49")

PETERS = """
name: peters-dome-2023
time: 2023-09-14T02:06:49Z
origin: {latitude: 63.1455, longitude: -151.1904}
attrs: {title: Peters Dome}
aois:
  regional: {bbox: [-152.2, 62.75, -150.2, 63.55]}
  local: {latitude: 63.1455, longitude: -151.1904, radius_km: 15}
aoi: regional
window: 2023-06-01/2023-09-21
output: s3://gaia/scratch/peters-dome
datasets:
  hrrr: {aoi: local, window: {before: 6h, after: 1d}, variable: temperature_2m}
  prism: {variables: [ppt, tmax]}
  nlcd: {year: 2021}
"""


def write(tmp_path, text):
    path = tmp_path / "stage.yaml"
    path.write_text(text)
    return str(path)


def problems(tmp_path, text):
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(write(tmp_path, text))
    return e.value.problems


# ---- window grammar ------------------------------------------------------


def test_date_only_end_includes_the_day_and_offsets_become_utc():
    (w,) = parse_window("2025-12-01/2025-12-20", None)
    assert (w.start, w.end) == (pd.Timestamp("2025-12-01"), pd.Timestamp("2025-12-21"))
    (w,) = parse_window({"start": "2025-12-01T00:00-08:00", "end": "2025-12-02"}, None)
    assert w.start == pd.Timestamp("2025-12-01T08:00")


def test_relative_and_mixed_windows():
    (w,) = parse_window({"before": "17min", "after": "23min"}, T)
    assert (w.start, w.end) == (
        T - pd.Timedelta(minutes=17),
        T + pd.Timedelta(minutes=23),
    )
    (w,) = parse_window({"start": "2023-06-01", "after": "1w"}, T)
    assert (w.start, w.end) == (pd.Timestamp("2023-06-01"), T + pd.Timedelta(days=7))


def test_pre_post_windows_default_to_time_on_the_missing_side():
    pre, post = parse_window(
        {"pre": {"before": "30d", "select": "last"}, "post": {"after": "30d"}}, T
    )
    assert (pre.name, pre.start, pre.end, pre.select) == (
        "pre",
        T - pd.Timedelta(days=30),
        T,
        "last",
    )
    assert (post.name, post.start) == ("post", T)


@pytest.mark.parametrize(
    "spec, message",
    [
        ({"before": "17min"}, "need a top-level 'time'"),
        ("2025-12-01", "not an interval"),
        ("2025-12-20/2025-12-01", "starts at"),
        ({"start": "2025-12-01", "before": "1d"}, "either start or before"),
        ({"strat": "2025-12-01"}, "did you mean 'start'"),
        ({"pre": {"before": "1d"}, "start": "2023-09-01"}, "unknown window key 'pre'"),
    ],
)
def test_bad_windows(spec, message):
    with pytest.raises(ValueError, match=message):
        parse_window(spec, None)


def test_durations():
    assert parse_duration("30d") == pd.Timedelta(days=30)
    assert parse_duration("PT17M") == pd.Timedelta(minutes=17)
    for bad in ("3M", "30m", "-1d"):  # no months; minutes are min
        with pytest.raises(ValueError):
            parse_duration(bad)


# ---- AOIs ------------------------------------------------------------------


def test_circle_aoi_has_true_radius():
    gdf = load_aoi({"latitude": 63.1455, "longitude": -151.1904, "radius_km": 400})
    utm = gdf.to_crs(gdf.estimate_utm_crs())
    assert utm.area.iloc[0] == pytest.approx(np.pi * 400e3**2, rel=0.01)


@pytest.mark.parametrize(
    "spec",
    [
        {"bbox": [10, 0, 5, 1]},
        {"latitude": 95, "longitude": 0, "radius_km": 1},
        {"intersects": {"type": "Blob"}},
        {
            "intersects": {
                "type": "MultiPolygon",
                "coordinates": [[[[0, 0], [1, 0], [1, 1], [0, 0]]]],
            }
        },
    ],
)
def test_bad_inline_aois(spec, tmp_path):
    with pytest.raises(ValueError):
        parse_aoi(spec, tmp_path)


def test_relative_aoi_paths_resolve_against_the_file(tmp_path):
    (tmp_path / "a.geojson").write_text("{}")
    cfg = stagefile.load_config(
        write(tmp_path, "name: x\naoi: a.geojson\ndatasets: {nlcd: {}}\n")
    )
    assert cfg.requests[0].aoi == str(tmp_path / "a.geojson")


# ---- stage files read from a URL -------------------------------------------

URL = "https://example.com/events/stage.yaml"


def serve(monkeypatch, text, status=200):
    """Answer the download of URL with *text*."""

    def get(url, timeout):
        assert url == URL
        r = requests.Response()
        r.status_code, r._content, r.url = status, text.encode(), url
        return r

    monkeypatch.setattr(stagefile.requests, "get", get)


def test_stage_file_from_a_url(monkeypatch, capsys):
    serve(monkeypatch, PETERS)
    runner.stage(config=URL, dry_run=True)
    assert "aoi_bbox: [-152.2, 62.75, -150.2, 63.55]" in capsys.readouterr().out


def test_stage_file_from_a_url_gives_aois_as_urls_or_inline(monkeypatch):
    serve(
        monkeypatch,
        "name: x\naoi: a.geojson\n"
        "aois: {b: /data/b.geojson, c: 'https://example.com/c.geojson'}\n"
        "datasets: {nlcd: {}}\n",
    )
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(URL)
    local = [p.split(":")[0] for p in e.value.problems if "local path" in p]
    assert local == ["aois.b", "aoi"]


def test_stage_file_url_that_fails(monkeypatch):
    serve(monkeypatch, "Not Found", status=404)
    with pytest.raises(stagefile.StageFileError, match="404"):
        stagefile.load_config(URL)


# ---- whole files -----------------------------------------------------------


def test_load_config_resolves_defaults_named_aois_and_options(tmp_path):
    cfg = stagefile.load_config(write(tmp_path, PETERS))
    reqs = {r.name: r for r in cfg.requests}
    assert reqs["hrrr"].aoi == {
        "latitude": 63.1455,
        "longitude": -151.1904,
        "radius_km": 15,
    }
    assert reqs["hrrr"].windows[0].start == T - pd.Timedelta(hours=6)
    assert reqs["hrrr"].options == {"variable": HRRRVariable.temperature_2m}
    assert reqs["prism"].aoi == {"bbox": [-152.2, 62.75, -150.2, 63.55]}
    assert str(reqs["prism"].windows[0]) == "2023-06-01T00:00:00/2023-09-22T00:00:00"
    assert reqs["prism"].options == {"variables": ("ppt", "tmax")}
    assert reqs["nlcd"].windows == []


def test_every_problem_is_reported_together(tmp_path):
    found = problems(
        tmp_path,
        """
tme: 2023-09-14
window: {before: 3M}
datasets:
  prsm: {}
  hrrr: {variable: temp, colour: red}
  nlcd: {year: "2021", window: 2021-01-01/2021-12-31}
""",
    )
    expected = [
        "unknown key 'tme'; did you mean 'time'?",
        "name: required",
        "window: '3M' is not a duration",
        "unknown dataset 'prsm'; did you mean 'prism'?",
        "hrrr: unknown option 'colour'",
        "hrrr: no aoi",
        "hrrr.variable: 'temp' is not one of",
        "nlcd: is static and takes no window",
        "nlcd.year: '2021' is not of type int",
    ]
    for text in expected:
        assert any(text in p for p in found), (text, found)


def test_prism_takes_only_its_variables(tmp_path):
    text = "name: x\naoi: {bbox: [0, 0, 1, 1]}\nwindow: 2025-12-01/2025-12-02\ndatasets:\n  prism: {variables: [ppt, tmx]}\n"
    assert any(
        "prism.variables: 'tmx' is not one of ppt, tmin, tmax" in p
        for p in problems(tmp_path, text)
    )
    path = write(tmp_path, text.replace("tmx", "tmin"))
    (req,) = stagefile.load_config(path).requests
    assert req.options == {"variables": ("ppt", "tmin")}


@pytest.mark.parametrize(
    "window, message",
    [
        ("{pre: {before: 1d}, post: {after: 1d}}", "takes one window, not pre/post"),
        ("{start: 2023-09-01, end: 2023-09-02, select: last}", "only for scene-based"),
        ("2023-09-01/..", "needs both a start and an end"),
    ],
)
def test_current_datasets_take_one_closed_window(tmp_path, window, message):
    text = f"name: x\ntime: 2023-09-14\naoi: {{bbox: [0, 0, 1, 1]}}\ndatasets:\n  prism: {{window: {window}}}\n"
    assert any(message in p for p in problems(tmp_path, text))


# ---- CLI and runner ----------------------------------------------------------


def test_data_flags_are_rejected_with_config(tmp_path):
    with pytest.raises(SystemExit, match="can't be combined with -s"):
        stage_app(["--config", write(tmp_path, PETERS), "-s", "2023-09-14"])
    with pytest.raises(SystemExit, match="can't be combined with --variables"):
        stage_app(["--config", write(tmp_path, PETERS), "--variables", "ppt"])


def test_dry_run_prints_resolved_config(tmp_path, capsys):
    runner.stage(config=write(tmp_path, PETERS), dry_run=True)
    out = capsys.readouterr().out
    assert "window: 2023-09-13T20:06:49/2023-09-15T02:06:49" in out
    assert "read: 2023-09-12T20:06:49/2023-09-16T02:06:49" in out
    assert "variable: temperature_2m" in out
    assert "aoi_bbox: [-152.2, 62.75, -150.2, 63.55]" in out
    # Peters Dome is in Alaska, outside the CONUS-only datasets
    assert "skipped: AOI [-152.2, 62.75, -150.2, 63.55] is outside" in out


# Around Bolt Creek, WA: inside every current dataset's extent
CONUS = """
name: bolt-creek
attrs: {title: Bolt Creek}
aoi: {bbox: [-121.6, 47.6, -121.0, 47.9]}
window: 2023-06-01/2023-09-21
datasets:
"""


def fake_load(name, calls, fail=False):
    """A load() for *name* that records its arguments in *calls*."""

    def load(aoi, start_date, end_date, variables: tuple[str, ...] = ("ppt",)):
        calls.append((name, start_date, end_date, variables))
        if fail:
            raise RuntimeError(f"{name} is down")
        return xr.Dataset({"v": ("time", [1.0])}, attrs={"title": name})

    return load


def test_run_writes_one_group_per_dataset(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    text = CONUS + "  prism: {variables: [ppt, tmax]}\n"
    out = tmp_path / "out"
    runner.stage(config=write(tmp_path, text), output_path=str(out))

    # Read a day beyond each end; the exclusive end becomes load()'s inclusive end
    assert calls == [
        ("prism", "2023-05-31T00:00:00", "2023-09-22T23:59:59", ("ppt", "tmax"))
    ]
    tree = xr.open_datatree(out / "bolt-creek.zarr", engine="zarr")
    assert tree.attrs["title"] == "Bolt Creek"
    assert tree["rasters/prism"].attrs["title"] == "prism"


# ---- registry and extent -----------------------------------------------------


@pytest.mark.parametrize("name", list(stagefile.DATASETS))
def test_registered_modules_declare_kind_and_extent(name):
    module = stagefile.DATASETS[name]
    assert module.KIND in stagefile.GROUPS
    # load() takes both window arguments (timed) or neither (static)
    params = inspect.signature(module.load).parameters
    assert stagefile.timed(module) == any(a in params for a in stagefile.WINDOW_ARGS)
    if hasattr(module, "REGRID_METHOD"):
        assert module.REGRID_METHOD in typing.get_args(regrid.Method)
    if hasattr(module, "VARIANT"):  # an option with a default, its values the keys
        assert module.VARIANT in stagefile.options(module.load)
    (overall, *_) = module.EXTENT["spatial"]["bbox"]
    assert len(overall) == 4
    for start, end in module.EXTENT["temporal"]["interval"]:
        for t in (start, end):
            assert t is None or parse_time(t)


def request(name, bbox, window="2023-06-01/2023-09-21"):
    aoi = {"bbox": bbox}
    return stagefile.Request(name, aoi, parse_window(window, None), {})


def test_outside_extent_checks_bbox_and_window():
    aoi = load_aoi
    bolt, peters = [-121.6, 47.6, -121.0, 47.9], [-152.2, 62.75, -150.2, 63.55]
    assert stagefile.outside_extent(request("prism", bolt), aoi({"bbox": bolt})) is None
    reason = stagefile.outside_extent(request("prism", peters), aoi({"bbox": peters}))
    assert reason is not None and "outside the dataset's extent" in reason
    assert (
        stagefile.outside_extent(request("synoptic", peters), aoi({"bbox": peters}))
        is None
    )
    early = request("hrrr", bolt, "2010-01-01/2010-01-31")
    reason = stagefile.outside_extent(early, aoi({"bbox": bolt}))
    assert reason is not None and "time extent 2014-10-01T00:00:00/.." in reason


def test_outside_extent_handles_antimeridian_and_sub_regions(monkeypatch):
    # An Aleutians-to-Alaska bbox crossing the antimeridian, plus two sub-regions
    extent = {
        "spatial": {
            "bbox": [
                [170.0, 50.0, -130.0, 72.0],
                [170.0, 50.0, -170.0, 56.0],
                [-160.0, 58.0, -140.0, 70.0],
            ]
        },
        "temporal": {"interval": [[None, None]]},
    }
    monkeypatch.setattr(prism, "EXTENT", extent)
    for bbox, inside in [
        ([175.0, 51.0, 176.0, 52.0], True),  # west of the antimeridian
        ([-151.0, 63.0, -150.0, 64.0], True),  # in the second sub-region
        ([-165.0, 59.0, -164.0, 60.0], False),  # in the overall bbox only
        ([-121.6, 47.6, -121.0, 47.9], False),
    ]:
        gdf = load_aoi({"bbox": bbox})
        assert (stagefile.outside_extent(request("prism", bbox), gdf) is None) == inside


def test_datasets_outside_their_extent_are_skipped(tmp_path, monkeypatch, caplog):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    monkeypatch.setattr(synoptic, "load", fake_load("synoptic", calls))
    monkeypatch.setattr(synoptic, "to_dataset", lambda ds: ds)
    text = CONUS.replace(
        "[-121.6, 47.6, -121.0, 47.9]", "[-152.2, 62.75, -150.2, 63.55]"
    )
    out = tmp_path / "out"
    runner.stage(
        config=write(tmp_path, text + "  prism: {}\n  synoptic: {}\n"),
        output_path=str(out),
    )
    assert [c[0] for c in calls] == ["synoptic"]
    assert "prism: skipped, AOI" in caplog.text
    tree = xr.open_datatree(out / "bolt-creek.zarr", engine="zarr")
    assert set(tree.children["stations"].children) == {"synoptic"}
    # The resolved file lists only what was staged
    doc = stagefile.read_config(str(out / "stage.yaml"))
    assert list(doc["datasets"]) == ["synoptic"]


def test_single_dataset_outside_its_extent_is_an_error(tmp_path):
    aoi = tmp_path / "alaska.geojson"
    load_aoi({"bbox": [-152.2, 62.75, -150.2, 63.55]}).to_file(aoi)
    args = ["prism", "-i", str(aoi), "-s", "2023-06-01", "-e", "2023-06-02"]
    with pytest.raises(SystemExit, match="nothing to stage .prism: AOI .* is outside"):
        stage_app([*args, "-o", str(tmp_path / "out")])


# ---- runner controls -----------------------------------------------------------


def test_only_selects_datasets_and_rejects_unknown_names(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    monkeypatch.setattr(hrrr, "load", fake_load("hrrr", calls))
    config = write(tmp_path, CONUS + "  prism: {}\n  nlcd: {}\n")
    with pytest.raises(SystemExit, match="unknown dataset in the file 'hrrr'"):
        runner.stage(config=config, only="prism,hrrr", dry_run=True)
    runner.stage(config=config, only="prism", output_path=str(tmp_path / "o"))
    assert [c[0] for c in calls] == ["prism"]


def test_failures_are_reported_and_retried_with_skip_existing(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls, fail=True))
    monkeypatch.setattr(hrrr, "load", fake_load("hrrr", calls))
    config = write(tmp_path, CONUS + "  prism: {}\n  hrrr: {}\n")
    out = str(tmp_path / "out")
    store = f"{out}/bolt-creek.zarr"
    with pytest.raises(SystemExit, match=r"1 of 2 dataset\(s\) failed; 1 written") as e:
        runner.stage(config=config, output_path=out)
    assert "prism: RuntimeError: prism is down" in str(e.value)
    assert set(xr.open_datatree(store, engine="zarr")["rasters"].children) == {"hrrr"}

    with pytest.raises(SystemExit, match="already exists"):
        runner.stage(config=config, output_path=out)
    with pytest.raises(SystemExit, match="can't be combined"):
        runner.stage(config=config, output_path=out, overwrite=True, skip_existing=True)

    calls.clear()
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    runner.stage(config=config, output_path=out, skip_existing=True)
    assert [c[0] for c in calls] == ["prism"]  # hrrr was already there
    tree = xr.open_datatree(store, engine="zarr")
    assert set(tree["rasters"].children) == {"hrrr", "prism"}
    assert tree.attrs["title"] == "Bolt Creek"


# ---- presets and flags -----------------------------------------------------------


def test_conus_precip_preset_takes_aoi_and_window_from_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"bbox": [-121.6, 47.6, -121.0, 47.9]}).to_file(aoi)
    stage_app(
        [
            "--config",
            "conus-precip",
            "-i",
            str(aoi),
            "-s",
            "2025-12-01",
            "-e",
            "2025-12-20",
            "--dry-run",
        ],
        result_action="return_value",
    )
    out = stagefile.yaml.safe_load(capsys.readouterr().out)
    assert list(out["datasets"]) == ["synoptic", "prism", "hrrr"]
    for entry in out["datasets"].values():
        assert entry["aoi"] == str(aoi)
        assert entry["window"] == "2025-12-01T00:00:00/2025-12-21T00:00:00"
        assert "skipped" not in entry


@pytest.mark.parametrize(
    "args, message",
    [
        (["--config", "conus-precip"], "the preset has no AOI"),
        (["--config", "conus-precipitation"], "not a preset .conus-precip"),
    ],
)
def test_preset_errors(args, message):
    with pytest.raises(SystemExit, match=message):
        stage_app([*args, "--dry-run"])


def test_flags_run_through_the_same_runner(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"bbox": [-121.6, 47.6, -121.0, 47.9]}).to_file(aoi)
    out = tmp_path / "skagit"
    stage_app(
        [
            "prism",
            "-i",
            str(aoi),
            "-s",
            "2025-12-01",
            "-e",
            "2025-12-20",
            "-v",
            "tmax",
            "-o",
            str(out),
        ],
        result_action="return_value",
    )
    assert calls == [("prism", "2025-11-30T00:00:00", "2025-12-21T23:59:59", ("tmax",))]
    # Named after the output folder
    tree = xr.open_datatree(out / "skagit.zarr", engine="zarr")
    assert tree["rasters/prism"].attrs["title"] == "prism"


# ---- output folder ---------------------------------------------------------------


def test_output_folder_holds_store_item_and_a_repeatable_stage_file(
    tmp_path, monkeypatch, capsys
):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"bbox": [-121.6, 47.6, -121.0, 47.9]}).to_file(aoi)
    text = CONUS.replace("{bbox: [-121.6, 47.6, -121.0, 47.9]}", "aoi.geojson")
    out = tmp_path / "event"
    runner.stage(
        config=write(tmp_path, text + "  prism: {variables: [tmax]}\n"),
        output_path=str(out),
    )
    assert sorted(p.name for p in out.iterdir()) == [
        "bolt-creek.json",
        "bolt-creek.zarr",
        "stage.yaml",
    ]
    item = stagefile.json.loads((out / "bolt-creek.json").read_text())
    assert item["assets"]["data"]["href"] == "./bolt-creek.zarr"
    assert item["assets"]["config"]["href"] == "./stage.yaml"
    assert item["assets"]["config"]["roles"] == ["metadata"]

    doc = stagefile.read_config(str(out / "stage.yaml"))
    assert doc["output"] == str(out)
    entry = doc["datasets"]["prism"]
    # The local AOI file is inlined, so the stage file stands on its own
    assert entry["aoi"]["intersects"]["type"] == "Polygon"
    assert entry["window"] == "2023-06-01T00:00:00/2023-09-22T00:00:00"
    assert entry["variables"] == ["tmax"]

    # Rerunning the resolved file resolves to the same configuration
    aoi.unlink()
    runner.stage(config=str(out / "stage.yaml"), dry_run=True)
    rerun = stagefile.yaml.safe_load(capsys.readouterr().out)["datasets"]["prism"]
    assert {k: rerun[k] for k in entry} == entry


def test_stage_file_keeps_url_aois():
    url = "https://example.com/basin.geojson"
    req = stagefile.Request("nlcd", url, [], {"year": 2021})
    cfg = stagefile.StageConfig("x", None, None, None, {}, [req])
    doc = stagefile.resolved(cfg, "out", {})
    assert doc["datasets"]["nlcd"] == {"aoi": url, "year": 2021}


def test_skip_existing_keeps_earlier_datasets_in_the_stage_file(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(prism, "load", fake_load("prism", calls))
    monkeypatch.setattr(hrrr, "load", fake_load("hrrr", calls))
    config = write(tmp_path, CONUS + "  prism: {}\n  hrrr: {}\n")
    out = str(tmp_path / "out")
    runner.stage(config=config, only="prism", output_path=out)
    runner.stage(config=config, only="hrrr", output_path=out, skip_existing=True)
    doc = stagefile.read_config(f"{out}/stage.yaml")
    assert list(doc["datasets"]) == ["prism", "hrrr"]


@pytest.mark.parametrize(
    "name, output, message",
    [
        ("bolt-creek", "out.zarr", "drop the .zarr from out.zarr"),
        ("bolt creek", "out", "'bolt creek' names the store"),
    ],
)
def test_output_and_name_must_make_a_folder_layout(tmp_path, name, output, message):
    config = write(tmp_path, CONUS.replace("bolt-creek", name) + "  prism: {}\n")
    with pytest.raises(SystemExit, match=message):
        runner.stage(config=config, output_path=output, dry_run=True)
