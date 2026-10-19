import numpy as np
import pandas as pd
import pytest
import xarray as xr

from gaia_cli import stagefile
from gaia_cli.cli import stage_app
from gaia_cli.hrrr import HRRRVariable

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
    (w,) = stagefile.parse_window("2025-12-01/2025-12-20", None)
    assert (w.start, w.end) == (pd.Timestamp("2025-12-01"), pd.Timestamp("2025-12-21"))
    (w,) = stagefile.parse_window(
        {"start": "2025-12-01T00:00-08:00", "end": "2025-12-02"}, None
    )
    assert w.start == pd.Timestamp("2025-12-01T08:00")


def test_relative_and_mixed_windows():
    (w,) = stagefile.parse_window({"before": "17min", "after": "23min"}, T)
    assert (w.start, w.end) == (
        T - pd.Timedelta(minutes=17),
        T + pd.Timedelta(minutes=23),
    )
    (w,) = stagefile.parse_window({"start": "2023-06-01", "after": "1w"}, T)
    assert (w.start, w.end) == (pd.Timestamp("2023-06-01"), T + pd.Timedelta(days=7))


def test_pre_post_windows_default_to_time_on_the_missing_side():
    pre, post = stagefile.parse_window(
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
        stagefile.parse_window(spec, None)


def test_durations():
    assert stagefile.parse_duration("30d") == pd.Timedelta(days=30)
    assert stagefile.parse_duration("PT17M") == pd.Timedelta(minutes=17)
    for bad in ("3M", "30m", "-1d"):  # no months; minutes are min
        with pytest.raises(ValueError):
            stagefile.parse_duration(bad)


# ---- AOIs ------------------------------------------------------------------


def test_circle_aoi_has_true_radius():
    gdf = stagefile.load_aoi(
        {"latitude": 63.1455, "longitude": -151.1904, "radius_km": 400}
    )
    utm = gdf.to_crs(gdf.estimate_utm_crs())
    assert utm.area.iloc[0] == pytest.approx(np.pi * 400e3**2, rel=0.01)


@pytest.mark.parametrize(
    "spec",
    [
        {"bbox": [10, 0, 5, 1]},
        {"latitude": 95, "longitude": 0, "radius_km": 1},
        {"intersects": {"type": "Blob"}},
    ],
)
def test_bad_inline_aois(spec, tmp_path):
    with pytest.raises(ValueError):
        stagefile.parse_aoi(spec, tmp_path)


def test_relative_aoi_paths_resolve_against_the_file(tmp_path):
    (tmp_path / "a.geojson").write_text("{}")
    cfg = stagefile.load_config(
        write(tmp_path, "name: x\naoi: a.geojson\ndatasets: {nlcd: {}}\n")
    )
    assert cfg.requests[0].aoi == str(tmp_path / "a.geojson")


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
    with pytest.raises(SystemExit, match="can't be combined with -s 2023-09-14"):
        stage_app(["--config", write(tmp_path, PETERS), "-s", "2023-09-14"])


def test_dry_run_prints_resolved_config(tmp_path, capsys):
    stagefile.stage(config=write(tmp_path, PETERS), dry_run=True)
    out = capsys.readouterr().out
    assert "window: 2023-09-13T20:06:49/2023-09-15T02:06:49" in out
    assert "variable: temperature_2m" in out
    assert "aoi_bbox: [-152.2, 62.75, -150.2, 63.55]" in out


def test_run_writes_one_group_per_dataset(tmp_path, monkeypatch):
    calls = {}

    def fake_prism(aoi, start_date, end_date, variables: tuple[str, ...] = ("ppt",)):
        calls["prism"] = (start_date, end_date, variables)
        return xr.Dataset({"precipitation": ("time", [1.0])}, attrs={"title": "PRISM"})

    monkeypatch.setitem(stagefile.DATASETS, "prism", (fake_prism, "rasters"))
    text = (
        PETERS.split("datasets:")[0] + "datasets:\n  prism: {variables: [ppt, tmax]}\n"
    )
    out = tmp_path / "out.zarr"
    stagefile.stage(config=write(tmp_path, text), output_path=str(out))

    # The window's exclusive end becomes load()'s inclusive end
    assert calls["prism"] == (
        "2023-06-01T00:00:00",
        "2023-09-21T23:59:59",
        ("ppt", "tmax"),
    )
    tree = xr.open_datatree(out, engine="zarr")
    assert tree.attrs["title"] == "Peters Dome"
    assert tree["rasters/prism"].attrs["title"] == "PRISM"
