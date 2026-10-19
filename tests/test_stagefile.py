import json

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from gaia_cli import stagefile
from gaia_cli.cli import stage_app
from gaia_cli.hrrr import HRRRVariable

T = pd.Timestamp("2023-09-14T02:06:49")


def write(tmp_path, text, name="stage.yaml"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def problems(tmp_path, text):
    with pytest.raises(stagefile.StageFileError) as e:
        stagefile.load_config(write(tmp_path, text))
    return e.value.problems


# ---- window grammar ------------------------------------------------------


def test_interval_string_date_only_end_includes_the_day():
    (w,) = stagefile.parse_window("2025-12-01/2025-12-20", None)
    assert w.start == pd.Timestamp("2025-12-01")
    assert w.end == pd.Timestamp("2025-12-21")


def test_interval_string_open_end():
    (w,) = stagefile.parse_window("2025-12-01/..", None)
    assert w.end is None


def test_start_end_mapping_converts_offsets_to_utc():
    (w,) = stagefile.parse_window(
        {"start": "2025-12-01T00:00-08:00", "end": "2025-12-02T06:00:00Z"}, None
    )
    assert w.start == pd.Timestamp("2025-12-01T08:00")
    assert w.end == pd.Timestamp("2025-12-02T06:00")


def test_yaml_dates_and_datetimes(tmp_path):
    import datetime as dt

    (w,) = stagefile.parse_window(
        {"start": dt.date(2025, 12, 1), "end": dt.date(2025, 12, 20)}, None
    )
    assert (w.start, w.end) == (pd.Timestamp("2025-12-01"), pd.Timestamp("2025-12-21"))


def test_relative_window():
    (w,) = stagefile.parse_window({"before": "17min", "after": "23min"}, T)
    assert w.start == T - pd.Timedelta("17min")
    assert w.end == T + pd.Timedelta("23min")


def test_relative_one_side_defaults_to_time():
    (w,) = stagefile.parse_window({"before": "30d"}, T)
    assert (w.start, w.end) == (T - pd.Timedelta(days=30), T)


def test_mixed_absolute_and_relative():
    (w,) = stagefile.parse_window({"start": "2023-06-01", "after": "1w"}, T)
    assert (w.start, w.end) == (pd.Timestamp("2023-06-01"), T + pd.Timedelta(days=7))


def test_pre_post_windows_with_select():
    pre, post = stagefile.parse_window(
        {
            "pre": {"before": "30d", "select": "last"},
            "post": {"after": "30d", "select": "first", "count": 2},
        },
        T,
    )
    assert (pre.name, pre.end, pre.select) == ("pre", T, "last")
    assert (post.name, post.start, post.select, post.count) == ("post", T, "first", 2)


@pytest.mark.parametrize(
    "spec, message",
    [
        ({"before": "17min"}, "need a top-level 'time'"),
        ("2025-12-01", "not an interval"),
        ("2025-12-20/2025-12-01", "starts at"),
        ({"start": "2025-12-01", "before": "1d"}, "either start or before"),
        ({"strat": "2025-12-01"}, "did you mean 'start'"),
        ("../..", "neither a start nor an end"),
        ({"start": "2025-12-01", "end": "2025-12-02", "select": "best"}, "select must"),
        (
            {"start": "2025-12-01", "end": "2025-12-02", "count": 2},
            "count needs select",
        ),
    ],
)
def test_bad_windows(spec, message):
    with pytest.raises(ValueError, match=message):
        stagefile.parse_window(spec, None)


@pytest.mark.parametrize("text", ["30d", "17min", "6h", "2w", "90s", "1.5h", "PT17M"])
def test_durations(text):
    assert stagefile.parse_duration(text) > pd.Timedelta(0)


@pytest.mark.parametrize("text", ["3M", "1y", "30m", "-1d", "soon"])
def test_bad_durations(text):
    with pytest.raises(ValueError):
        stagefile.parse_duration(text)


# ---- AOIs ------------------------------------------------------------------


def test_bbox_aoi():
    gdf = stagefile.load_aoi({"bbox": [-152.2, 62.75, -150.2, 63.55]})
    np.testing.assert_allclose(gdf.total_bounds, [-152.2, 62.75, -150.2, 63.55])


def test_circle_aoi_has_true_radius():
    gdf = stagefile.load_aoi(
        {"latitude": 63.1455, "longitude": -151.1904, "radius_km": 400}
    )
    utm = gdf.to_crs(gdf.estimate_utm_crs())
    assert utm.area.iloc[0] == pytest.approx(np.pi * 400e3**2, rel=0.01)


def test_intersects_aoi():
    geom = {"type": "Point", "coordinates": [-151.2, 63.1]}
    assert stagefile.load_aoi({"intersects": geom}).geometry.iloc[0].x == -151.2


@pytest.mark.parametrize(
    "spec",
    [
        {"bbox": [1, 2, 3]},
        {"bbox": [10, 0, 5, 1]},
        {"latitude": 95, "longitude": 0, "radius_km": 1},
        {"latitude": 0, "longitude": 0},
        {"intersects": {"type": "Blob"}},
    ],
)
def test_bad_inline_aois(spec, tmp_path):
    with pytest.raises(ValueError):
        stagefile.parse_aoi(spec, tmp_path)


def test_relative_aoi_paths_resolve_against_the_file(tmp_path):
    (tmp_path / "aoi").mkdir()
    (tmp_path / "aoi" / "a.geojson").write_text("{}")
    path = write(tmp_path, "name: x\naoi: aoi/a.geojson\ndatasets: {nlcd: {}}\n")
    cfg = stagefile.load_config(path)
    assert cfg.requests[0].aoi == str(tmp_path / "aoi" / "a.geojson")


# ---- whole files -----------------------------------------------------------

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


def test_load_config_resolves_defaults_named_aois_and_options(tmp_path):
    cfg = stagefile.load_config(write(tmp_path, PETERS))
    reqs = {r.name: r for r in cfg.requests}
    assert cfg.time == T
    assert cfg.output == "s3://gaia/scratch/peters-dome"
    assert reqs["hrrr"].aoi == {
        "latitude": 63.1455,
        "longitude": -151.1904,
        "radius_km": 15,
    }
    assert reqs["hrrr"].windows[0].start == T - pd.Timedelta("6h")
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


def test_timed_dataset_needs_a_window(tmp_path):
    text = "name: x\naoi: {bbox: [0, 0, 1, 1]}\ndatasets: {prism: {}}\n"
    assert any("no window" in p for p in problems(tmp_path, text))


def test_unreadable_file(tmp_path):
    with pytest.raises(stagefile.StageFileError):
        stagefile.load_config(str(tmp_path / "nope.yaml"))


# ---- CLI and runner ----------------------------------------------------------


def test_stage_config_is_the_default_command():
    assert stage_app.default_command is stagefile.stage


def test_flags_and_file_are_exclusive(tmp_path):
    path = write(tmp_path, PETERS)
    with pytest.raises(SystemExit):
        stage_app(["--config", path, "-s", "2023-09-14"], exit_on_error=True)


def test_dry_run_prints_resolved_config(tmp_path, capsys):
    stagefile.stage(config=write(tmp_path, PETERS), dry_run=True)
    out = capsys.readouterr().out
    assert "group: rasters/hrrr" in out
    assert "window: 2023-09-13T20:06:49/2023-09-15T02:06:49" in out
    assert "variable: temperature_2m" in out
    assert "aoi_bbox: [-152.2, 62.75, -150.2, 63.55]" in out


def test_dry_run_exits_with_problems(tmp_path, capsys):
    with pytest.raises(SystemExit, match="did you mean 'prism'"):
        stagefile.stage(
            config=write(tmp_path, "name: x\ndatasets: {prsm: {}}\n"), dry_run=True
        )


def test_run_writes_one_group_per_dataset(tmp_path, monkeypatch):
    calls = {}

    def fake_load(
        aoi, start_date, end_date, variable=HRRRVariable.precipitation_surface
    ):
        calls["hrrr"] = (start_date, end_date, variable)
        time = pd.date_range(start_date, periods=3, freq="h")
        return xr.Dataset({"t2m": ("time", np.arange(3.0))}, coords={"time": time})

    def fake_prism(aoi, start_date, end_date, variables: tuple[str, ...] = ("ppt",)):
        calls["prism"] = (start_date, end_date, variables)
        return xr.Dataset({"precipitation": ("time", [1.0])}, attrs={"title": "PRISM"})

    monkeypatch.setitem(stagefile.DATASETS, "hrrr", (fake_load, "rasters"))
    monkeypatch.setitem(stagefile.DATASETS, "prism", (fake_prism, "rasters"))
    text = PETERS.replace("  nlcd: {year: 2021}\n", "")
    out = tmp_path / "out.zarr"
    stagefile.stage(config=write(tmp_path, text), output_path=str(out))

    # Exclusive window ends become the loaders' inclusive end
    assert calls["prism"] == (
        "2023-06-01T00:00:00",
        "2023-09-21T23:59:59",
        ("ppt", "tmax"),
    )
    assert calls["hrrr"][:2] == ("2023-09-13T20:06:49", "2023-09-15T02:06:48")
    tree = xr.open_datatree(out, engine="zarr")
    assert set(tree.groups) >= {"/rasters/hrrr", "/rasters/prism"}
    assert tree.attrs["title"] == "Peters Dome"
    assert tree["rasters/prism"].attrs["title"] == "PRISM"
    assert json.loads((tmp_path / "out.json").read_text())["type"] == "Feature"


def test_named_and_plain_window_keys_dont_mix():
    with pytest.raises(ValueError, match="unknown window key 'pre'"):
        stagefile.parse_window({"pre": {"before": "1d"}, "start": "2023-09-01"}, T)
