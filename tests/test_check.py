import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from shapely.geometry import Point

from gaia_cli import check, io, regrid, stations
from gaia_cli.cli import app

SOURCE = {"title": "t", "license": "CC0-1.0", "providers": [], "links": []}


def _daily(**attrs):
    times = pd.date_range("2025-12-01", periods=3, freq="D")
    da = xr.DataArray(
        np.ones((3, 2, 3), "float32"),
        coords={
            "time": times,
            "y": [5e6 + 15, 5e6 + 5],
            "x": [5e5 + 5, 5e5 + 15, 5e5 + 25],
        },
        dims=("time", "y", "x"),
        name="ppt",
        attrs={"units": "mm", regrid.METHOD_ATTR: "conservative"},
    ).rio.write_crs("EPSG:32610")
    return io.with_time_bounds(da, -12, 12, "time: sum").assign_attrs(SOURCE | attrs)


def _stations():
    sites = xr.Dataset(
        coords={
            "station": ["A", "B"],
            "geometry": ("station", [Point(0, 0), Point(1, 1)]),
        }
    )
    reports = [
        pd.DataFrame(
            {"tmp": [1.0, 2.0]}, index=pd.date_range("2025-12-01", periods=2, freq="h")
        ).rename_axis("time"),
        pd.DataFrame(
            {"tmp": [3.0]}, index=pd.DatetimeIndex(["2025-12-01"], name="time")
        ),
    ]
    ds = stations.ragged(sites, reports).assign_attrs(SOURCE)
    ds["tmp"].attrs["units"] = "degC"
    return ds


def _findings(tree):
    return [f for _, _, fs in check.check(tree) for f in fs]


def _messages(groups):
    return [f.message for f in _findings(xr.DataTree.from_dict(groups))]


def test_staged_store_has_no_findings(tmp_path, capsys):
    tree = xr.DataTree.from_dict(
        {"rasters/prism": _daily(), "stations/synoptic": _stations()}
    )
    store = str(tmp_path / "s.zarr")
    io.save_datatree(tree, store)

    app(["check", store], result_action="return_value")
    out = capsys.readouterr().out
    assert "rasters/prism: raster, 2 x 3 cells of 10m x 10m, EPSG:32610" in out
    assert "every 1D, each covering (-12h, 12h] around its time" in out
    assert "stations/synoptic: stations, 2 station(s), 3 report(s)" in out
    assert out.rstrip().endswith("0 error(s), 0 warning(s)")


def test_period_without_cell_methods_is_an_error():
    ds = _daily()
    del ds["ppt"].attrs["cell_methods"]
    assert _messages({"rasters/prism": ds}) == [
        "ppt has time bounds but no cell_methods: gaia align averages it (time: mean)"
    ]


def test_instantaneous_values_need_no_cell_methods():
    ds = _daily().drop_vars("time_bounds")
    del ds.time.attrs["bounds"]
    del ds["ppt"].attrs["cell_methods"]
    assert _messages({"rasters/hrrr": ds}) == []


def test_missing_bounds_variable_and_unknown_methods():
    ds = _daily().drop_vars("time_bounds")
    ds["ppt"].attrs |= {"cell_methods": "time: median", regrid.METHOD_ATTR: "spline"}
    messages = _messages({"rasters/prism": ds})
    assert any("names bounds 'time_bounds', which is missing" in m for m in messages)
    assert any("'time: median'" in m for m in messages)
    assert any("'spline' is not one of" in m for m in messages)


def test_grid_problems():
    no_crs = _daily().drop_vars("spatial_ref")
    irregular = _daily().assign_coords(x=[5e5 + 5, 5e5 + 15, 5e5 + 50])
    assert _messages({"a": no_crs})[0].startswith("no spatial_ref coordinate")
    assert _messages({"b": irregular}) == [
        "x spacing is irregular; gaia align expects a regular grid"
    ]


def test_metadata_warnings():
    ds = _daily()
    ds.attrs = {}
    del ds["ppt"].attrs["units"]
    findings = _findings(xr.DataTree.from_dict({"rasters/prism": ds}))
    assert {f.level for f in findings} == {check.WARNING}
    assert [f.message for f in findings] == [
        "ppt has no units",
        "no title, license, providers, links (source provenance)",
    ]


def test_ragged_array_must_account_for_every_report():
    ds = _stations()
    ds["row_size"] = ds.row_size + 1
    assert _messages({"stations/x": ds}) == [
        "row_size adds up to 5, but obs has 3 reports"
    ]


def test_store_not_written_by_gaia(tmp_path, capsys):
    # e.g. written by hand, with integer times and no CRS
    ds = xr.Dataset(
        {"v": (("time", "y", "x"), np.zeros((2, 2, 2)))},
        coords={"time": [0, 1], "y": [1.0, 0.0], "x": [0.0, 1.0]},
    )
    store = str(tmp_path / "other.zarr")
    xr.DataTree.from_dict({"grid": ds}).to_zarr(store, zarr_format=3)

    with pytest.raises(SystemExit) as e:
        app(["check", "-i", store], result_action="return_value")
    assert e.value.code == 1
    out = capsys.readouterr().out
    assert "no gaia:cli-version: not written by gaia stage or gaia align" in out
    assert "ERROR: time is int64, not datetime64" in out
    assert "ERROR: no spatial_ref coordinate: gaia align skips this group" in out


def test_strict_fails_on_warnings(tmp_path):
    ds = _daily()
    del ds["ppt"].attrs["units"]
    store = str(tmp_path / "s.zarr")
    io.save_datatree(xr.DataTree.from_dict({"rasters/prism": ds}), store)

    app(["check", store], result_action="return_value")  # warnings only: exit 0
    with pytest.raises(SystemExit):
        app(["check", store, "--strict"], result_action="return_value")


def test_unreadable_store(tmp_path):
    with pytest.raises(SystemExit, match="could not open"):
        app(["check", str(tmp_path / "missing.zarr")], result_action="return_value")


def _waveforms(times):
    rate = xr.Dataset(
        {"counts": (("channel", "time"), np.zeros((2, len(times)), "int32"))},
        coords={"channel": ["AK.A..HHZ", "AK.A..HHN"], "time": times},
    )
    return xr.DataTree.from_dict(
        {
            "waveforms/seismic": xr.Dataset(attrs=SOURCE),
            "waveforms/seismic/sites": xr.Dataset({"site": ("geometry", ["AK.A."])}),
            "waveforms/seismic/channels": xr.Dataset({"cha": ("channel", ["Z", "N"])}),
            "waveforms/seismic/50hz": rate,
        }
    )


def test_waveforms_are_described_with_their_children():
    times = pd.date_range("2023-09-14", periods=4, freq="20ms")
    (path, summary, findings), *rest = check.check(_waveforms(times))
    assert path == "waveforms/seismic" and findings == [] and rest == []
    assert summary == [
        "waveforms/seismic: waveforms, 1 site(s), 2 channel(s); 50hz (2 channel(s))"
    ]
    assert [(f.group, f.message) for f in _findings(_waveforms(times[::-1]))] == [
        ("waveforms/seismic/50hz", "time is not increasing")
    ]
