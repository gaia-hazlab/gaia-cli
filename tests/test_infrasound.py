import numpy as np
import pandas as pd
import yaml

from gaia_cli import io, runner, waveforms
from gaia_cli.aoi import load_aoi
from gaia_cli.datasets import infrasound

T0 = pd.Timestamp("2023-09-14T02:06:00")
STAGE = """
name: pd
time: 2023-09-14T02:06:49Z
origin: {latitude: 63.1455, longitude: -151.1904}
aoi: {latitude: 63.5, longitude: -151.0, radius_km: 100}
datasets:
  infrasound: {window: {before: 49s, after: 11s}}
"""


def add_sensors(fdsn, seconds=60):
    """An array of three elements by location (100 Hz), and a sensor with a
    20 Hz and a 40 Hz stream, only the 40 Hz one with data."""
    for i in range(1, 4):
        sid = f"AV.ARR.0{i}.HDF"
        fdsn.channel(
            sid, lat=63.2 + i * 1e-3, lon=-150.5, rate=100.0, scale=11000.0, units="PA"
        )
        fdsn.trace(sid, 100.0, T0, np.full(100 * seconds, 11000 * i))
    fdsn.channel("AK.CUT.20.BDF", lat=63.8, rate=20.0, scale=64704.1, units="Pa")
    fdsn.channel("AK.CUT.EP.BDF", lat=63.8, rate=40.0, scale=64704.1, units="Pa")
    fdsn.trace("AK.CUT.EP.BDF", 40.0, T0, np.full(40 * seconds, 64704))


def test_load_keeps_array_elements_and_the_fastest_stream(fdsn):
    add_sensors(fdsn)
    aoi = load_aoi({"latitude": 63.5, "longitude": -151.0, "radius_km": 100})
    end = T0 + pd.Timedelta("60s") - pd.Timedelta("1us")
    tree = infrasound.load(aoi, T0.isoformat(), end.isoformat())

    assert set(tree.children) == {"sites", "channels", "40hz", "100hz"}
    assert tree["100hz"].channel.values.tolist() == [
        f"AV.ARR.0{i}.HDF" for i in (1, 2, 3)
    ]
    assert tree["40hz"].channel.values.tolist() == ["AK.CUT.EP.BDF"]
    assert not any("CUT 20" in p for p in fdsn.posts)
    assert tree["sites"].site.values.tolist() == [
        "AK.CUT.EP",
        "AV.ARR.01",
        "AV.ARR.02",
        "AV.ARR.03",
    ]
    assert tree.attrs["networks"] == ["AK", "AV"]
    assert tree.attrs["gaia:query"]["channel"] == "BDF,HDF"
    assert tree.attrs["title"].startswith("Infrasound")
    pressure = waveforms.to_physical(tree, "100hz")
    assert pressure.name == "pressure" and pressure.attrs["units"] == "Pa"
    np.testing.assert_allclose(pressure.isel(time=0).values, [1.0, 2.0, 3.0])


def test_stage_writes_waveforms_infrasound(fdsn, tmp_path, capsys):
    add_sensors(fdsn)
    config = tmp_path / "stage.yaml"
    config.write_text(STAGE)
    runner.stage(config=str(config), dry_run=True)
    entry = yaml.safe_load(capsys.readouterr().out)["datasets"]["infrasound"]
    assert entry["group"] == "waveforms/infrasound"
    assert entry["read"] == entry["window"] == "2023-09-14T02:06:00/2023-09-14T02:07:00"
    assert entry["channel"] == "BDF,HDF"

    out = tmp_path / "out"
    runner.stage(config=str(config), output_path=str(out))
    node = io.open_datatree(str(out / "pd.zarr"))["waveforms/infrasound"]
    grid = waveforms.to_grid(node, rate=20)
    assert grid.pressure.attrs["units"] == "Pa"
    assert grid.sizes == {"channel": 4, "time": 1200}
    np.testing.assert_allclose(
        grid.pressure.isel(time=600), [1.0, 1.0, 2.0, 3.0], rtol=1e-4
    )
    assert node.attrs["gaia:origin"]["source"] == "stage file origin"
