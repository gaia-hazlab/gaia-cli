import json

import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import yaml

from gaia_cli import align, catalog, io, runner, stagefile, waveforms
from gaia_cli.aoi import load_aoi
from gaia_cli.cli import stage_app
from gaia_cli.datasets import seismic

T0 = pd.Timestamp("2023-09-14T02:06:00")
ORIGIN = {"latitude": 63.1455, "longitude": -151.1904}
STAGE = """
name: pd
time: 2023-09-14T02:06:49Z
origin: {latitude: 63.1455, longitude: -151.1904}
aoi: {latitude: 63.5, longitude: -151.0, radius_km: 100}
datasets:
  seismic: {window: {before: 49s, after: 11s}}
"""


def add_stations(fdsn, seconds=60, gap=True):
    """BPAW (BH and HH, 50 and 100 Hz), and BB (BH only, at a jittered rate)."""
    for c in "ZNE":
        fdsn.channel(f"AK.BPAW..BH{c}", lat=63.8, lon=-151.0)
        fdsn.channel(f"AK.BPAW..HH{c}", lat=63.8, lon=-151.0, rate=100.0)
        fdsn.channel(f"AK.BB..BH{c}", lat=63.2, lon=-150.5, units="M/S")
        fdsn.trace(f"AK.BPAW..BH{c}", 50.0, T0, np.arange(50 * seconds))
        fdsn.trace(f"AK.BPAW..HH{c}", 100.0, T0, np.arange(100 * seconds))
    for c in "ZN":
        fdsn.trace(
            f"AK.BB..BH{c}", 50.0, T0 + pd.Timedelta("8.4ms"), np.ones(50 * seconds)
        )
    if gap:  # BHE has a 10 s gap
        fdsn.trace("AK.BB..BHE", 50.0, T0, np.ones(500))
        fdsn.trace("AK.BB..BHE", 50.0, T0 + pd.Timedelta("20s"), np.ones(50 * 40))


def test_load_returns_one_child_per_rate(fdsn):
    add_stations(fdsn)
    aoi = load_aoi({"latitude": 63.5, "longitude": -151.0, "radius_km": 100})
    end = T0 + pd.Timedelta("60s") - pd.Timedelta("1us")
    tree = seismic.load(aoi, T0.isoformat(), end.isoformat())

    assert set(tree.children) == {"50hz", "100hz"}
    g50, g100 = (tree[k].to_dataset() for k in ("50hz", "100hz"))
    # HH preferred at BPAW, so its BH channels are neither requested nor kept
    assert g100.channel.values.tolist() == [f"AK.BPAW..HH{c}" for c in "ENZ"]
    assert g50.channel.values.tolist() == [f"AK.BB..BH{c}" for c in "ENZ"]
    assert not any("BPAW -- BH" in p for p in fdsn.posts)
    assert g50.sizes["time"] == 3000 and g100.sizes["time"] == 6000
    assert g50.sensitivity_units.values.tolist() == ["m/s"] * 3
    assert g50.time_offset.values[1] == np.timedelta64(8_400_000, "ns")
    assert (g50.counts.values[0, 500:1000] == waveforms.FILL).all()
    # Without a stage file origin, distances are from the AOI centroid
    assert g50.attrs["gaia:origin"]["source"] == "AOI centroid"
    assert g50.attrs["gaia:origin"]["latitude"] == pytest.approx(63.5, abs=0.01)
    assert g50.attrs["networks"] == ["AK"]
    assert g50.attrs["start_datetime"] == "2023-09-14T02:06:00Z"
    assert g50.attrs["gaia:query"]["channel"] == "HH?,BH?"
    # A second's margin around the window
    assert "2023-09-14T02:05:59.000000 2023-09-14T02:07:00.999999" in fdsn.posts[0]

    tree = seismic.load(aoi, T0.isoformat(), end.isoformat(), origin=ORIGIN)
    attrs = tree["50hz"].attrs["gaia:origin"]
    assert attrs == {**ORIGIN, "source": "stage file origin"}


def test_load_errors(fdsn):
    aoi = load_aoi({"latitude": 63.5, "longitude": -151.0, "radius_km": 100})
    with pytest.raises(ValueError, match="no HH.,BH. channels"):
        seismic.load(aoi, "2023-09-14T02:06", "2023-09-14T02:07")
    fdsn.channel("AK.BB..BHZ", lat=63.2, lon=-150.5)
    with pytest.raises(RuntimeError, match="no waveform data"):
        seismic.load(aoi, "2023-09-14T02:06", "2023-09-14T02:07")
    with pytest.raises(ValueError, match="is not a time"):
        seismic._time("NaT")


def test_stage_writes_waveform_groups(fdsn, tmp_path, monkeypatch, capsys):
    add_stations(fdsn)
    config = tmp_path / "stage.yaml"
    config.write_text(STAGE)
    runner.stage(config=str(config), dry_run=True)
    entry = yaml.safe_load(capsys.readouterr().out)["datasets"]["seismic"]
    assert entry["group"] == "waveforms/seismic"
    # Read exactly the window (PAD = 0)
    assert entry["read"] == entry["window"] == "2023-09-14T02:06:00/2023-09-14T02:07:00"
    assert "origin" not in entry

    # Large groups get one-channel chunks
    monkeypatch.setattr(io, "_CHUNK_TARGET_BYTES", 10_000)
    monkeypatch.setattr(waveforms, "CHUNK_SECONDS", 20)
    out = tmp_path / "out"
    runner.stage(config=str(config), output_path=str(out))

    tree = io.open_datatree(str(out / "pd.zarr"))
    assert set(tree["waveforms/seismic"].children) == {"50hz", "100hz"}
    stored = tree["waveforms/seismic/50hz"].to_dataset()
    assert stored.counts.encoding["chunks"] == (1, 1000)
    assert stored.counts.encoding["_FillValue"] == waveforms.FILL
    g = waveforms.decoded(stored)
    assert g.xindexes["geometry"].crs.to_epsg() == 4326
    # The gap reads back as NaN, data as counts
    assert g.counts.isel(channel=0, time=slice(500, 1000)).isnull().all()
    assert float(g.counts.isel(channel=1, time=0)) == 1.0
    assert g.location.values.tolist() == ["", "", ""]
    assert g.time_offset.values[1] == np.timedelta64(8_400_000, "ns")
    assert g.attrs["gaia:origin"] == {**ORIGIN, "source": "stage file origin"}
    # Distances from the stage file's origin
    assert g.distance_km.values[0] == pytest.approx(35.3, abs=0.1)

    item = json.loads((out / "pd.json").read_text())
    assert item["bbox"] == pytest.approx([-151.0, 63.2, -150.5, 63.8])
    assert item["properties"]["start_datetime"] == "2023-09-14T02:06:00Z"
    assert "EarthScope Consortium" in [
        p["name"] for p in item["properties"]["providers"]
    ]

    # The store is found by the catalog, and an existing group is skipped
    catalog.rebuild(str(tmp_path))
    links = json.loads((tmp_path / "catalog.json").read_text())["links"]
    assert "./out/pd.json" in [lk["href"] for lk in links]
    posts = len(fdsn.posts)
    runner.stage(config=str(config), output_path=str(out), skip_existing=True)
    assert len(fdsn.posts) == posts


def test_align_passes_waveform_groups_through(fdsn, tmp_path):
    add_stations(fdsn, gap=False)
    config = tmp_path / "stage.yaml"
    config.write_text(STAGE)
    runner.stage(config=str(config), output_path=str(tmp_path / "out"))
    tree = io.open_datatree(str(tmp_path / "out" / "pd.zarr"))
    raster = xr.DataArray(
        np.ones((2, 2)),
        coords={"y": [63.5, 63.0], "x": [-151.0, -150.5]},
        dims=("y", "x"),
        name="elevation",
    ).rio.write_crs("EPSG:4326")
    nodes = {p: n.to_dataset() for p, n in tree.subtree_with_keys}
    nodes["dem/x"] = raster.to_dataset()
    aligned = align.align(xr.DataTree.from_dict(nodes))
    before = tree["waveforms/seismic/50hz"].to_dataset()
    after = aligned["waveforms/seismic/50hz"].to_dataset()
    xr.testing.assert_identical(before, after)
    io.save_datatree(aligned, str(tmp_path / "aligned.zarr"))


def test_seismic_command_takes_fdsn_codes_but_not_origin(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    load_aoi({"latitude": 64.0997, "longitude": -150.9873, "radius_km": 5}).to_file(aoi)
    args = [
        "seismic",
        "-i",
        str(aoi),
        "-s",
        "2023-09-14T02:06",
        "-e",
        "2023-09-14T02:08",
    ]
    stage_app(
        [*args, "-n", "AK", "--station", "BPAW", "-c", "BH?", "-o", "x", "--dry-run"],
        result_action="return_value",
    )
    entry = yaml.safe_load(capsys.readouterr().out)["datasets"]["seismic"]
    assert {k: entry[k] for k in ("network", "station", "location", "channel")} == {
        "network": "AK",
        "station": "BPAW",
        "location": "*",
        "channel": "BH?",
    }
    # A time-only end is that instant, not the whole day
    assert entry["window"] == "2023-09-14T02:06:00/2023-09-14T02:08:00"
    assert "origin" not in stagefile.options(seismic.load)
    with pytest.raises(SystemExit):
        stage_app([*args, "--origin", "x", "-o", "x", "--dry-run"])
