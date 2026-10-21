import logging

import numpy as np
import pandas as pd
import pytest
import requests
import xarray as xr
from conftest import STATION_HEADER, mseed

from gaia_cli import io, waveforms
from gaia_cli.aoi import load_aoi

T0 = pd.Timestamp("2023-09-14T02:06:00")
S = pd.Timedelta(seconds=1)


def table(fdsn):
    return waveforms.parse_channel_text("\n".join([STATION_HEADER, *fdsn.rows]))


# ---- channels ----------------------------------------------------------------


def test_parse_channel_text_keeps_empty_locations_and_normalizes_units(fdsn):
    fdsn.channel("AK.BPAW..BHZ", units="M/S")
    fdsn.channel("AK.BPAW.00.HNZ", units="M/S**2", end="2024-01-01T00:00:00")
    fdsn.channel("AV.ACH.EP.BDF", units="PA")
    fdsn.channel("AV.ACH.EP.BDF", units="Pa", azimuth="")
    t = table(fdsn)
    assert t.location.tolist() == ["", "00", "EP", "EP"]
    assert t.sensitivity_units.tolist() == ["m/s", "m/s**2", "Pa", "Pa"]
    assert t.sample_rate.tolist() == [50.0] * 4
    assert t.end.isna().tolist() == [True, False, True, True]
    assert t.azimuth.isna().tolist() == [False, False, False, True]
    assert waveforms.seed_ids(t).tolist()[0] == "AK.BPAW..BHZ"
    assert waveforms.parse_channel_text("").empty


def test_select_epochs_keeps_the_epoch_covering_the_start(fdsn, caplog):
    fdsn.channel("AK.A..BHZ", end="2020-01-01T00:00:00")  # ended before
    fdsn.channel("AK.A..BHZ", start="2020-01-01T00:00:00", end="2023-09-14T02:07:00")
    fdsn.channel("AK.A..BHZ", start="2023-09-14T02:07:00", scale=5e8)
    fdsn.channel("AK.B..BHZ", start="2023-09-14T02:07:00")  # starts in the window
    caplog.set_level(logging.WARNING)
    out = waveforms.select_epochs(table(fdsn), T0, T0 + 120 * S)
    assert waveforms.seed_ids(out).tolist() == ["AK.A..BHZ", "AK.B..BHZ"]
    assert out.end.iloc[0] == pd.Timestamp("2023-09-14T02:07:00")
    assert "AK.A..BHZ: sensitivity changes in the window" in caplog.text


def test_prefer_bands_keeps_one_band_and_location_per_station(fdsn):
    for c in "ZNE":
        fdsn.channel(f"AK.BOTH..BH{c}")
        fdsn.channel(f"AK.BOTH..HH{c}")
        fdsn.channel(f"AK.BB..BH{c}")
        fdsn.channel(f"AK.LOCS.10.HH{c}")
        fdsn.channel(f"AK.LOCS.00.HH{c}")
    out = waveforms.prefer_bands(table(fdsn), "HH?,BH?")
    ids = sorted(waveforms.seed_ids(out))
    assert ids == sorted(
        [f"AK.BOTH..HH{c}" for c in "ZNE"]
        + [f"AK.BB..BH{c}" for c in "ZNE"]
        + [f"AK.LOCS.00.HH{c}" for c in "ZNE"]
    )


def test_query_channels_keeps_stations_inside_the_aoi_polygon(fdsn):
    aoi = load_aoi({"latitude": 64.0, "longitude": -151.0, "radius_km": 50})
    fdsn.channel("AK.IN..BHZ", lat=64.0, lon=-151.0)
    # Inside the circle's bounding box, outside the circle
    fdsn.channel("AK.CORNER..BHZ", lat=64.4, lon=-150.1)
    out = waveforms.query_channels(aoi, T0, T0 + 60 * S, channel="BH?", location="")
    assert waveforms.seed_ids(out).tolist() == ["AK.IN..BHZ"]
    (params,) = fdsn.gets
    assert params["location"] == "--"
    assert params["starttime"] == "2023-09-14T02:06:00.000000"
    assert params["minlatitude"] == pytest.approx(64.0 - 50 / 111.2, abs=0.01)
    fdsn.rows.clear()  # no channels: the service answers 404
    assert waveforms.query_channels(aoi, T0, T0 + 60 * S).empty


def test_check_size_refuses_large_requests(fdsn):
    for i in range(3):
        fdsn.channel(f"AK.S{i}..HHZ", rate=100.0)
    t = table(fdsn)
    # 3 channels x 100 Hz x 1 day x 4 bytes
    assert waveforms.estimated_bytes(t, T0, T0 + pd.Timedelta(days=1)) == 103_680_000
    waveforms.check_size(t, T0, T0 + pd.Timedelta(days=1))
    with pytest.raises(ValueError, match="3 channels .* over the 0.1 GB limit"):
        waveforms.check_size(t, T0, T0 + pd.Timedelta(days=1), limit=1e8)


# ---- dataselect --------------------------------------------------------------


def test_dataselect_posts_in_chunks_and_skips_failed_ones(fdsn, monkeypatch, caplog):
    monkeypatch.setattr(waveforms, "STATIONS_PER_REQUEST", 2)
    for i in range(6):
        fdsn.channel(f"AK.S{i}..BHZ")
        fdsn.trace(f"AK.S{i}..BHZ", 50.0, T0, np.arange(50))
    # Chunk 1 succeeds after a retry, chunk 2 has no data, chunk 3 always fails
    fdsn.statuses = [503, 200, 204, 500, 500, 500]
    caplog.set_level(logging.INFO)
    buffers = waveforms.dataselect(table(fdsn), T0, T0 + S)
    assert len(fdsn.posts) == 6
    assert fdsn.posts[0].splitlines() == [
        "AK S0 -- BHZ 2023-09-14T02:06:00.000000 2023-09-14T02:06:01.000000",
        "AK S1 -- BHZ 2023-09-14T02:06:00.000000 2023-09-14T02:06:01.000000",
    ]
    assert [len(b) > 0 for b in buffers] == [True, False]
    assert sorted(waveforms.decode(buffers)) == ["AK.S0..BHZ", "AK.S1..BHZ"]
    assert "request 2: no data" in caplog.text
    assert "request 3 (2 stations) failed, skipped: 500" in caplog.text


def test_post_retries_timeouts_and_raises_client_errors(fdsn, monkeypatch):
    calls = []

    def post(url, data, timeout):
        calls.append(url)
        if len(calls) == 1:
            raise requests.Timeout("slow")
        return fdsn.post(url, data, timeout)

    monkeypatch.setattr(waveforms.requests, "post", post)
    assert waveforms._post("u", "AK S0 -- BHZ a b\n") == b""
    assert len(calls) == 2
    fdsn.statuses = [400]
    with pytest.raises(requests.HTTPError):
        waveforms._post("u", "")
    assert len(calls) == 3  # no retry


# ---- decoding and the grid -----------------------------------------------------


def test_decode_skips_float_samples(caplog):
    ints = mseed(("AK.A..BHZ", 50.0, T0, np.arange(10)))
    floats = mseed(("AK.A..BHN", 50.0, T0, np.arange(10)), sample_type="f")
    segments = waveforms.decode([ints, floats, b""])
    assert list(segments) == ["AK.A..BHZ"]
    (seg,) = segments["AK.A..BHZ"]
    assert (seg.rate, seg.start, seg.data.dtype) == (50.0, T0.value, np.int32)
    assert "AK.A..BHN: f samples skipped" in caplog.text


@pytest.mark.parametrize(
    "rate, nominal, name",
    [(50.0, 50.0, "50hz"), (49.99999, 50.0, "50hz"), (0.1, 0.1, "0.1hz")],
)
def test_nominal_rate_and_group_name(rate, nominal, name):
    assert waveforms.nominal_rate(rate) == nominal
    assert waveforms.group_name(waveforms.nominal_rate(rate)) == name


def test_place_snaps_segments_and_records_the_offset():
    seg = waveforms.Segment
    start = T0.value
    # 8.4 ms after the window's start: step 0, shifted +8.4 ms
    out, offset = waveforms.place(
        [seg(50.0, start + 8_400_000, np.arange(5))], start, 6, 50.0
    )
    assert out.tolist() == [0, 1, 2, 3, 4, waveforms.FILL]
    assert offset == 8_400_000
    # 12 ms is nearer step 1 (20 ms): shifted -8 ms
    out, offset = waveforms.place(
        [seg(50.0, start + 12_000_000, np.arange(5))], start, 6, 50.0
    )
    assert out.tolist() == [waveforms.FILL, 0, 1, 2, 3, 4]
    assert offset == -8_000_000


def test_place_keeps_gaps_and_lets_later_overlaps_win(caplog):
    seg = waveforms.Segment
    start = T0.value
    step = 20_000_000
    segs = [
        seg(50.0, start, np.full(4, 1)),
        seg(50.0, start + 6 * step, np.full(4, 2)),  # a gap at steps 4, 5
        seg(50.0, start + 8 * step, np.full(4, 3)),  # overlaps steps 8, 9
        seg(50.0, start - 2 * step, np.full(3, 4)),  # partly before, overwritten
    ]
    caplog.set_level(logging.INFO)
    out, _ = waveforms.place(segs, start, 12, 50.0, "AK.A..BHZ")
    f = waveforms.FILL
    assert out.tolist() == [1, 1, 1, 1, f, f, 2, 2, 3, 3, 3, 3]
    assert "AK.A..BHZ: 2 overlapping segment(s)" in caplog.text


def test_to_groups_splits_rates_and_keeps_the_last_second(fdsn, caplog):
    for c in "ZNE":
        fdsn.channel(f"AK.A..BH{c}", lat=64.0, lon=-151.0)
    fdsn.channel("AK.B..HHZ", lat=63.0, lon=-151.0, rate=100.0)
    fdsn.channel("AK.C..BHZ", lat=63.0, lon=-150.0)  # no data
    end = T0 + 10 * S - pd.Timedelta(microseconds=1)  # the runner's inclusive end
    segments = {
        # Data at a jittered rate, with metadata saying 50 Hz
        f"AK.A..BH{c}": [waveforms.Segment(49.99999, T0.value, np.arange(600))]
        for c in "ZNE"
    }
    segments["AK.B..HHZ"] = [waveforms.Segment(100.0, T0.value, np.arange(1000))]
    segments["XX.NOT..REQ"] = [waveforms.Segment(1.0, T0.value, np.arange(3))]
    caplog.set_level(logging.INFO)
    groups = waveforms.to_groups(table(fdsn), segments, T0, end, (-151.0, 63.0))

    assert list(groups) == ["50hz", "100hz"]
    g50, g100 = groups["50hz"], groups["100hz"]
    # Every step of [T0, T0 + 10 s): the last second is kept
    assert g50.sizes == {"channel": 3, "time": 500}
    assert g100.sizes == {"channel": 1, "time": 1000}
    assert g50.counts.values[0, -1] == 499
    assert g50.time.values[-1] == np.datetime64("2023-09-14T02:06:09.980")
    assert g50.channel.values.tolist() == ["AK.A..BHE", "AK.A..BHN", "AK.A..BHZ"]
    assert g50.code.values.tolist() == ["BHE", "BHN", "BHZ"]
    assert g50.location.values.tolist() == ["", "", ""]
    assert g50.sample_rate.values.tolist() == [50.0] * 3
    assert g50.azimuth.values.tolist() == [90.0, 0.0, 0.0]
    assert g50.xindexes["geometry"].crs.to_epsg() == 4326
    # One degree of latitude north of the origin; the origin is due south
    assert g50.distance_km.values[0] == pytest.approx(111.4, abs=0.2)
    assert g50.back_azimuth.values[0] == pytest.approx(180.0)
    assert g50.counts.dtype == np.int32
    assert g50.counts.encoding["_FillValue"] == waveforms.FILL
    assert g50.counts.encoding["chunks"] == (1, 500)
    assert "1 channel(s) returned no data" in caplog.text


# ---- derived views ---------------------------------------------------------------


def groups(fdsn, channels, seconds=10):
    """Groups of *channels*: ``{seed id: (rate, samples, kwargs for channel())}``,
    each starting at T0, over ``[T0, T0 + seconds)``."""
    segments = {}
    for sid, (rate, data, kw) in channels.items():
        fdsn.channel(sid, rate=rate, **kw)
        segments[sid] = [waveforms.Segment(rate, T0.value, np.asarray(data))]
    end = T0 + seconds * S - pd.Timedelta(microseconds=1)
    return waveforms.to_groups(table(fdsn), segments, T0, end, (-151.0, 63.0))


def test_to_physical_divides_by_sensitivity_and_masks_gaps(fdsn):
    data = np.full(500, 1000)
    data[10:20] = waveforms.FILL
    g = groups(fdsn, {"AK.A..BHZ": (50.0, data, {"scale": 500.0})})["50hz"]
    v = waveforms.to_physical(g)
    assert v.name == "velocity" and v.attrs["units"] == "m/s"
    assert np.isnan(v.values[0, 10:20]).all()
    assert v.values[0, 0] == 2.0
    # Also as stored: CF-encoded and read back
    assert waveforms.to_physical(waveforms.to_dataset(g)).values[0, 0] == 2.0

    g = g.assign_coords(sensitivity_units=("channel", ["m/s**2"]))
    assert waveforms.to_physical(g).name == "acceleration"
    both = groups(
        fdsn,
        {"AK.B..BHZ": (50.0, data, {}), "AK.B..BDF": (50.0, data, {"units": "Pa"})},
    )["50hz"]
    mixed = waveforms.to_physical(both)
    assert mixed.name == "physical" and "units" not in mixed.attrs


def test_to_grid_filters_only_when_downsampling(fdsn, tmp_path):
    t50, t10 = np.arange(500) / 50, np.arange(100) / 10
    data = {
        # 15 Hz is above the 8 Hz corner for a 20 Hz grid
        "AK.FAST..HHZ": (50.0, 1e6 * np.sin(2 * np.pi * 15 * t50), {}),
        # 1 Hz passes; at 10 Hz no filter is applied, so a ramp stays exact
        "AK.FAST..HHN": (50.0, 1e6 * np.sin(2 * np.pi * 1 * t50), {}),
        "AK.SLOW..LHZ": (10.0, 1000 * t10, {}),
    }
    gs = groups(fdsn, data)
    assert set(gs) == {"10hz", "50hz"}
    # Store and read back the whole tree, as a user would
    tree = xr.DataTree.from_dict(
        {f"waveforms/x/{k}": waveforms.to_dataset(g) for k, g in gs.items()}
    )
    io.save_datatree(tree, str(tmp_path / "w.zarr"))
    stored = io.open_datatree(str(tmp_path / "w.zarr"))["waveforms/x"]

    grid = waveforms.to_grid(stored, rate=20.0, physical=False)
    # 0 to 9.95 s at 20 Hz, to the last 50 Hz sample
    assert grid.sizes == {"channel": 3, "time": 200}
    assert grid.channel.values.tolist() == [
        "AK.SLOW..LHZ",
        "AK.FAST..HHN",
        "AK.FAST..HHZ",
    ]
    assert grid.sample_rate.values.tolist() == [10.0, 50.0, 50.0]
    assert grid.xindexes["geometry"].crs.to_epsg() == 4326
    assert "20 Hz" in grid.attrs["gaia:to_grid"]
    middle = slice(40, 160)  # away from the filter's edge effects
    fast = grid.counts.sel(channel="AK.FAST..HHZ").values[middle]
    passed = grid.counts.sel(channel="AK.FAST..HHN").values[middle]
    assert np.abs(fast).max() < 0.02 * 1e6
    assert np.abs(passed).max() == pytest.approx(1e6, rel=0.05)
    ramp = grid.counts.sel(channel="AK.SLOW..LHZ").values
    np.testing.assert_allclose(ramp[:199], 1000 * np.arange(199) / 20, atol=1e-6)
    assert np.isnan(ramp[199])  # after its last sample, at 9.9 s


def test_to_grid_uses_time_offset_and_leaves_gaps_nan(fdsn):
    data = 100 * np.arange(100.0)  # one count per ms
    data[40:50] = waveforms.FILL
    (g,) = groups(fdsn, {"AK.S..LHZ": (10.0, data.astype(int), {})}).values()
    g = g.assign_coords(time_offset=("channel", np.array([25], "timedelta64[ms]")))
    grid = waveforms.to_grid(g, rate=20.0, physical=False)
    v = grid.counts.values[0]
    # Sample i was recorded at 0.1 i + 0.025 s, so grid time 0.5 s holds 475
    assert v[10] == pytest.approx(475.0)
    assert np.isnan(v[0])  # before the first sample
    assert np.isnan(v[79:100]).all() and not np.isnan(v[[78, 101]]).any()


def test_by_component_maps_horizontals_by_azimuth(fdsn, caplog):
    data = np.ones(500)
    gs = groups(
        fdsn,
        {
            "AK.A..BHZ": (50.0, data, {}),
            "AK.A..BH1": (50.0, data, {"azimuth": 2.0}),
            "AK.A..BH2": (50.0, data, {"azimuth": 92.0}),
            "AK.B.00.HHZ": (50.0, 2 * data, {"lat": 63.5}),
            "AK.B.00.HH1": (50.0, 2 * data, {"lat": 63.5, "azimuth": 30.0}),
        },
    )
    grid = waveforms.to_grid(gs["50hz"], rate=50.0)
    out = waveforms.by_component(grid)
    assert out.velocity.dims == ("station", "component", "time")
    assert out.station.values.tolist() == ["AK.A", "AK.B.00"]
    assert out.component.values.tolist() == ["Z", "N", "E", "1"]
    assert out.channel.values.tolist() == [
        ["AK.A..BHZ", "AK.A..BH1", "AK.A..BH2", ""],
        ["AK.B.00.HHZ", "", "", "AK.B.00.HH1"],
    ]
    assert np.isnan(out.velocity.sel(station="AK.B.00", component="N")).all()
    assert out.distance_km.dims == ("station",)
    assert out.xindexes["geometry"].crs.to_epsg() == 4326
    assert waveforms.component("BHN", np.nan) == "N"
    assert waveforms.component("BH1", np.nan) == "1"
    with pytest.raises(ValueError, match="no waveform groups"):
        waveforms.to_grid(xr.DataTree())
