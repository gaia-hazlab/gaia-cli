import logging

import numpy as np
import pandas as pd
import pytest
import requests
from conftest import STATION_HEADER, mseed

from gaia_cli import waveforms
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
