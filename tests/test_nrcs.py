import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

from gaia_cli import fetch, stations
from gaia_cli.datasets import nrcs

AOI = gpd.GeoDataFrame(geometry=[box(-122, 46, -117, 49)], crs="EPSG:4326")


def _site(triplet, lon, lat, begin="2000-01-01 00:00", end="2100-01-01 00:00"):
    return {
        "stationTriplet": triplet,
        "name": triplet.split(":")[0],
        "networkCode": triplet.split(":")[2],
        "elevation": 1000.0,
        "longitude": lon,
        "latitude": lat,
        "dataTimeZone": -8.0,
        "beginDate": begin,
        "endDate": end,
    }


def _series(code, unit, values, depth=None, ordinal=1):
    element = {"elementCode": code, "ordinal": ordinal, "storedUnitCode": unit}
    if depth is not None:
        element["heightDepth"] = depth
    return {
        "stationElement": element,
        "values": [{"date": d, "value": v} for d, v in values],
    }


class FakeResponse:
    def __init__(self, data):
        self.data = data

    def raise_for_status(self):
        pass

    def json(self):
        return self.data


SITES = [
    _site("1:WA:SCAN", -118.0, 47.0),
    _site("2:WA:SNTL", -121.0, 47.5),
    _site("3:OR:SNTL", -121.0, 44.0),  # outside the AOI
    _site("4:WA:SNTL", -120.0, 48.0, end="2010-01-01 00:00"),  # closed before
]


def test_get_stations_clips_to_aoi_and_window(monkeypatch):
    calls = []

    def fake_get(url, params):
        calls.append(params)
        return FakeResponse(SITES)

    monkeypatch.setattr(fetch, "get", fake_get)
    sites = nrcs.get_stations(AOI, "2026-10-01", "2026-10-31", ["SMS", "STO"])

    assert sites.stationTriplet.tolist() == ["1:WA:SCAN", "2:WA:SNTL"]
    assert calls[0]["elements"] == "SMS:*,STO:*"
    assert calls[0]["activeOnly"] == "false"  # closed stations keep their records
    assert "*:*:SCAN" in calls[0]["stationTriplets"]


def test_get_data_starts_early_and_chunks(monkeypatch):
    calls = []

    def fake_get(url, params):
        calls.append(params)
        return FakeResponse(
            [{"stationTriplet": t} for t in params["stationTriplets"].split(",")]
        )

    monkeypatch.setattr(fetch, "get", fake_get)
    triplets = [f"{i}:WA:SNTL" for i in range(nrcs.CHUNK + 1)]
    data = nrcs.get_data(triplets, "2026-10-01", "2026-10-31T23:59:59.999999")

    assert len(calls) == 2 and len(data) == nrcs.CHUNK + 1
    # Local standard time is behind UTC, so the request starts earlier
    assert calls[0]["beginDate"] == "2026-09-30 12:00"
    assert calls[0]["endDate"] == "2026-10-31 23:59"
    assert calls[0]["duration"] == "HOURLY"


def test_name_by_depth_height_and_ordinal():
    assert nrcs._name({"elementCode": "SMS", "heightDepth": -2}) == "SMS_5cm"
    assert nrcs._name({"elementCode": "SMS", "heightDepth": -20}) == "SMS_51cm"
    assert nrcs._name({"elementCode": "STO", "heightDepth": -8, "ordinal": 2}) == (
        "STO_20cm_2"
    )
    assert nrcs._name({"elementCode": "WSPDV", "heightDepth": 120}) == (
        "WSPDV_305cm_above"
    )
    assert nrcs._name({"elementCode": "WTEQ", "ordinal": 1}) == "WTEQ"


def test_construct_ragged_converts_units_and_times():
    sites = gpd.GeoDataFrame(
        SITES[:2],
        geometry=gpd.points_from_xy([-118.0, -121.0], [47.0, 47.5]),
        crs="EPSG:4326",
    )
    data = [
        {
            "stationTriplet": "1:WA:SCAN",
            "data": [
                _series(
                    "SMS",
                    "pct",
                    # 15:00 PST is 23:00 UTC, the day before the window
                    [("2026-09-30 15:00", 10.0), ("2026-09-30 16:00", 20.0)],
                    depth=-2,
                ),
                _series(
                    "STO",
                    "degF",
                    [("2026-09-30 16:00", 50.0), ("2026-09-30 17:00", None)],
                    depth=-2,
                ),
            ],
        },
        {
            "stationTriplet": "2:WA:SNTL",
            "data": [_series("WTEQ", "in", [("2026-10-01 00:00", 1.0)])],
        },
    ]

    ds = nrcs.construct_ragged(sites, data, "2026-10-01", "2026-10-31")

    assert ds.station.values.tolist() == ["1:WA:SCAN", "2:WA:SNTL"]
    assert ds.row_size.values.tolist() == [1, 1]
    assert ds.time.values.astype("datetime64[h]").astype(str).tolist() == [
        "2026-10-01T00",
        "2026-10-01T08",
    ]
    assert ds.SMS_5cm.values[0] == pytest.approx(0.20)
    assert ds.STO_5cm.values[0] == pytest.approx(10.0)
    assert ds.WTEQ.values[1] == pytest.approx(25.4)
    assert np.isnan(ds.SMS_5cm.values[1])
    assert ds.SMS_5cm.attrs["units"] == "m3 m-3"
    assert ds.SMS_5cm.attrs["depth"] == "5.08 cm"
    assert ds.STO_5cm.attrs["units"] == "degC"
    assert ds.WTEQ.attrs["units"] == "mm"
    assert ds.elevation.values[0] == pytest.approx(304.8)
    # The stations helpers read it
    assert stations.to_grid(ds, freq="1h").SMS_5cm.sizes == {"geometry": 2, "time": 9}


def test_construct_ragged_raises_without_reports():
    sites = gpd.GeoDataFrame(
        SITES[:1], geometry=gpd.points_from_xy([-118.0], [47.0]), crs="EPSG:4326"
    )
    data = [{"stationTriplet": "1:WA:SCAN", "data": []}]
    with pytest.raises(ValueError, match="no NRCS station"):
        nrcs.construct_ragged(sites, data, "2026-10-01", "2026-10-31")


def test_load(monkeypatch):
    data = [
        {
            "stationTriplet": "1:WA:SCAN",
            "data": [_series("SMS", "pct", [("2026-10-01 00:00", 30.0)], depth=-8)],
        }
    ]

    def fake_get(url, params):
        return FakeResponse(SITES if url.endswith("/stations") else data)

    monkeypatch.setattr(fetch, "get", fake_get)
    ds = nrcs.load(AOI, "2026-10-01", "2026-10-31")

    assert ds.SMS_20cm.values.tolist() == [pytest.approx(0.30)]
    assert ds.attrs["license"] == "other"
    assert ds.attrs["start_datetime"].startswith("2026-10-01T08")


def test_load_raises_without_stations(monkeypatch):
    monkeypatch.setattr(fetch, "get", lambda url, params: FakeResponse(SITES[2:3]))
    with pytest.raises(ValueError, match="no NRCS station with SMS"):
        nrcs.load(AOI, "2026-10-01", "2026-10-31")
