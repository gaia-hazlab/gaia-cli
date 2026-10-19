import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

from gaia_cli import synoptic

OK = {"RESPONSE_CODE": 1, "RESPONSE_MESSAGE": "OK"}
FORBIDDEN = {"RESPONSE_CODE": 403, "RESPONSE_MESSAGE": "no access"}


def _station(stid, observations):
    return {
        "STID": stid,
        "NAME": stid,
        "ELEVATION": "1000",
        "LONGITUDE": "-121.5",
        "LATITUDE": "48.5",
        "OBSERVATIONS": observations,
    }


class FakeResponse:
    def __init__(self, data):
        self.data = data

    def raise_for_status(self):
        pass

    def json(self):
        return self.data


def test_check_response_raises_on_api_error():
    with pytest.raises(RuntimeError, match="403: no access"):
        synoptic._check_response({"SUMMARY": FORBIDDEN})


def test_accumulate_raw_precip_handles_jitter_and_reset():
    raw = pd.Series([485.14, 485.14, 482.60, 485.14, 487.68, 30.0, 32.54])
    expected = [0.0, 0.0, 0.0, 0.0, 2.54, 2.54, 5.08]
    result = synoptic.accumulate_raw_precip(raw)
    assert result.round(2).tolist() == expected


def test_accumulate_raw_precip_keeps_missing_values():
    raw = pd.Series([10.0, None, 12.0])
    result = synoptic.accumulate_raw_precip(raw)
    assert result.isna().tolist() == [False, True, False]
    assert result.iloc[-1] == 2.0


def test_get_station_timeseries_falls_back_without_derived_precip(monkeypatch):
    calls = []

    def fake_get(url, params):
        calls.append(dict(params))
        if "precip" in params:
            return FakeResponse({"SUMMARY": FORBIDDEN})
        return FakeResponse({"SUMMARY": OK, "STATION": []})

    monkeypatch.setattr(synoptic.requests, "get", fake_get)
    data = synoptic.get_station_timeseries(["A"], "2025-12-01", "2025-12-02")

    assert data["STATION"] == []
    assert len(calls) == 2
    assert "precip" not in calls[1]


def test_get_station_timeseries_raises_on_other_errors(monkeypatch):
    error = {"RESPONSE_CODE": 2, "RESPONSE_MESSAGE": "bad request"}
    monkeypatch.setattr(
        synoptic.requests, "get", lambda url, params: FakeResponse({"SUMMARY": error})
    )
    with pytest.raises(RuntimeError, match="bad request"):
        synoptic.get_station_timeseries(["A"], "2025-12-01", "2025-12-02")


def test_construct_uniform_datarray_from_raw_precip():
    hourly = _station(
        "HOURLY",
        {
            "date_time": [f"2025-12-01T0{h}:00:00Z" for h in range(4)],
            "precip_accum_set_1": [100.0, 101.0, 101.0, 103.0],
        },
    )
    quarter_hourly = _station(
        "QUARTER",
        {
            "date_time": pd.date_range("2025-12-01", periods=13, freq="15min")
            .strftime("%Y-%m-%dT%H:%M:%SZ")
            .tolist(),
            "precip_accum_set_1": [50.0 + 0.5 * i for i in range(13)],
        },
    )

    da = synoptic.construct_uniform_datarray({"STATION": [hourly, quarter_hourly]})

    assert da.sizes == {"station": 2, "time": 4}
    assert da.sel(station="HOURLY").values.tolist() == [0.0, 1.0, 1.0, 3.0]
    assert da.sel(station="QUARTER").values.tolist() == [0.0, 2.0, 4.0, 6.0]


def test_construct_uniform_datarray_off_hour_reports_never_include_later_precip():
    station = _station(
        "OFFHOUR",
        {
            "date_time": [f"2025-12-01T0{h}:28:00Z" for h in range(4)],
            "precip_accumulated_set_1d": [0.0, 1.0, 3.0, 6.0],
        },
    )
    da = synoptic.construct_uniform_datarray({"STATION": [station]}).sel(
        station="OFFHOUR"
    )

    # 01:00 holds the 00:28 report, not the 01:28 one
    assert da.time.dt.strftime("%H:%M").values.tolist() == [
        "01:00",
        "02:00",
        "03:00",
        "04:00",
    ]
    assert da.values.tolist() == [0.0, 1.0, 3.0, 6.0]
    assert da.report_time.dt.strftime("%H:%M").values.tolist() == [
        "00:28",
        "01:28",
        "02:28",
        "03:28",
    ]


def test_construct_uniform_datarray_leaves_hours_without_a_report_missing():
    station = _station(
        "GAP",
        {
            "date_time": [
                "2025-12-01T00:00:00Z",
                "2025-12-01T01:00:00Z",
                "2025-12-01T04:00:00Z",
            ],
            "precip_accumulated_set_1d": [0.0, 1.0, 9.0],
        },
    )
    da = synoptic.construct_uniform_datarray({"STATION": [station]}).sel(station="GAP")

    # The 8 mm reported at 04:00 is not assigned to 02:00 or 03:00
    assert da.sizes["time"] == 5
    assert da.isnull().values.tolist() == [False, False, True, True, False]
    assert da.sel(time="2025-12-01T04:00").item() == 9.0


def test_construct_uniform_datarray_prefers_derived_precip():
    station = _station(
        "DERIVED",
        {
            "date_time": ["2025-12-01T00:00:00Z", "2025-12-01T01:00:00Z"],
            "precip_accum_set_1": [100.0, 200.0],
            "precip_accumulated_set_1d": [0.0, 1.5],
        },
    )
    da = synoptic.construct_uniform_datarray({"STATION": [station]})
    assert da.sel(station="DERIVED").values.tolist() == [0.0, 1.5]


def test_get_station_metadata_raises_on_api_error(monkeypatch):
    error = {"RESPONSE_CODE": 2, "RESPONSE_MESSAGE": "Invalid token"}
    monkeypatch.setattr(
        synoptic.requests, "get", lambda url, params: FakeResponse({"SUMMARY": error})
    )
    aoi = gpd.GeoDataFrame(geometry=[box(-122, 48, -121, 49)], crs="EPSG:4326")
    with pytest.raises(RuntimeError, match="Invalid token"):
        synoptic.get_station_metadata(aoi)
