import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

from gaia_cli.datasets import synoptic

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


def test_get_station_timeseries_requests_derived_precip_only_with_precip(monkeypatch):
    calls = []

    def fake_get(url, params):
        calls.append(dict(params))
        return FakeResponse({"SUMMARY": OK, "STATION": []})

    monkeypatch.setattr(synoptic.requests, "get", fake_get)
    synoptic.get_station_timeseries(["A"], "2025-12-01", "2025-12-02", ["air_temp"])
    synoptic.get_station_timeseries(["A"], "2025-12-01", "2025-12-02")

    assert calls[0]["vars"] == "air_temp" and "precip" not in calls[0]
    assert calls[1]["vars"] == "precip_accum" and calls[1]["precip"] == 1


def test_construct_ragged_from_raw_precip():
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

    ds = synoptic.construct_ragged({"STATION": [hourly, quarter_hourly]})

    # Every report is kept
    assert ds.station.values.tolist() == ["HOURLY", "QUARTER"]
    assert ds.row_size.values.tolist() == [4, 13]
    assert ds.precipitation.values[:4].tolist() == [0.0, 1.0, 1.0, 3.0]
    assert ds.precipitation.values[-1] == 6.0
    assert ds.precipitation.attrs["units"] == "mm"
    assert ds.elevation.values.tolist() == pytest.approx([304.8, 304.8])
    assert ds.precip_source.values.tolist() == ["precip_accum", "precip_accum"]


def test_construct_ragged_prefers_derived_precip():
    station = _station(
        "DERIVED",
        {
            "date_time": ["2025-12-01T00:00:00Z", "2025-12-01T01:00:00Z"],
            "precip_accum_set_1": [100.0, 200.0],
            "precip_accumulated_set_1d": [0.0, 1.5],
        },
    )
    station["SENSOR_VARIABLES"] = {
        "precip_accumulated": {
            "precip_accumulated_set_1d": {"derived_from": "precip_accum_one_hour"}
        }
    }
    ds = synoptic.construct_ragged({"STATION": [station]})
    assert ds.precipitation.values.tolist() == [0.0, 1.5]
    assert ds.precip_source.values.tolist() == ["precip_accum_one_hour"]


def test_construct_ragged_keeps_other_variables():
    both = _station(
        "BOTH",
        {
            "date_time": ["2025-12-01T00:00:00Z", "2025-12-01T00:15:00Z"],
            "air_temp_set_1": [2.5, 3.0],
            "air_temp_set_2": [9.0, 9.0],
            "precip_accumulated_set_1d": [0.0, None],
        },
    )
    temp_only = _station(
        "TEMP",
        {"date_time": ["2025-12-01T00:00:00Z"], "air_temp_set_1": [1.0]},
    )
    neither = _station(
        "NEITHER",
        {"date_time": ["2025-12-01T00:00:00Z"], "wind_speed_set_1": [1.0]},
    )
    data = {"UNITS": {"air_temp": "Celsius"}, "STATION": [both, temp_only, neither]}

    ds = synoptic.construct_ragged(data, ["precip_accum", "air_temp"])

    # Stations without any of the variables are dropped
    assert ds.station.values.tolist() == ["BOTH", "TEMP"]
    assert ds.row_size.values.tolist() == [2, 1]
    # The first sensor set
    assert ds.air_temp.values.tolist() == [2.5, 3.0, 1.0]
    assert ds.air_temp.attrs["units"] == "Celsius"
    assert np.isnan(ds.precipitation.values[1:]).all()


def test_construct_ragged_raises_without_any_reports():
    station = _station("NONE", {"date_time": ["2025-12-01T00:00:00Z"]})
    with pytest.raises(ValueError, match="no station reported"):
        synoptic.construct_ragged({"STATION": [station]}, ["air_temp"])


def test_load_sets_provenance_and_encodes_for_the_store(monkeypatch):
    monkeypatch.setattr(synoptic, "TOKEN", "token")
    monkeypatch.setattr(
        synoptic,
        "get_station_metadata",
        lambda aoi, variables: gpd.GeoDataFrame({"stid": ["A"]}),
    )
    station = _station(
        "A",
        {
            "date_time": ["2025-12-01T00:00:00Z", "2025-12-01T01:00:00Z"],
            "precip_accumulated_set_1d": [0.0, 1.5],
        },
    )
    monkeypatch.setattr(
        synoptic,
        "get_station_timeseries",
        lambda ids, start, end, variables: {"STATION": [station]},
    )
    aoi = gpd.GeoDataFrame(geometry=[box(-122, 48, -121, 49)], crs="EPSG:4326")

    ds = synoptic.load(aoi, "2025-12-01", "2025-12-01T23:59:59")
    assert ds.attrs["start_datetime"] == "2025-12-01T00:00:00Z"
    assert ds.attrs["end_datetime"] == "2025-12-01T01:00:00Z"

    encoded = synoptic.to_dataset(ds)
    assert "geometry_container" in encoded
    assert encoded.attrs == ds.attrs


def test_get_station_metadata_raises_on_api_error(monkeypatch):
    error = {"RESPONSE_CODE": 2, "RESPONSE_MESSAGE": "Invalid token"}
    monkeypatch.setattr(
        synoptic.requests, "get", lambda url, params: FakeResponse({"SUMMARY": error})
    )
    aoi = gpd.GeoDataFrame(geometry=[box(-122, 48, -121, 49)], crs="EPSG:4326")
    with pytest.raises(RuntimeError, match="Invalid token"):
        synoptic.get_station_metadata(aoi)
