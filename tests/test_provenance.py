import json

import numpy as np
import pandas as pd
import requests
import xarray as xr

from gaia_cli import hrrr, io, nlcd, provenance, synoptic

COLLECTION = {
    "type": "Collection",
    "id": "noaa-hrrr-analysis",
    "title": "NOAA HRRR analysis",
    "license": "CC-BY-4.0",
    "version": "0.3.0",
    "providers": [{"name": "dynamical.org", "roles": ["processor", "host"]}],
    "links": [
        {"rel": "self", "href": "https://example.com/collection.json"},
        {"rel": "license", "href": "https://creativecommons.org/licenses/by/4.0/"},
        {"rel": "about", "href": "https://dynamical.org/catalog/noaa-hrrr-analysis/"},
    ],
}


class FakeResponse:
    def __init__(self, doc):
        self.doc = doc

    def raise_for_status(self):
        pass

    def json(self):
        return self.doc


def _hourly():
    times = pd.date_range("2025-12-01T01:00", periods=3, freq="h")
    da = xr.DataArray(np.ones(3), coords={"time": times}, dims="time", name="p")
    return io.with_time_bounds(da, -1, 0, "time: sum")


def test_cli_version_is_this_checkout_commit():
    version = provenance.cli_version()
    assert len(version) == 40 and int(version, 16) >= 0


def test_run_attrs(monkeypatch):
    monkeypatch.setattr(provenance, "command", None)
    attrs = provenance.run_attrs()
    assert attrs["created"].endswith("Z")
    pd.Timestamp(attrs["created"])  # parses as RFC 3339
    assert attrs["gaia:cli-version"] == provenance.cli_version()
    assert "gaia:command" not in attrs

    provenance.set_command(["stage", "hrrr", "-i", "my aoi.json"])
    assert provenance.run_attrs()["gaia:command"] == "gaia stage hrrr -i 'my aoi.json'"


def test_time_extent_uses_bounds():
    extent = provenance.time_extent(_hourly())
    assert extent == {
        "start_datetime": "2025-12-01T00:00:00Z",
        "end_datetime": "2025-12-01T03:00:00Z",
    }


def test_collection_attrs_copies_stac_fields_and_links():
    attrs = provenance.collection_attrs(COLLECTION, "https://example.com/c.json")
    assert attrs["license"] == "CC-BY-4.0"
    assert attrs["providers"] == COLLECTION["providers"]
    assert {lk["rel"] for lk in attrs["links"]} == {"license", "about", "derived_from"}
    assert "id" not in attrs and "type" not in attrs


def test_from_stac_collection_survives_fetch_failure(monkeypatch):
    def fail(*args, **kwargs):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(provenance.requests, "get", fail)
    attrs = provenance.from_stac_collection("https://example.com/c.json")
    assert attrs == {
        "links": [
            {
                "rel": "derived_from",
                "href": "https://example.com/c.json",
                "type": "application/json",
            }
        ]
    }


def test_hrrr_source_attrs_prefer_store_version(monkeypatch):
    monkeypatch.setattr(
        provenance.requests, "get", lambda *a, **k: FakeResponse(COLLECTION)
    )
    ds = xr.Dataset(attrs={"dataset_version": "0.2.0", "license": "CC-BY-4.0"})
    ds.attrs["gaia:icechunk-snapshot"] = "SNAP"
    ds.encoding["source"] = "s3://bucket/prefix/"

    attrs = hrrr.source_attrs(ds)

    assert attrs["version"] == "0.2.0"  # the store read, not the collection's latest
    assert attrs["gaia:icechunk-snapshot"] == "SNAP"
    rels = {lk["rel"]: lk["href"] for lk in attrs["links"]}
    assert rels["derived_from"] == hrrr.STAC_COLLECTION
    assert rels["via"] == "s3://bucket/prefix/"


def test_save_datatree_root_writes_run_and_source_attrs(tmp_path):
    ds = _hourly()
    ds.attrs = nlcd.source_attrs(2021)
    path = str(tmp_path / "t.zarr")

    io.save_datatree(xr.DataTree(ds), path)
    attrs = xr.open_zarr(path).attrs

    assert attrs["gaia:cli-version"] == provenance.cli_version()
    assert attrs["license"] == "other"
    assert attrs["providers"][0]["name"] == "U.S. Geological Survey"
    assert attrs["links"][0]["rel"] == "via"


def test_save_datatree_notes_multiple_licenses(tmp_path):
    a, b = _hourly(), _hourly()
    a.attrs = {"license": "CC-BY-4.0"}
    b.attrs = {"license": "other"}
    path = str(tmp_path / "t.zarr")

    io.save_datatree(xr.DataTree.from_dict({"rasters/a": a, "rasters/b": b}), path)
    tree = io.open_datatree(path)

    assert "license" not in tree.attrs
    assert tree.attrs["gaia:license-note"] == provenance.MULTIPLE_LICENSES_NOTE
    assert "created" in tree.attrs
    assert tree["rasters/a"].attrs == {"license": "CC-BY-4.0"}


def test_derived_from_keeps_input_run():
    root = {
        "license": "other",
        **provenance.run_attrs(),
        "links": [{"rel": "via", "href": "x"}],
    }

    attrs = provenance.derived_from(root, "in.zarr")

    assert "created" not in attrs and attrs["license"] == "other"
    assert attrs["gaia:inputs"][0]["href"] == "in.zarr"
    assert attrs["gaia:inputs"][0]["created"] == root["created"]
    assert attrs["links"][-1] == {"rel": "derived_from", "href": "in.zarr"}
    # a second derivation appends rather than replacing
    again = provenance.derived_from(provenance.stamp(attrs), "aligned.zarr")
    assert [i["href"] for i in again["gaia:inputs"]] == ["in.zarr", "aligned.zarr"]


def test_synoptic_attrs_never_include_token(monkeypatch):
    monkeypatch.setattr(synoptic, "TOKEN", "secret-token-123")
    ds = _hourly()
    assert "secret-token-123" not in json.dumps(synoptic.source_attrs(ds))


def test_time_extent_of_station_reports():
    times = pd.to_datetime(["2025-12-01T00:28", "2025-12-01T03:28"])
    ds = xr.Dataset({"p": ("obs", [0.0, 1.0])}, coords={"time": ("obs", times)})
    assert provenance.time_extent(ds) == {
        "start_datetime": "2025-12-01T00:28:00Z",
        "end_datetime": "2025-12-01T03:28:00Z",
    }
