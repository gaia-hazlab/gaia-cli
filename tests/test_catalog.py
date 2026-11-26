import json
import shutil

import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from botocore.exceptions import ClientError

from gaia_cli import catalog, io, provenance


def _raster(license="CC-BY-4.0"):
    times = pd.date_range("2025-12-01T01:00", periods=2, freq="h")
    da = xr.DataArray(
        np.ones((2, 2, 2)),
        coords={"time": times, "y": [48.5, 48.0], "x": [-122.0, -121.5]},
        dims=("time", "y", "x"),
        name="precipitation",
    ).rio.write_crs("EPSG:4326")
    ds = io.with_time_bounds(da, -1, 0, "time: sum")
    ds.attrs = {
        "title": "Test raster",
        "license": license,
        "processing:lineage": "test",
        "links": [{"rel": "license", "href": f"https://example.com/{license}"}],
        **provenance.time_extent(ds),
    }
    return ds


def _read(path):
    return json.loads(path.read_text())


def test_save_datatree_root_writes_item_next_to_store(tmp_path):
    io.save_datatree(xr.DataTree(_raster()), str(tmp_path / "r.zarr"))

    item = _read(tmp_path / "r.json")
    assert item["id"] == "r"
    assert item["assets"]["data"]["href"] == "./r.zarr"
    assert item["properties"]["license"] == "CC-BY-4.0"
    assert item["properties"]["start_datetime"] == "2025-12-01T00:00:00Z"
    assert item["properties"]["datetime"] is None
    # bbox comes from the cell edges, not the cell centres
    np.testing.assert_allclose(item["bbox"], [-122.25, 47.75, -121.25, 48.75])
    assert any("processing" in e for e in item["stac_extensions"])
    assert not any("scientific" in e for e in item["stac_extensions"])


def test_tree_item_merges_groups(tmp_path):
    tree = xr.DataTree.from_dict(
        {"rasters/a": _raster("CC-BY-4.0"), "rasters/b": _raster("other")}
    )
    io.save_datatree(tree, str(tmp_path / "t.zarr"))

    item = _read(tmp_path / "t.json")
    props = item["properties"]
    assert "license" not in props
    assert props["gaia:license-note"] == provenance.MULTIPLE_LICENSES_NOTE
    assert props["description"] == "gaia store with groups: rasters/a, rasters/b"
    # a description from the stage file's attrs keeps its text
    tree.attrs["description"] = "Two rasters"
    assert provenance.stac_item(tree, "t.zarr")["properties"]["description"] == (
        "gaia store: Two rasters"
    )
    assert {lk["href"] for lk in item["links"]} == {
        "https://example.com/CC-BY-4.0",
        "https://example.com/other",
    }


class _Lifecycle:
    def __init__(self, rules=None, error=None):
        self.rules, self.error = rules, error

    def client(self, name):
        return self

    def get_bucket_lifecycle_configuration(self, Bucket):
        if self.error:
            raise self.error
        return {"Rules": self.rules}


def test_expiration_days_from_lifecycle_rules(monkeypatch):
    rules = [
        {"ID": "all", "Prefix": "", "Status": "Enabled", "Expiration": {"Days": 30}},
        {
            "ID": "tmp",
            "Filter": {"Prefix": "gaia-cli/tmp/"},
            "Status": "Enabled",
            "Expiration": {"Days": 7},
        },
        {
            "ID": "off",
            "Filter": {"Prefix": ""},
            "Status": "Disabled",
            "Expiration": {"Days": 1},
        },
        {
            "ID": "tagged",
            "Filter": {"Tag": {"Key": "k", "Value": "v"}},
            "Status": "Enabled",
            "Expiration": {"Days": 2},
        },
        {
            "ID": "uploads",
            "Prefix": "",
            "Status": "Enabled",
            "AbortIncompleteMultipartUpload": {"DaysAfterInitiation": 3},
        },
    ]
    monkeypatch.setattr(io.boto3, "Session", lambda: _Lifecycle(rules))
    assert io.expiration_days("s3://b/gaia-cli/pioneer/pioneer.zarr") == 30
    assert io.expiration_days("s3://b/gaia-cli/tmp/x/x.zarr") == 7
    assert io.expiration_days("./pioneer/pioneer.zarr") is None

    error = ClientError(
        {"Error": {"Code": "NoSuchLifecycleConfiguration"}}, "GetLifecycle"
    )
    monkeypatch.setattr(io.boto3, "Session", lambda: _Lifecycle(error=error))
    assert io.expiration_days("s3://b/gaia-cli/pioneer/pioneer.zarr") is None


def test_item_expires_with_the_bucket_rule(tmp_path, monkeypatch):
    store = str(tmp_path / "r.zarr")
    io.save_datatree(xr.DataTree(_raster()), store)
    assert "expires" not in _read(tmp_path / "r.json")["properties"]

    monkeypatch.setattr(io, "expiration_days", lambda path: 30)
    io.save_datatree(xr.DataTree(_raster()), store, overwrite=True)
    item = _read(tmp_path / "r.json")
    props = item["properties"]
    assert pd.Timestamp(props["expires"]) - pd.Timestamp(props["created"]) == (
        pd.Timedelta(days=30)
    )
    assert any("timestamps" in e for e in item["stac_extensions"])

    # Appending keeps the earlier groups' expiry
    first = "2026-01-01T00:00:00Z"
    item["properties"]["expires"] = first
    (tmp_path / "r.json").write_text(json.dumps(item))
    tree = xr.DataTree.from_dict({"rasters/b": _raster()})
    io.save_datatree(tree, store, append=True)
    assert _read(tmp_path / "r.json")["properties"]["expires"] == first


def test_rebuild_links_live_items_only(tmp_path):
    io.save_datatree(xr.DataTree(_raster()), str(tmp_path / "a.zarr"))
    (tmp_path / "event").mkdir()
    io.save_datatree(xr.DataTree(_raster()), str(tmp_path / "event" / "b.zarr"))
    io.save_datatree(xr.DataTree(_raster()), str(tmp_path / "gone.zarr"))
    shutil.rmtree(tmp_path / "gone.zarr")
    (tmp_path / "notes.json").write_text('{"not": "an item"}')

    catalog.rebuild(str(tmp_path))

    links = _read(tmp_path / "catalog.json")["links"]
    items = sorted(lk["href"] for lk in links if lk["rel"] == "item")
    assert items == ["./a.json", "./event/b.json"]


def test_find_json_skips_zarr_stores(tmp_path):
    io.save_datatree(xr.DataTree(_raster()), str(tmp_path / "a.zarr"))
    assert catalog.find_json(io.object_store(str(tmp_path))) == ["a.json"]


def test_rebuild_requires_aws_credentials_for_s3(no_aws_credentials):
    with pytest.raises(SystemExit, match="no AWS credentials"):
        catalog.rebuild("s3://bucket/scratch")
