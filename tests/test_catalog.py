import json
import shutil

import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr

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


def test_save_zarr_writes_item_next_to_store(tmp_path):
    io.save_zarr(_raster(), str(tmp_path / "r.zarr"))

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
    assert props["description"].endswith("rasters/a, rasters/b")
    assert {lk["href"] for lk in item["links"]} == {
        "https://example.com/CC-BY-4.0",
        "https://example.com/other",
    }


def test_rebuild_links_live_items_only(tmp_path):
    io.save_zarr(_raster(), str(tmp_path / "a.zarr"))
    (tmp_path / "event").mkdir()
    io.save_zarr(_raster(), str(tmp_path / "event" / "b.zarr"))
    io.save_zarr(_raster(), str(tmp_path / "gone.zarr"))
    shutil.rmtree(tmp_path / "gone.zarr")
    (tmp_path / "notes.json").write_text('{"not": "an item"}')

    catalog.rebuild(str(tmp_path))

    links = _read(tmp_path / "catalog.json")["links"]
    items = sorted(lk["href"] for lk in links if lk["rel"] == "item")
    assert items == ["./a.json", "./event/b.json"]


def test_find_json_skips_zarr_stores(tmp_path):
    io.save_zarr(_raster(), str(tmp_path / "a.zarr"))
    assert catalog.find_json(io.object_store(str(tmp_path))) == ["a.json"]


def test_rebuild_requires_aws_profile_for_s3(monkeypatch):
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    with pytest.raises(ValueError):
        catalog.rebuild("s3://bucket/scratch")
