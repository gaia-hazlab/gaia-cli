import datetime as dt
import json
import typing

import geopandas as gpd
import numpy as np
import pystac
import pytest
import rasterio
from obstore.store import LocalStore
from pystac.extensions.projection import ProjectionExtension
from rasterio.transform import from_origin
from shapely.geometry import box

from gaia_cli import prism
from gaia_cli.prism_cache import (
    CacheMiss,
    PrismCache,
    cache,
    items,
    release,
    upstream,
)

D1 = dt.date(2025, 12, 1)


def _tif_bytes(tmp_path, value=1.0) -> bytes:
    path = tmp_path / "src.tif"
    profile = {
        "driver": "GTiff",
        "width": 4,
        "height": 3,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4269",
        "transform": from_origin(-125.0, 50.0, 0.5, 0.5),
        "nodata": -9999.0,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(np.full((1, 3, 4), value, dtype="float32"))
    return path.read_bytes()


@pytest.fixture
def fake_upstream(tmp_path, monkeypatch):
    """Replace PRISM downloads and release queries with in-memory fakes."""
    state = {"release": 4, "downloads": [], "released": {}, "missing": set()}
    tif = _tif_bytes(tmp_path)

    def download(var, date, res="800m", session=None):
        if date in state["missing"]:
            raise upstream.NotReleased(date)
        state["downloads"].append((var, date))
        meta = {
            "variable": var,
            "date": date.isoformat(),
            "resolution": res,
            "release": state["release"],
            "upstream_url": upstream.zip_url(var, date, res),
            "info": {"PRISM_DATASET_VERSION": "D2"},
        }
        return upstream.Grid(tif=tif, meta=meta)

    def releases(var, start, end, res="800m", session=None):
        return {
            d: upstream.Release(d, "2026-06-25", r)
            for d, r in state["released"].items()
            if start <= d <= end
        }

    monkeypatch.setattr(upstream, "download", download)
    monkeypatch.setattr(upstream, "releases", releases)
    return state


@pytest.fixture
def store_dir(tmp_path):
    d = tmp_path / "cache"
    d.mkdir()
    return d


def _cache(store_dir, writable=True):
    return PrismCache(LocalStore(store_dir), str(store_dir), writable=writable)


def test_key_mirrors_upstream_layout():
    assert (
        PrismCache.key("ppt", D1)
        == "us/an/800m/ppt/daily/2025/prism_ppt_us_30s_20251201"
    )
    assert upstream.zip_url("tmean", D1, "4km").endswith(
        "/us/an/4km/tmean/daily/2025/prism_tmean_us_25m_20251201.zip"
    )


def test_parse_info():
    text = "PRISM_DATASET_VERSION: D2\nPRISM_DATASET_RELEASE_NUMBER: 4\nnoise\n"
    info = upstream.parse_info(text)
    assert info == {"PRISM_DATASET_VERSION": "D2", "PRISM_DATASET_RELEASE_NUMBER": "4"}


def test_ensure_fills_missing(store_dir, fake_upstream):
    c = _cache(store_dir)
    entries = c.ensure(["ppt"], "2025-12-01", "2025-12-02")
    assert sorted(entries["ppt"]) == [D1, D1 + dt.timedelta(days=1)]
    assert len(fake_upstream["downloads"]) == 2
    key = c.key("ppt", D1)
    assert (store_dir / f"{key}.tif").exists()
    item = pystac.Item.from_dict(json.loads((store_dir / f"{key}.json").read_text()))
    assert item.id == "prism_ppt_us_30s_20251201"
    assert item.properties["prism:release"] == 4
    assert item.properties["prism:dataset_version"] == "D2"
    assert "updated" in item.properties
    assert item.properties["start_datetime"] == "2025-11-30T12:00:00Z"
    assert item.assets["ppt"].href == f"{store_dir}/{key}.tif"
    assert item.get_single_link("via").href.endswith(".zip")

    # Second call is served from the cache
    c.ensure(["ppt"], "2025-12-01", "2025-12-02")
    assert len(fake_upstream["downloads"]) == 2


def test_ensure_skips_unreleased(store_dir, fake_upstream):
    fake_upstream["missing"].add(D1)
    entries = _cache(store_dir).ensure(["ppt"], D1, D1)
    assert entries["ppt"] == {}


def test_read_only_miss_raises(store_dir, fake_upstream):
    with pytest.raises(CacheMiss):
        _cache(store_dir, writable=False).ensure(["ppt"], D1, D1)
    assert fake_upstream["downloads"] == []


def _age_entry(c, store_dir, date, hours):
    path = store_dir / f"{c.key('ppt', date)}.json"
    item = json.loads(path.read_text())
    then = cache._now() - dt.timedelta(hours=hours)
    item["properties"]["updated"] = then.isoformat()
    path.write_text(json.dumps(item))


def test_revalidate_refetches_newer_release(store_dir, fake_upstream):
    c = _cache(store_dir)
    c.ensure(["ppt"], D1, D1)
    _age_entry(c, store_dir, D1, hours=48)

    fake_upstream["released"][D1] = 8
    fake_upstream["release"] = 8
    entries = c.ensure(["ppt"], D1, D1)
    assert release(entries["ppt"][D1]) == 8
    assert len(fake_upstream["downloads"]) == 2

    # Final releases are never checked again
    fake_upstream["released"][D1] = 9
    c.ensure(["ppt"], D1, D1)
    assert len(fake_upstream["downloads"]) == 2


def test_revalidate_respects_24h_limit(store_dir, fake_upstream):
    c = _cache(store_dir)
    c.ensure(["ppt"], D1, D1)
    fake_upstream["released"][D1] = 5
    c.ensure(["ppt"], D1, D1)
    assert len(fake_upstream["downloads"]) == 1


def test_read_only_stale_warns(store_dir, fake_upstream):
    c = _cache(store_dir)
    c.ensure(["ppt"], D1, D1)
    _age_entry(c, store_dir, D1, hours=48)
    fake_upstream["released"][D1] = 8
    with pytest.warns(UserWarning, match="newer releases"):
        entries = _cache(store_dir, writable=False).ensure(["ppt"], D1, D1)
    assert release(entries["ppt"][D1]) == 4


class _Resp:
    def __init__(self, body):
        self.body = body

    def raise_for_status(self):
        pass

    def json(self):
        return json.loads(self.body)


class _Session:
    def __init__(self, bodies):
        self.bodies = list(bodies)
        self.urls = []

    def get(self, url, timeout=None):
        self.urls.append(url)
        return _Resp(self.bodies.pop(0))


def test_releases_steps_back_past_newest_grid():
    end = dt.datetime.now(dt.UTC).date() - dt.timedelta(days=1)
    start = end - dt.timedelta(days=1)
    row = [start.isoformat(), "2026-10-13", "ppt", "1", "url"]
    session = _Session(["Invalid date", json.dumps(row)])
    out = upstream.releases("ppt", start, end, session=session)
    assert out == {start: upstream.Release(start, "2026-10-13", 1)}
    assert len(session.urls) == 2


def test_items(store_dir, fake_upstream):
    c = _cache(store_dir)
    entries = c.ensure(["ppt", "tmean"], D1, D1 + dt.timedelta(days=1))
    out = items(entries)
    assert [i.id for i in out] == ["prism_us_30s_20251201", "prism_us_30s_20251202"]
    item = out[0]
    assert set(item.assets) == {"ppt", "tmean"}
    assert ProjectionExtension.ext(item).epsg == 4269
    assert item.properties["proj:shape"] == [3, 4]
    assert item.assets["ppt"].extra_fields["raster:bands"][0]["nodata"] == -9999.0
    assert item.datetime == dt.datetime(2025, 12, 1, tzinfo=dt.UTC)


def _read(path):
    return json.loads(path.read_text())


def test_fill_updates_catalogs(store_dir, fake_upstream):
    c = _cache(store_dir)
    c.ensure(["ppt", "tmean"], D1, D1 + dt.timedelta(days=1))

    coll = _read(store_dir / "collection.json")
    children = [lk["href"] for lk in coll["links"] if lk["rel"] == "child"]
    assert children == [
        f"{store_dir}/us/an/800m/ppt/daily/2025/catalog.json",
        f"{store_dir}/us/an/800m/tmean/daily/2025/catalog.json",
    ]
    assert coll["extent"]["temporal"]["interval"] == [
        ["2025-11-30T12:00:00Z", "2025-12-02T12:00:00Z"]
    ]
    assert coll["summaries"]["prism:variable"] == ["ppt", "tmean"]

    # The whole tree is walkable from the collection
    root = pystac.Collection.from_file(str(store_dir / "collection.json"))
    found = sorted(i.id for i in root.get_items(recursive=True))
    assert found == [
        "prism_ppt_us_30s_20251201",
        "prism_ppt_us_30s_20251202",
        "prism_tmean_us_30s_20251201",
        "prism_tmean_us_30s_20251202",
    ]


def test_rebuild_catalog_matches_incremental(store_dir, fake_upstream):
    c = _cache(store_dir)
    c.ensure(["ppt", "tmean"], D1, D1 + dt.timedelta(days=1))
    files = ["collection.json", "us/an/800m/ppt/daily/2025/catalog.json"]
    before = [_read(store_dir / f) for f in files]
    for f in files:
        (store_dir / f).unlink()

    assert c.rebuild_catalog() == 4
    assert [_read(store_dir / f) for f in files] == before


@pytest.fixture
def aws_config(tmp_path, monkeypatch):
    """Point boto3 at a throwaway config with a Kopah and a B2 profile."""
    config = tmp_path / "config"
    config.write_text(
        "[profile kopah]\n"
        "endpoint_url = https://s3.kopah.uw.edu\n"
        "[profile b2]\n"
        "region = us-west-004\n"
        "endpoint_url = https://s3.us-west-004.backblazeb2.com\n"
    )
    creds = tmp_path / "credentials"
    creds.write_text(
        "[kopah]\naws_access_key_id = k\naws_secret_access_key = s\n"
        "[b2]\naws_access_key_id = k\naws_secret_access_key = s\n"
    )
    monkeypatch.setenv("AWS_CONFIG_FILE", str(config))
    monkeypatch.setenv("AWS_SHARED_CREDENTIALS_FILE", str(creds))
    for var in (
        "AWS_PROFILE",
        "AWS_ACCESS_KEY_ID",
        "AWS_ENDPOINT_URL",
        "PRISM_CACHE_URL",
        "PRISM_CACHE_ENDPOINT",
        "PRISM_CACHE_PROFILE",
    ):
        monkeypatch.delenv(var, raising=False)
    return monkeypatch


def test_from_env_defaults_to_anonymous_kopah(aws_config):
    c = PrismCache.from_env()
    assert c.href_base == "https://s3.kopah.uw.edu/gaia/public/prism"
    assert not c.writable


def test_from_env_output_profile_for_other_store_leaves_cache_read_only(aws_config):
    aws_config.setenv("AWS_PROFILE", "b2")
    c = PrismCache.from_env()
    assert c.href_base == "https://s3.kopah.uw.edu/gaia/public/prism"
    assert not c.writable


def test_from_env_output_profile_for_kopah_writes_cache(aws_config):
    aws_config.setenv("AWS_PROFILE", "kopah")
    assert PrismCache.from_env().writable


def test_from_env_cache_profile_supplies_endpoint(aws_config):
    aws_config.setenv("AWS_PROFILE", "kopah")
    aws_config.setenv("PRISM_CACHE_PROFILE", "b2")
    aws_config.setenv("PRISM_CACHE_URL", "s3://gaia-scratch/prism-cache")
    c = PrismCache.from_env()
    assert c.href_base == (
        "https://s3.us-west-004.backblazeb2.com/gaia-scratch/prism-cache"
    )
    assert c.writable


def test_prism_load_temperatures(store_dir, fake_upstream, monkeypatch):
    monkeypatch.setattr(PrismCache, "from_env", lambda: _cache(store_dir))
    aoi = gpd.GeoDataFrame(geometry=[box(-125.0, 48.5, -123.0, 50.0)], crs="EPSG:4326")
    ds = prism.load(aoi, "2025-12-01", "2025-12-02", ("ppt", "tmin", "tmax"))

    assert set(ds.data_vars) == {"precipitation", "tmin", "tmax", "time_bounds"}
    assert ds.sizes["time"] == 2
    methods = {
        v: ds[v].attrs["cell_methods"] for v in ("precipitation", "tmin", "tmax")
    }
    assert methods == {
        "precipitation": "time: sum",
        "tmin": "time: minimum",
        "tmax": "time: maximum",
    }
    assert ds.tmin.attrs["units"] == "degC"
    assert ds.tmax.attrs["standard_name"] == "air_temperature"
    assert ds.tmax.attrs["long_name"] == "PRISM daily maximum temperature"
    assert ds.tmin.attrs["prism_release"] == [4, 4]


def test_prism_variables_are_published_and_described():
    variables = set(typing.get_args(prism.Variable))
    assert variables == {"ppt", "tmin", "tmax", "tmean"}
    assert variables <= set(upstream.VARIABLES)
    assert set(prism.STANDARD_NAMES) == set(prism.LONG_NAMES) == variables
