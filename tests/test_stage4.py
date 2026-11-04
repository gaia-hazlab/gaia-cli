import datetime as dt
import gzip
import io
import json
import tarfile

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from shapely.geometry import box

from gaia_cli.datasets import stage4

BOLT_CREEK = gpd.GeoDataFrame(geometry=[box(-121.6, 47.6, -121.0, 47.9)], crs=4326)
ANCHORAGE = gpd.GeoDataFrame(geometry=[box(-150.0, 61.0, -149.0, 61.5)], crs=4326)


@pytest.fixture(autouse=True)
def cache(monkeypatch, tmp_path):
    monkeypatch.setenv("STAGE4_CACHE", str(tmp_path / "cache"))
    return tmp_path / "cache"


# ---- domain -------------------------------------------------------------------


def test_domain_follows_the_aoi():
    assert stage4.domain(BOLT_CREEK) == "conus"
    assert stage4.domain(ANCHORAGE) == "ak"


@pytest.mark.parametrize(
    "bbox",
    [
        (-158.0, 19.0, -155.0, 22.0),  # Hawaii
        (-150.0, 40.0, -100.0, 65.0),  # spans both grids
    ],
)
def test_domain_outside_both_grids(bbox):
    with pytest.raises(ValueError, match="not within"):
        stage4.domain(gpd.GeoDataFrame(geometry=[box(*bbox)], crs=4326))


def test_extent_bboxes_contain_the_grids():
    union, conus, ak = stage4.EXTENT["spatial"]["bbox"]
    for aoi, bbox in [(BOLT_CREEK, conus), (ANCHORAGE, ak)]:
        w, s, e, n = aoi.total_bounds
        assert bbox[0] <= w and e <= bbox[2] and bbox[1] <= s and n <= bbox[3]
    assert union == [min(conus[0], ak[0]), conus[1], conus[2], ak[3]]


def test_cache_defaults_to_the_working_directory(monkeypatch, tmp_path):
    monkeypatch.delenv("STAGE4_CACHE")
    monkeypatch.chdir(tmp_path)
    assert stage4.cache_dir().resolve() == tmp_path / ".cache" / "gaia" / "stage4"


# ---- archive ------------------------------------------------------------------


def test_member_names_cover_each_era():
    noon = dt.datetime(2016, 1, 1, 12)
    assert stage4.member_names(noon, "conus") == [
        "st4_conus.2016010112.24h.grb2",
        "ST4.2016010112.24h.gz",
        "ST4.2016010112.24h.Z",
    ]
    assert stage4.member_names(noon, "ak", 6) == [
        "st4_ak.2016010112.06h.grb2",
        "st4_ak.2016010112.06h.gz",
    ]
    assert stage4.member_names(dt.datetime(2016, 1, 1, 7), "conus", 1)[1] == (
        "ST4.2016010107.01h.gz"
    )


def test_period_ends_are_those_in_the_daily_tar():
    day = dt.date(2016, 1, 1)
    hours = lambda h: [t.hour for t in stage4.period_ends(day, h)]
    assert hours(24) == [12]
    assert hours(6) == [0, 6, 12, 18]
    assert hours(1) == list(range(24))


def _tar(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


class FakeResponse:
    def __init__(self, status_code, content=b"", headers=None):
        self.status_code, self.content, self.headers = (
            status_code,
            content,
            headers or {},
        )

    def raise_for_status(self):
        if self.status_code >= 400:
            raise OSError(self.status_code)


class FakeArchive:
    """A requests.Session serving monthly tars by HTTP range."""

    def __init__(self, files: dict[str, bytes]):
        self.files, self.gets = files, 0

    def _file(self, url):
        return self.files.get(url.rsplit("/", 1)[-1])

    def head(self, url, timeout=None):
        data = self._file(url)
        if data is None:
            return FakeResponse(404)
        return FakeResponse(200, headers={"Content-Length": str(len(data))})

    def get(self, url, headers, timeout=None):
        self.gets += 1
        start, end = map(int, headers["Range"].removeprefix("bytes=").split("-"))
        return FakeResponse(206, self._file(url)[start : end + 1])


def archive():
    day = _tar(
        {
            "ST4.2016010100.06h.gz": gzip.compress(b"6h 00"),
            "ST4.2016010111.01h.gz": gzip.compress(b"hourly 11"),
            "ST4.2016010112.24h.gz": gzip.compress(b"daily grib"),
            "ST4.2016010112.06h.gz": gzip.compress(b"6h 12"),
        }
    )
    other = _tar({"ST4.2016010212.24h.gz": gzip.compress(b"next day")})
    # Members are out of date order, with paths, like the archive's
    return FakeArchive(
        {"stage4.201601.tar": _tar({"./ST4.20160102": other, "./ST4.20160101": day})}
    )


def test_fetch_day_extracts_one_grid_by_range(cache):
    session = archive()
    noon = pd.Timestamp("2016-01-01T12")
    paths = stage4.fetch_day(dt.date(2016, 1, 1), "conus", session=session)
    assert paths == {noon: cache / "2016" / "ST4.2016010112.24h"}
    assert paths[noon].read_bytes() == b"daily grib"
    index = json.loads((cache / "index" / "stage4.201601.json").read_text())
    assert set(index) == {"ST4.20160101", "ST4.20160102"}

    # Cached: no more requests
    gets = session.gets
    assert stage4.fetch_day(dt.date(2016, 1, 1), "conus", session=session) == paths
    assert session.gets == gets


def test_fetch_day_sub_daily_keeps_the_grids_it_has(cache):
    session = archive()
    six = stage4.fetch_day(dt.date(2016, 1, 1), "conus", 6, session)
    assert sorted(t.hour for t in six) == [0, 12]  # 06 and 18 aren't in the tar
    assert six[pd.Timestamp("2016-01-01T12")].read_bytes() == b"6h 12"
    hourly = stage4.fetch_day(dt.date(2016, 1, 1), "conus", 1, session)
    assert list(hourly) == [pd.Timestamp("2016-01-01T11")]


def test_fetch_day_missing(cache):
    session = archive()
    with pytest.raises(FileNotFoundError, match="ST4.20160103 is not in"):
        stage4.fetch_day(dt.date(2016, 1, 3), "conus", session=session)
    assert stage4.fetch_day(dt.date(2016, 1, 1), "ak", session=session) == {}
    with pytest.raises(FileNotFoundError, match="stage4.201602.tar"):
        stage4.fetch_day(dt.date(2016, 2, 1), "conus", session=session)


def test_month_index_skips_loose_files(cache):
    # Like 2019-01: every daily tar, each followed by loose copies of its grids
    members = {}
    for d in range(1, 32):
        members[f"./ST4.201901{d:02d}"] = _tar({"x": b"day"})
        for h in range(5):
            members[f"./ST4.201901{d:02d}{h:02d}.01h.gz"] = b"g" * 3000
    session = FakeArchive({"stage4.201901.tar": _tar(members)})
    index = stage4.month_index("201901", session)
    assert sorted(index) == [f"ST4.201901{d:02d}" for d in range(1, 32)]
    # Headers are read ahead in 1 MB blocks, not one request each
    assert session.gets < 10


def test_decompress():
    assert stage4._decompress("a.24h.gz", gzip.compress(b"x")) == b"x"
    assert stage4._decompress("a.24h.grb2", b"x") == b"x"
    # 'abc' compressed with Unix compress (LZW)
    lzw = bytes.fromhex("1f9d9061c48c01")
    assert stage4._decompress("a.24h.Z", lzw) == b"abc"


# ---- reading ------------------------------------------------------------------


def _geotiff(path, shape, value):
    with rasterio.open(
        path, "w", driver="GTiff", width=shape[1], height=shape[0], count=1,
        dtype="float64", nodata=9999, crs=stage4.CRS,
        transform=rasterio.transform.from_origin(0, 0, 4762, 4762),
    ) as dst:  # fmt: skip
        data = np.full(shape, value, "float64")
        data[0, 0] = 9999
        dst.write(data, 1)
    return path


def test_read_grib_puts_values_on_the_canonical_grid(tmp_path):
    path = _geotiff(tmp_path / "g.tif", (530, 460), 2.5)
    da = stage4.read_grib(path, dt.date(2024, 8, 1), "ak")
    x, y = stage4.grid("ak")
    np.testing.assert_array_equal(da.x, x)
    np.testing.assert_array_equal(da.y, y)
    assert da.dtype == np.float32 and da.attrs["units"] == "mm"
    assert np.isnan(da[0, 0, 0]) and float(da[0, 1, 1]) == 2.5
    assert str(da.time.values[0])[:10] == "2024-08-01"


def test_read_grib_rejects_another_grid(tmp_path):
    path = _geotiff(tmp_path / "g.tif", (10, 10), 1.0)
    with pytest.raises(ValueError, match="isn't the Stage IV conus grid"):
        stage4.read_grib(path, dt.date(2024, 8, 1), "conus")


# ---- load ---------------------------------------------------------------------


@pytest.fixture
def fake_days(monkeypatch):
    """Serve each grid as its period's end hour + day/100, except missing ones.

    Add dates to ``missing["days"]`` (no daily tar), period ends to
    ``missing["ends"]`` (grid not in the tar), or period ends to
    ``missing["empty"]`` (grid all NaN, as for hourly over the West).
    """
    missing = {"days": set(), "ends": set(), "empty": set()}

    def fetch_day(day, domain, hours=24, session=None):
        if day in missing["days"]:
            raise FileNotFoundError(f"ST4.{day:%Y%m%d} is not in the archive")
        return {
            e: e for e in stage4.period_ends(day, hours) if e not in missing["ends"]
        }

    def read_grib(path, label, domain):
        x, y = stage4.grid(domain)
        value = np.nan if path in missing["empty"] else path.hour + path.day / 100
        return stage4.empty_grid(label, domain, np.full((y.size, x.size), value))

    monkeypatch.setattr(stage4, "fetch_day", fetch_day)
    monkeypatch.setattr(stage4, "read_grib", read_grib)
    return missing


def test_load_clips_and_bounds_each_day(fake_days, caplog):
    fake_days["days"].add(dt.date(2025, 12, 9))
    ds = stage4.load(BOLT_CREEK, "2025-12-07T00:00:00", "2025-12-10T23:59:59")
    p = ds.precipitation
    assert p.sizes["time"] == 4 and 0 < p.sizes["x"] < 20 and 0 < p.sizes["y"] < 20
    assert [str(t)[:13] for t in ds.time_bounds.values[0]] == [
        "2025-12-06T12",
        "2025-12-07T12",
    ]
    assert p.attrs["cell_methods"] == "time: sum" and p.attrs["units"] == "mm"
    assert p.attrs["long_name"] == "Stage IV 24 h precipitation"
    assert np.isnan(p.sel(time="2025-12-09")).all()
    assert float(p.sel(time="2025-12-10").max()) == pytest.approx(12.10)
    assert ds.attrs["stage4:domain"] == "conus"
    assert ds.attrs["start_datetime"] == "2025-12-06T12:00:00Z"
    assert p.rio.crs is not None
    assert "1 of 4 daily grid(s) are missing" in caplog.text


@pytest.mark.parametrize(
    "accumulation, hours, steps",
    [
        (stage4.Accumulation.hourly, 1, 24),
        (stage4.Accumulation.sixhourly, 6, 4),
    ],
)
def test_load_sub_daily_is_stamped_at_period_end(fake_days, accumulation, hours, steps):
    ds = stage4.load(
        BOLT_CREEK, "2025-12-07T00:00:00", "2025-12-07T23:59:59", accumulation
    )
    p = ds.precipitation
    assert p.sizes["time"] == steps
    first, last = (pd.Timestamp(t) for t in ds.time.values[[0, -1]])
    assert first == pd.Timestamp("2025-12-07T00") and last.hour == 24 - hours
    bounds = ds.time_bounds.values[0]
    assert pd.Timestamp(bounds[1]) - pd.Timestamp(bounds[0]) == pd.Timedelta(
        hours=hours
    )
    assert pd.Timestamp(bounds[1]) == first
    assert p.attrs["long_name"] == f"Stage IV {hours} h precipitation"
    assert float(p.isel(time=-1).max()) == pytest.approx(24 - hours + 0.07)


def test_load_date_only_end_includes_the_day(fake_days):
    ds = stage4.load(BOLT_CREEK, "2025-12-07", "2025-12-07", "hourly")
    assert ds.sizes["time"] == 24


def test_load_hourly_without_data_over_the_aoi_warns(fake_days, caplog):
    # Before 2018, hourly grids have no data over the Northwest
    fake_days["empty"].update(pd.date_range("2015-01-01T00", "2015-01-01T05", freq="h"))
    fake_days["ends"].add(pd.Timestamp("2015-01-01T06"))
    ds = stage4.load(BOLT_CREEK, "2015-01-01T00:00", "2015-01-01T07:00", "hourly")
    p = ds.precipitation
    assert p.sizes["time"] == 8
    assert np.isnan(p.isel(time=slice(0, 7))).all() and p.isel(time=7).notnull().any()
    assert "6 of 8 hourly step(s) have no data over the AOI" in caplog.text
    assert "western US before 2018" in caplog.text
    assert "1 of 8 hourly grid(s) are missing from the archive" in caplog.text


def test_load_alaska(fake_days):
    ds = stage4.load(ANCHORAGE, "2020-06-30", "2020-07-02")
    assert ds.attrs["stage4:domain"] == "ak"
    assert "Alaska" in ds.attrs["title"]
    six = stage4.load(ANCHORAGE, "2024-08-01", "2024-08-01", "sixhourly")
    assert six.sizes["time"] == 4


def test_load_alaska_has_no_hourly(fake_days):
    with pytest.raises(ValueError, match="no hourly Alaska grids"):
        stage4.load(ANCHORAGE, "2024-08-01", "2024-08-01", "hourly")


def test_load_alaska_before_its_archive(fake_days):
    # The runner reads a day early, so a window starting 2020-07-01 is fine
    # (above); one starting earlier is not
    with pytest.raises(ValueError, match="Alaska grids start 2020-07-01"):
        stage4.load(ANCHORAGE, "2020-06-29", "2020-07-02")


def test_load_without_any_data(fake_days):
    fake_days["days"].update({dt.date(2025, 12, 7), dt.date(2025, 12, 8)})
    with pytest.raises(ValueError, match="No Stage IV daily data"):
        stage4.load(BOLT_CREEK, "2025-12-07", "2025-12-08")
