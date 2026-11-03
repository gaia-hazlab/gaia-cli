import datetime as dt
import gzip
import io
import json
import tarfile

import geopandas as gpd
import numpy as np
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
    day = dt.date(2016, 1, 1)
    assert stage4.member_names(day, "conus") == [
        "st4_conus.2016010112.24h.grb2",
        "ST4.2016010112.24h.gz",
        "ST4.2016010112.24h.Z",
    ]
    assert stage4.member_names(day, "ak") == [
        "st4_ak.2016010112.24h.grb2",
        "st4_ak.2016010112.24h.gz",
    ]


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
            "ST4.2016010111.01h.gz": gzip.compress(b"hourly"),
            "ST4.2016010112.24h.gz": gzip.compress(b"daily grib"),
        }
    )
    other = _tar({"ST4.2016010212.24h.gz": gzip.compress(b"next day")})
    # Members are out of date order, with paths, like the archive's
    return FakeArchive(
        {"stage4.201601.tar": _tar({"./ST4.20160102": other, "./ST4.20160101": day})}
    )


def test_fetch_day_extracts_one_grid_by_range(cache):
    session = archive()
    path = stage4.fetch_day(dt.date(2016, 1, 1), "conus", session)
    assert path == cache / "2016" / "ST4.2016010112.24h"
    assert path.read_bytes() == b"daily grib"
    index = json.loads((cache / "index" / "stage4.201601.json").read_text())
    assert set(index) == {"ST4.20160101", "ST4.20160102"}

    # Cached: no more requests
    gets = session.gets
    assert stage4.fetch_day(dt.date(2016, 1, 1), "conus", session) == path
    assert session.gets == gets


def test_fetch_day_missing(cache):
    session = archive()
    with pytest.raises(FileNotFoundError, match="ST4.20160103 is not in"):
        stage4.fetch_day(dt.date(2016, 1, 3), "conus", session)
    with pytest.raises(FileNotFoundError, match="no ak 24 h grid"):
        stage4.fetch_day(dt.date(2016, 1, 1), "ak", session)
    with pytest.raises(FileNotFoundError, match="stage4.201602.tar"):
        stage4.fetch_day(dt.date(2016, 2, 1), "conus", session)


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
    """Serve each day's grid as its day of month, except *missing* days."""
    missing = set()

    def fetch_day(day, domain, session=None):
        if day in missing:
            raise FileNotFoundError(f"ST4.{day:%Y%m%d} is not in the archive")
        return day

    def read_grib(path, day, domain):
        x, y = stage4.grid(domain)
        return stage4.empty_day(day, domain, np.full((y.size, x.size), day.day))

    monkeypatch.setattr(stage4, "fetch_day", fetch_day)
    monkeypatch.setattr(stage4, "read_grib", read_grib)
    return missing


def test_load_clips_and_bounds_each_day(fake_days):
    fake_days.add(dt.date(2025, 12, 9))
    ds = stage4.load(BOLT_CREEK, "2025-12-07T00:00:00", "2025-12-10T23:59:59")
    p = ds.precipitation
    assert p.sizes["time"] == 4 and 0 < p.sizes["x"] < 20 and 0 < p.sizes["y"] < 20
    assert [str(t)[:13] for t in ds.time_bounds.values[0]] == [
        "2025-12-06T12",
        "2025-12-07T12",
    ]
    assert p.attrs["cell_methods"] == "time: sum" and p.attrs["units"] == "mm"
    assert np.isnan(p.sel(time="2025-12-09")).all()
    assert float(p.sel(time="2025-12-10").max()) == 10
    assert ds.attrs["stage4:domain"] == "conus"
    assert ds.attrs["start_datetime"] == "2025-12-06T12:00:00Z"
    assert p.rio.crs is not None


def test_load_alaska(fake_days):
    ds = stage4.load(ANCHORAGE, "2020-06-30", "2020-07-02")
    assert ds.attrs["stage4:domain"] == "ak"
    assert "Alaska" in ds.attrs["title"]


def test_load_alaska_before_its_archive(fake_days):
    # The runner reads a day early, so a window starting 2020-07-01 is fine
    # (above); one starting earlier is not
    with pytest.raises(ValueError, match="Alaska grids start 2020-07-01"):
        stage4.load(ANCHORAGE, "2020-06-29", "2020-07-02")


def test_load_without_any_data(fake_days):
    fake_days.update({dt.date(2025, 12, 7), dt.date(2025, 12, 8)})
    with pytest.raises(ValueError, match="No Stage IV data"):
        stage4.load(BOLT_CREEK, "2025-12-07", "2025-12-08")
