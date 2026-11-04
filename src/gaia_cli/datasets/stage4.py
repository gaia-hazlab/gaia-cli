"""NCEP Stage IV precipitation (~4.8 km): daily, 6-hourly or hourly, clipped to an AOI.

CONUS from 2002, Alaska from July 2020 (no hourly); the domain follows the AOI.
Hourly grids have no data over parts of the western US before 2018; those
steps are NaN, with a warning. Grids are
fetched from the NCAR GDEX archive (d507005) into a local cache, by default
./.cache/gaia/stage4 under the directory gaia runs in (set STAGE4_CACHE to
change it). See dev/stage4.md.

Examples
--------
    $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
    $ gaia stage stage4 -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-stage4
    $ gaia stage stage4 -i $AOI -s 2025-12-08 -e 2025-12-11 -a hourly -o skagit-stage4-hourly
"""

from __future__ import annotations

import datetime as dt
import gzip
import io as _io
import json
import logging
import os
import re
import tarfile
from enum import Enum
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import rioxarray
import unlzw3
import xarray as xr
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .. import io, provenance, window

log = logging.getLogger(__name__)

ARCHIVE_URL = "https://data.gdex.ucar.edu/d507005/stage4"
DATASET_URL = "https://gdex.ucar.edu/datasets/d507005/"

# Polar stereographic on a sphere, shared by the CONUS (HRAP) and Alaska grids
CRS = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-105 +x_0=0 +y_0=0 +R=6371200 +units=m +no_defs"
DX = 4762.5

# Cell centers of the first column and row (top), as GDAL reads the GRIB2 files.
# GRIB1 files store DX rounded to whole metres (4762 or 4763), so their
# coordinates drift by up to a few hundred metres; every file gets these.
GRIDS = {
    "conus": {
        "x0": -1902531.6743569847,
        "y0": -3426604.93005048,
        "nx": 1121,
        "ny": 881,
    },
    "ak": {"x0": -3288463.388124831, "y0": -1007361.1363076605, "nx": 460, "ny": 530},
}
# Alaska grids are archived from 2020-07-01
AK_START = pd.Timestamp("2020-07-01")

# Written by hand: the grids' outlines in EPSG:4326 (from GRIDS), the CONUS
# archive's first day, and its union with Alaska. The domain test in load() is
# exact, in the grids' CRS; this bbox check is only the runner's coarse one.
EXTENT = {
    "spatial": {
        "bbox": [
            [-178.02, 19.78, -59.91, 75.72],
            [-134.09, 19.78, -59.91, 57.87],  # CONUS
            [-178.02, 45.81, -122.31, 75.72],  # Alaska
        ]
    },
    "temporal": {"interval": [["2002-01-01T00:00:00Z", None]]},
}
KIND = "raster"
# Each accumulation is its own group, stage4/<accumulation>
VARIANT = "accumulation"
# `gaia align` default: area-weighted, so totals are preserved
REGRID_METHOD = "conservative"


class Accumulation(str, Enum):
    """Stage IV accumulation periods."""

    daily = "daily"
    sixhourly = "sixhourly"
    hourly = "hourly"


# Hours per accumulation: the archive's 24h, 06h and 01h files
HOURS = {Accumulation.daily: 24, Accumulation.sixhourly: 6, Accumulation.hourly: 1}


def domain(aoi: gpd.GeoDataFrame) -> str:
    """The Stage IV grid (``conus`` or ``ak``) that contains all of *aoi*."""
    minx, miny, maxx, maxy = aoi.to_crs(CRS).total_bounds
    for name, g in GRIDS.items():
        x0, x1 = g["x0"] - DX / 2, g["x0"] + (g["nx"] - 0.5) * DX
        y1, y0 = g["y0"] + DX / 2, g["y0"] - (g["ny"] - 0.5) * DX
        if x0 <= minx and maxx <= x1 and y0 <= miny and maxy <= y1:
            return name
    raise ValueError(
        "AOI is not within the Stage IV CONUS or Alaska grid "
        f"(bounds {[round(float(v), 4) for v in aoi.to_crs(4326).total_bounds]})"
    )


def period_ends(day: dt.date, hours: int) -> list:
    """End times of the *hours* grids in *day*'s daily tar.

    A daily tar holds the grids ending 00-23 UTC on that date: one 24 h grid
    (ending 12 UTC), 6 h grids ending 00, 06, 12 and 18 UTC, and 24 1 h grids.
    """
    first = 12 if hours == 24 else 0
    return [pd.Timestamp(day) + pd.Timedelta(hours=h) for h in range(first, 24, hours)]


def cache_dir() -> Path:
    # Under the working directory, not ~: home directories on shared HPC
    # systems are small
    return Path(os.environ.get("STAGE4_CACHE", Path(".cache") / "gaia" / "stage4"))


def member_names(end: dt.datetime, domain: str, hours: int = 24) -> list[str]:
    """Names the *hours* grid ending at *end* has had in the archive.

    CONUS: ``ST4.*.Z`` (GRIB1, to ~2014), ``ST4.*.gz`` (GRIB1, to mid-2020),
    ``st4_conus.*.grb2``; Alaska: ``st4_ak.*.gz`` (GRIB1, July 2020),
    ``st4_ak.*.grb2``.
    """
    stamp = f"{end:%Y%m%d%H}.{hours:02d}h"
    if domain == "ak":
        return [f"st4_ak.{stamp}.grb2", f"st4_ak.{stamp}.gz"]
    return [f"st4_conus.{stamp}.grb2", f"ST4.{stamp}.gz", f"ST4.{stamp}.Z"]


# (connect, read) seconds; a stalled read is retried rather than waited out
TIMEOUT = (10, 60)


def new_session() -> requests.Session:
    """A session that retries GDEX's transient failures."""
    retry = Retry(
        total=4,
        backoff_factor=2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("HEAD", "GET"),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


class _HTTPRangeFile(_io.RawIOBase):
    """A read-only, seekable file over HTTP range requests, for ``tarfile``."""

    def __init__(self, url: str, session: requests.Session):
        self.url, self.session, self.pos = url, session, 0
        # Bytes fetched beyond a read and kept for the next ones; 0 reads
        # exactly what is asked (one request per tar header)
        self.readahead = 0
        self._buf_start, self._buf = 0, b""
        r = session.head(url, timeout=TIMEOUT)
        if r.status_code == 404:
            raise FileNotFoundError(url)
        r.raise_for_status()
        self.size = int(r.headers["Content-Length"])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: 0, 1: self.pos, 2: self.size}[whence] + offset
        return self.pos

    def read(self, n=-1):
        end = self.size if n is None or n < 0 else min(self.pos + n, self.size)
        if end <= self.pos:
            return b""
        buf_end = self._buf_start + len(self._buf)
        if not (self._buf_start <= self.pos and end <= buf_end):
            fetch_end = min(max(end, self.pos + self.readahead), self.size)
            self._buf_start = self.pos
            self._buf = range_get(
                self.session, self.url, self.pos, fetch_end - self.pos
            )
        data = self._buf[self.pos - self._buf_start : end - self._buf_start]
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def range_get(session: requests.Session, url: str, offset: int, size: int) -> bytes:
    r = session.get(
        url, headers={"Range": f"bytes={offset}-{offset + size - 1}"}, timeout=TIMEOUT
    )
    r.raise_for_status()
    if r.status_code != 206:
        raise OSError(f"{url}: server ignored the range request")
    return r.content


def month_index(month: str, session: requests.Session) -> dict[str, list[int]]:
    """``{daily tar name: [offset, size]}`` for ``stage4.<month>.tar``.

    Built by reading only the tar's headers, then cached as JSON. Some months
    (2019-01 to 2019-03) also hold every grid as a loose file next to the
    daily tars; their many small headers are read 1 MB at a time.
    """
    path = cache_dir() / "index" / f"stage4.{month}.json"
    if path.exists():
        return json.loads(path.read_text())
    url = f"{ARCHIVE_URL}/stage4.{month}.tar"
    log.info("Stage IV: indexing %s", url)
    days = pd.Period(f"{month[:4]}-{month[4:]}").days_in_month
    f = _HTTPRangeFile(url, session)
    index = {}
    with tarfile.open(fileobj=f, mode="r:") as tar:
        for m in tar:
            if not m.isfile():
                continue
            name = Path(m.name).name
            if re.fullmatch(r"ST4\.\d{8}", name):
                index[name] = [m.offset_data, m.size]
                if len(index) == days:
                    break
            else:
                f.readahead = 1 << 20
    _write_atomic(path, json.dumps(index).encode())
    return index


def fetch_day(
    day: dt.date,
    domain: str,
    hours: int = 24,
    session: requests.Session | None = None,
) -> dict[pd.Timestamp, Path]:
    """Local paths of *day*'s *hours* grids, by end time, fetched if needed.

    Grids the day's tar lacks are left out. Raises ``FileNotFoundError`` if
    the archive has no tar for *day*.
    """
    year_dir = cache_dir() / f"{day:%Y}"
    names = {end: member_names(end, domain, hours) for end in period_ends(day, hours)}

    def cached() -> dict[pd.Timestamp, Path]:
        found = {}
        for end, candidates in names.items():
            for name in candidates:
                if (local := year_dir / _decompressed(name)).exists():
                    found[end] = local
                    break
        return found

    if len(found := cached()) == len(names):
        return found

    session = session or new_session()
    index = month_index(f"{day:%Y%m}", session)
    daily = f"ST4.{day:%Y%m%d}"
    if daily not in index:
        raise FileNotFoundError(f"{daily} is not in the Stage IV archive")
    offset, size = index[daily]
    url = f"{ARCHIVE_URL}/stage4.{day:%Y%m}.tar"
    log.info("Stage IV: fetching %s (%.1f MB)", daily, size / 1e6)
    with tarfile.open(
        fileobj=_io.BytesIO(range_get(session, url, offset, size))
    ) as tar:
        members = {Path(m.name).name: m for m in tar if m.isfile()}
        for candidates in names.values():
            for name in candidates:
                if name in members:
                    f = tar.extractfile(members[name])
                    assert f is not None
                    local = year_dir / _decompressed(name)
                    _write_atomic(local, _decompress(name, f.read()))
                    break
    return cached()


def _decompressed(name: str) -> str:
    return name.removesuffix(".gz").removesuffix(".Z")


def _decompress(name: str, data: bytes) -> bytes:
    if name.endswith(".gz"):
        return gzip.decompress(data)
    if name.endswith(".Z"):
        return unlzw3.unlzw(data)
    return data


def _write_atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_bytes(data)
    tmp.replace(path)


def grid(domain: str) -> tuple[np.ndarray, np.ndarray]:
    """Cell-center ``x`` and ``y`` (top to bottom) of *domain*'s grid."""
    g = GRIDS[domain]
    return g["x0"] + DX * np.arange(g["nx"]), g["y0"] - DX * np.arange(g["ny"])


def read_grib(path: Path, label: dt.datetime, domain: str) -> xr.DataArray:
    """One grid as ``precipitation`` (mm), at time *label*, on *domain*'s grid."""
    da = rioxarray.open_rasterio(path, masked=True)
    assert isinstance(da, xr.DataArray)
    da = da.squeeze("band", drop=True).load()
    x, y = grid(domain)
    if da.shape != (y.size, x.size):
        raise ValueError(f"{path}: grid {da.shape} isn't the Stage IV {domain} grid")
    return empty_grid(label, domain, da.values)


def empty_grid(label: dt.datetime, domain: str, values=None) -> xr.DataArray:
    """*values* (NaN if None) as one time step of ``precipitation`` on *domain*'s grid."""
    x, y = grid(domain)
    if values is None:
        values = np.full((y.size, x.size), np.nan)
    da = xr.DataArray(
        np.asarray(values, dtype="float32")[np.newaxis],
        dims=("time", "y", "x"),
        coords={"time": [pd.Timestamp(label)], "y": y, "x": x},
        name="precipitation",
        attrs={
            "units": "mm",
            "standard_name": "lwe_thickness_of_precipitation_amount",
        },
    )
    da = da.rio.write_crs(CRS)
    return da.rio.write_nodata(np.nan, encoded=False)


def source_attrs(domain: str) -> dict:
    region = {"conus": "conterminous US", "ak": "Alaska"}[domain]
    return {
        "title": f"NCEP Stage IV precipitation ({region})",
        "description": "NCEP/EMC Stage IV quantitative precipitation estimates: "
        "radar and gauge analyses mosaicked from the River Forecast Centers, "
        f"~4.8 km polar stereographic grid, {region}.",
        # GDEX's license for d507005 ("NCEP/CPC Four Kilometer Precipitation Set")
        "license": "CC-BY-4.0",
        "providers": [
            {
                "name": "NOAA NCEP Environmental Modeling Center",
                "roles": ["producer"],
                "url": "https://www.emc.ncep.noaa.gov/mmb/ylin/pcpanl/stage4/",
            },
            {
                "name": "NSF NCAR Geoscience Data Exchange",
                "roles": ["host", "licensor"],
                "url": DATASET_URL,
            },
        ],
        "sci:doi": "10.5065/D69Z93M3",
        "stage4:domain": domain,
        "links": [
            provenance.link("via", DATASET_URL, type="text/html"),
            provenance.link(
                "license",
                "https://creativecommons.org/licenses/by/4.0/",
                type="text/html",
            ),
        ],
    }


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    accumulation: Accumulation = Accumulation.daily,
) -> xr.Dataset:
    """Load Stage IV grids clipped to *aoi* on their native grid.

    The grid (CONUS or Alaska) is the one containing the AOI. The daily
    value for date D is the 24 h total ending D 12:00 UTC, i.e. it covers
    ``(D-1 12:00, D 12:00]`` UTC, like a PRISM day. 6-hourly and hourly
    values are stamped at the end of their period, like HRRR: the value at
    T covers ``(T-6 h, T]`` or ``(T-1 h, T]`` (see dev/stage4.md). Grids
    missing from the archive, and hourly grids without data over the AOI
    (parts of the western US before 2018), are NaN, with a warning.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    accumulation : Accumulation
        Accumulation period (daily, sixhourly or hourly), each written to its
        own group, e.g. stage4/daily. Alaska has no hourly grids.

    Returns
    -------
    xarray.Dataset
        ``precipitation`` (mm) plus CF ``time_bounds``, with provenance
        attributes.
    """
    dom = domain(aoi)
    accumulation = Accumulation(accumulation)
    hours = HOURS[accumulation]
    if dom == "ak" and accumulation == Accumulation.hourly:
        raise ValueError("Stage IV has no hourly Alaska grids; use daily or sixhourly")
    start, end = pd.Timestamp(start_date), pd.Timestamp(end_date)
    if window.DATE_ONLY.fullmatch(str(end_date)):  # a date includes its whole day
        end = end + pd.Timedelta(days=1) - pd.Timedelta(seconds=1)
    assert isinstance(start, pd.Timestamp) and isinstance(end, pd.Timestamp)
    # The runner reads a day beyond each end of the requested window
    if dom == "ak" and start + window.PAD < AK_START:
        raise ValueError(
            f"Stage IV Alaska grids start {AK_START:%Y-%m-%d}; "
            f"the window starts {start + window.PAD:%Y-%m-%d}"
        )

    # (end of period, time label) of each step: days are labelled like PRISM
    # (D 00:00 for the day ending D 12 UTC), shorter periods by their end
    if hours == 24:
        days = pd.date_range(start.normalize(), end.normalize(), freq="D")
        steps = [(d + pd.Timedelta(hours=12), d) for d in days]
    else:
        freq = f"{hours}h"
        ends = pd.date_range(start.ceil(freq), end.floor(freq), freq=freq)
        steps = [(t, t) for t in ends]
    if not steps:
        raise ValueError(
            f"No {accumulation.value} Stage IV steps in {start_date} to {end_date}"
        )

    session = new_session()
    paths: dict[pd.Timestamp, Path] = {}
    for day in sorted({e.date() for e, _ in steps}):
        try:
            paths |= fetch_day(day, dom, hours, session)
        except FileNotFoundError as e:
            log.warning("Stage IV: %s", e)

    parts, missing = [], []
    for period_end, label in steps:
        if period_end in paths:
            da = read_grib(paths[period_end], label, dom)
        else:
            missing.append(period_end)
            da = empty_grid(label, dom)
        parts.append(io.clip_to_aoi(da, aoi))
    if len(missing) == len(steps):
        raise ValueError(
            f"No Stage IV {accumulation.value} data for {start_date} to {end_date}"
        )
    if missing:
        log.warning(
            "Stage IV: %d of %d %s grid(s) are missing from the archive (first "
            "ending %s); filled with NaN",
            len(missing),
            len(steps),
            accumulation.value,
            f"{missing[0]:%Y-%m-%d %H:%M}",
        )

    da = xr.concat(parts, dim="time", combine_attrs="override")
    da.attrs["long_name"] = f"Stage IV {hours} h precipitation"
    # Grids that exist but have no data anywhere over the AOI
    empty = (
        da.isnull().all(["x", "y"]).values
        & ~da.time.isin([label for e, label in steps if e in missing]).values
    )
    if empty.any():
        hint = (
            " Hourly Stage IV has no data over parts of the western US before "
            "2018 (e.g. the Northwest); daily and sixhourly cover it."
            if accumulation == Accumulation.hourly
            else ""
        )
        log.warning(
            "Stage IV: %d of %d %s step(s) have no data over the AOI and are "
            "NaN (first %s).%s",
            int(empty.sum()),
            len(steps),
            accumulation.value,
            pd.Timestamp(da.time.values[empty.argmax()]).strftime("%Y-%m-%d %H:%M"),
            hint,
        )
    log.info(
        "Stage IV %s %s: %d step(s) x %d x %d cells",
        dom,
        accumulation.value,
        da.sizes["time"],
        da.sizes["y"],
        da.sizes["x"],
    )
    log.info(
        "Stage IV source grids are kept locally in %s; feel free to delete it "
        "(they are re-fetched as needed). Set STAGE4_CACHE to use another directory.",
        cache_dir().resolve(),
    )
    if hours == 24:
        bounds, period = (-12, 12), "24 h grids ending 12 UTC"
    else:
        bounds, period = (-hours, 0), f"{hours} h grids, stamped at their end"
    out = io.with_time_bounds(
        da, start_hours=bounds[0], end_hours=bounds[1], cell_methods="time: sum"
    )
    out.attrs = {
        **source_attrs(dom),
        "processing:lineage": f"{period}; clipped to the AOI on the native grid; "
        "9999 and grids missing from the archive → NaN",
        **provenance.time_extent(out),
    }
    return out
