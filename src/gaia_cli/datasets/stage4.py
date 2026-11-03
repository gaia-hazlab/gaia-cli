"""NCEP Stage IV daily precipitation (~4.8 km), clipped to an AOI.

CONUS from 2002, Alaska from July 2020; the domain follows the AOI. Grids are
fetched from the NCAR GDEX archive (d507005) into a local cache, by default
./.cache/gaia/stage4 under the directory gaia runs in (set STAGE4_CACHE to
change it). See dev/stage4.md.

Examples
--------
    $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
    $ gaia stage stage4 -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-stage4
"""

from __future__ import annotations

import datetime as dt
import gzip
import io as _io
import json
import logging
import os
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


def cache_dir() -> Path:
    # Under the working directory, not ~: home directories on shared HPC
    # systems are small
    return Path(os.environ.get("STAGE4_CACHE", Path(".cache") / "gaia" / "stage4"))


def member_names(day: dt.date, domain: str) -> list[str]:
    """Names the 24 h grid ending at 12 UTC on *day* has had in the archive.

    CONUS: ``ST4.*.Z`` (GRIB1, to ~2014), ``ST4.*.gz`` (GRIB1, to mid-2020),
    ``st4_conus.*.grb2``; Alaska: ``st4_ak.*.gz`` (GRIB1, July 2020),
    ``st4_ak.*.grb2``.
    """
    stamp = f"{day:%Y%m%d}12.24h"
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
        data = range_get(self.session, self.url, self.pos, end - self.pos)
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

    Built by reading only the tar's headers, then cached as JSON.
    """
    path = cache_dir() / "index" / f"stage4.{month}.json"
    if path.exists():
        return json.loads(path.read_text())
    url = f"{ARCHIVE_URL}/stage4.{month}.tar"
    log.info("Stage IV: indexing %s", url)
    with tarfile.open(fileobj=_HTTPRangeFile(url, session), mode="r:") as tar:
        index = {Path(m.name).name: [m.offset_data, m.size] for m in tar if m.isfile()}
    _write_atomic(path, json.dumps(index).encode())
    return index


def fetch_day(
    day: dt.date, domain: str, session: requests.Session | None = None
) -> Path:
    """Local path of the 24 h grid ending 12 UTC on *day*, fetched if needed.

    Raises ``FileNotFoundError`` if the archive doesn't have it.
    """
    names = member_names(day, domain)
    year_dir = cache_dir() / f"{day:%Y}"
    for name in names:
        local = year_dir / _decompressed(name)
        if local.exists():
            return local

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
        for name in names:
            if name in members:
                f = tar.extractfile(members[name])
                assert f is not None
                local = year_dir / _decompressed(name)
                _write_atomic(local, _decompress(name, f.read()))
                return local
    raise FileNotFoundError(f"{daily} has no {domain} 24 h grid ({' or '.join(names)})")


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


def read_grib(path: Path, day: dt.date, domain: str) -> xr.DataArray:
    """One 24 h grid as ``precipitation`` (mm), labelled *day*, on *domain*'s grid."""
    da = rioxarray.open_rasterio(path, masked=True)
    assert isinstance(da, xr.DataArray)
    da = da.squeeze("band", drop=True).load()
    x, y = grid(domain)
    if da.shape != (y.size, x.size):
        raise ValueError(f"{path}: grid {da.shape} isn't the Stage IV {domain} grid")
    return empty_day(day, domain, da.values)


def empty_day(day: dt.date, domain: str, values=None) -> xr.DataArray:
    """*values* (NaN if None) as one day of ``precipitation`` on *domain*'s grid."""
    x, y = grid(domain)
    if values is None:
        values = np.full((y.size, x.size), np.nan)
    da = xr.DataArray(
        np.asarray(values, dtype="float32")[np.newaxis],
        dims=("time", "y", "x"),
        coords={"time": [pd.Timestamp(day)], "y": y, "x": x},
        name="precipitation",
        attrs={
            "units": "mm",
            "standard_name": "lwe_thickness_of_precipitation_amount",
            "long_name": "Stage IV 24 h precipitation",
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
    ``(D-1 12:00, D 12:00]`` UTC, like a PRISM day (see dev/stage4.md). Days
    missing from the archive are NaN.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    accumulation : Accumulation
        Accumulation period, each written to its own group (stage4/daily).

    Returns
    -------
    xarray.Dataset
        ``precipitation`` (mm) plus CF ``time_bounds``, with provenance
        attributes.
    """
    dom = domain(aoi)
    days = pd.date_range(str(start_date)[:10], str(end_date)[:10], freq="D")
    # The runner reads a day beyond each end of the requested window
    if dom == "ak" and days[0] + window.PAD < AK_START:
        raise ValueError(
            f"Stage IV Alaska grids start {AK_START:%Y-%m-%d}; "
            f"the window starts {days[0] + window.PAD:%Y-%m-%d}"
        )

    session = new_session()
    parts, missing = [], []
    for day in (t.date() for t in days):
        try:
            da = read_grib(fetch_day(day, dom, session), day, dom)
        except FileNotFoundError as e:
            log.warning("Stage IV: %s; filling %s with NaN", e, day)
            missing.append(day)
            da = empty_day(day, dom)
        parts.append(io.clip_to_aoi(da, aoi))
    if len(missing) == len(days):
        raise ValueError(
            f"No Stage IV data for {days[0]:%Y-%m-%d} to {days[-1]:%Y-%m-%d}"
        )

    da = xr.concat(parts, dim="time", combine_attrs="override")
    log.info(
        "Stage IV %s: %d day(s) x %d x %d cells",
        dom,
        da.sizes["time"],
        da.sizes["y"],
        da.sizes["x"],
    )
    log.info(
        "Stage IV source grids are kept locally in %s; feel free to delete it "
        "(they are re-fetched as needed). Set STAGE4_CACHE to use another directory.",
        cache_dir().resolve(),
    )
    out = io.with_time_bounds(
        da, start_hours=-12, end_hours=12, cell_methods="time: sum"
    )
    out.attrs = {
        **source_attrs(dom),
        "processing:lineage": "24 h grids ending 12 UTC; clipped to the AOI on "
        "the native grid; 9999 and days missing from the archive → NaN",
        **provenance.time_extent(out),
    }
    return out
