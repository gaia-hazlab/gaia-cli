"""Access to the PRISM servers at Oregon State.

Two endpoints are used:

- ``data.prism.oregonstate.edu`` serves one ``.zip`` "grid package" per
  variable and day (a COG plus ancillary files).
- ``services.nacse.org`` answers "release date" queries for a whole date range
  in one request. These do not count against the download limit.

PRISM blocks a file after it is downloaded twice in 24 hours and may block IPs
for excessive activity, so downloads are serialized with a pause between them.
"""

from __future__ import annotations

import datetime as dt
import io
import threading
import time
import zipfile
from dataclasses import dataclass

import requests

DATA_URL = "https://data.prism.oregonstate.edu/time_series/us/an"
SERVICE_URL = "https://services.nacse.org/prism/data/get"

VARIABLES = ("ppt", "tmin", "tmax", "tmean", "tdmean", "vpdmin", "vpdmax")
# Resolution name used in URLs -> resolution code used in file names
RES_CODES = {"800m": "30s", "4km": "25m"}
# Daily grids are re-released up to 8 times over ~6 months; 8 is final.
FINAL_RELEASE = 8
# Pause between downloads, as in PRISM's sample bulk download script.
REQUEST_PAUSE = 2.0

_throttle_lock = threading.Lock()
_last_download = 0.0


class NotReleased(Exception):
    """The requested grid is not (yet) available upstream."""


@dataclass(frozen=True)
class Release:
    """Upstream release info for one grid."""

    date: dt.date
    release_date: str
    release: int


@dataclass(frozen=True)
class Grid:
    """A downloaded grid: COG bytes plus PRISM metadata."""

    tif: bytes
    meta: dict


def stem(var: str, date: dt.date, res: str = "800m") -> str:
    """File name without extension, e.g. ``prism_ppt_us_30s_20251201``."""
    return f"prism_{var}_us_{RES_CODES[res]}_{date:%Y%m%d}"


def zip_url(var: str, date: dt.date, res: str = "800m") -> str:
    return f"{DATA_URL}/{res}/{var}/daily/{date:%Y}/{stem(var, date, res)}.zip"


def parse_info(text: str) -> dict:
    """Parse a PRISM ``.info.txt`` file into a dict of ``KEY: value`` pairs."""
    info = {}
    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep and key.startswith("PRISM_"):
            info[key.strip()] = value.strip()
    return info


def releases(
    var: str,
    start: dt.date,
    end: dt.date,
    res: str = "800m",
    session: requests.Session | None = None,
) -> dict[dt.date, Release]:
    """Return upstream release info for each available day in ``[start, end]``.

    The service rejects end dates past the newest grid, so *end* is clamped
    to yesterday (UTC) and stepped back a few days if still rejected.
    """
    http = session or requests
    yesterday = dt.datetime.now(dt.UTC).date() - dt.timedelta(days=1)
    end = min(end, yesterday)
    for _ in range(4):
        if end < start:
            return {}
        url = (
            f"{SERVICE_URL}/releaseDate/us/{res}/{var}/"
            f"{start:%Y%m%d}/{end:%Y%m%d}?json=true"
        )
        resp = http.get(url, timeout=60)
        resp.raise_for_status()
        try:
            rows = resp.json()
        except ValueError:
            # Plain-text error such as "Invalid date: 20261014"
            end -= dt.timedelta(days=1)
            continue
        # A single-day response is a flat row rather than a list of rows
        if rows and not isinstance(rows[0], list):
            rows = [rows]
        out = {}
        for row in rows:
            if len(row) < 4 or not row[0] or not row[3]:
                continue
            day = dt.date.fromisoformat(row[0])
            out[day] = Release(day, row[1], int(row[3]))
        return out
    return {}


def _throttle():
    global _last_download
    with _throttle_lock:
        wait = _last_download + REQUEST_PAUSE - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        _last_download = time.monotonic()


def download(
    var: str,
    date: dt.date,
    res: str = "800m",
    session: requests.Session | None = None,
) -> Grid:
    """Download one grid package and return its COG and metadata.

    Raises
    ------
    NotReleased
        If the grid does not exist upstream (HTTP 404).
    """
    http = session or requests
    url = zip_url(var, date, res)
    _throttle()
    resp = http.get(url, timeout=300)
    if resp.status_code == 404:
        raise NotReleased(url)
    resp.raise_for_status()

    name = stem(var, date, res)
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        tif = zf.read(f"{name}.tif")
        info = parse_info(zf.read(f"{name}.info.txt").decode())

    meta = {
        "variable": var,
        "date": date.isoformat(),
        "resolution": res,
        "release": int(info.get("PRISM_DATASET_RELEASE_NUMBER", 0)),
        "upstream_url": url,
        "upstream_etag": resp.headers.get("ETag"),
        "upstream_last_modified": resp.headers.get("Last-Modified"),
        "info": info,
    }
    return Grid(tif=tif, meta=meta)
