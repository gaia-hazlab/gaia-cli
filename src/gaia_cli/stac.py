"""Search and load STAC items from Planetary Computer.

Items are found with pystac-client and their assets signed with an anonymous
SAS token, fetched once per collection and reused until it nears expiry. Rasters
load lazily with odc-stac on the items' native grid (see dev/sentinel1-rtc.md).
"""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Callable, Sequence

import geopandas as gpd
import odc.stac
import pandas as pd
import pystac
import pystac_client
import requests
import xarray as xr

from . import io
from .window import Window

log = logging.getLogger(__name__)

PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
PC_SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/token"
# A token is renewed when it expires sooner than this: lazy reads happen while
# the store is written, after load() returns
TOKEN_MARGIN = pd.Timedelta(hours=1)

# GDAL settings for efficient remote reading
odc.stac.configure_rio(cloud_defaults=True)

# collection -> (token, expiry)
_tokens: dict[str, tuple[str, pd.Timestamp]] = {}


def token(collection: str) -> str:
    """An anonymous SAS token for *collection*'s assets."""
    cached = _tokens.get(collection)
    if cached and cached[1] - pd.Timestamp.now(tz="UTC") > TOKEN_MARGIN:
        return cached[0]
    resp = requests.get(f"{PC_SAS}/{collection}", timeout=30)
    resp.raise_for_status()
    doc = resp.json()
    sas, expiry = str(doc["token"]), pd.Timestamp(doc["msft:expiry"])
    assert isinstance(expiry, pd.Timestamp)  # not NaT
    _tokens[collection] = (sas, expiry)
    return sas


def sign(items: Sequence[pystac.Item], collection: str) -> list[pystac.Item]:
    """Copies of *items* whose blob-storage asset hrefs carry a SAS token."""
    sas = token(collection)
    signed = []
    for item in items:
        item = item.clone()
        for asset in item.assets.values():
            if ".blob.core.windows.net/" in asset.href and "?" not in asset.href:
                asset.href = f"{asset.href}?{sas}"
        signed.append(item)
    return signed


def search(
    collection: str,
    aoi: gpd.GeoDataFrame,
    start: str,
    end: str,
    query: dict | None = None,
) -> list[pystac.Item]:
    """Signed items of *collection* intersecting *aoi* between *start* and *end*
    (inclusive UTC times), filtered by a STAC API *query*."""
    client = pystac_client.Client.open(PC_STAC)
    found = client.search(
        collections=[collection],
        intersects=aoi.to_crs("EPSG:4326").geometry.union_all(),
        datetime=f"{start}/{end}",
        query=query,
    )
    items = list(found.items())
    log.info("%s: %d item(s) from %s to %s", collection, len(items), start, end)
    return sign(items, collection)


def item_crs(item: pystac.Item) -> str:
    """*item*'s CRS from ``proj:code`` or ``proj:epsg``."""
    return item.properties.get("proj:code") or f"EPSG:{item.properties['proj:epsg']}"


def native_crs(items: Sequence[pystac.Item]) -> str:
    """The CRS most *items* are in."""
    return Counter(map(item_crs, items)).most_common(1)[0][0]


def load(
    items: Sequence[pystac.Item],
    bands: Sequence[str],
    aoi: gpd.GeoDataFrame,
    groupby: str | Callable = "id",
) -> xr.Dataset:
    """*bands* of *items*, lazily, clipped to *aoi* on the native grid.

    Items in another CRS than most are reprojected onto that one; *groupby*
    (an odc-stac groupby) decides which items share a time step.
    """
    crs = native_crs(items)
    native = [it for it in items if item_crs(it) == crs]
    if others := len(items) - len(native):
        log.warning("%d of %d item(s) reprojected to %s", others, len(items), crs)
    # The pixel size of the items already in that CRS (proj:transform's a)
    transform = native[0].properties.get("proj:transform")
    ds = odc.stac.load(
        items,
        bands=bands,
        crs=crs,
        resolution=abs(transform[0]) if transform else None,
        bbox=tuple(aoi.to_crs("EPSG:4326").total_bounds),
        groupby=groupby,
        chunks={},
    )
    return io.clip_to_aoi(ds, aoi)  # ty: ignore[invalid-argument-type, invalid-return-type]


def select(
    scenes: Sequence[Sequence[pystac.Item]],
    windows: Sequence[Window],
    match: str | None = None,
) -> list[tuple[str | None, Sequence[pystac.Item]]]:
    """The scenes (each a group of items, such as a pass's frames) that
    *windows* choose, with each one's window name, in time order.

    A scene's time is its first item's. Each window keeps the scenes in
    ``[start, end)`` that its ``select`` picks (``all`` by default; ``first``,
    ``last`` or ``nearest`` to its anchor time, ``count`` of them, default 1).
    With *match*, an item property such as ``sat:relative_orbit``, the choice
    is made for each value separately, and only values with scenes in every
    window are kept: a pre and a post scene on each track that has both.
    """

    def when(scene):
        t = pd.Timestamp(scene[0].datetime)
        return t.tz_convert(None) if t.tzinfo else t

    def pick(found, w):
        if w.select in (None, "all"):
            return found
        n = w.count or 1
        if w.select == "first":
            return found[:n]
        if w.select == "last":
            return found[-n:]
        assert w.anchor is not None  # parse_window requires `time` for nearest
        return sorted(found, key=lambda s: (abs(when(s) - w.anchor), when(s)))[:n]

    ordered = sorted(scenes, key=lambda s: (when(s), s[0].id))
    chosen: dict[object, list[tuple[str | None, Sequence[pystac.Item]]]] = {}
    for w in windows:
        inside = [
            s
            for s in ordered
            if (w.start is None or when(s) >= w.start)
            and (w.end is None or when(s) < w.end)
        ]
        keys = dict.fromkeys(
            s[0].properties.get(match) if match else None for s in inside
        )
        for key in keys:
            group = [
                s for s in inside if not match or s[0].properties.get(match) == key
            ]
            chosen.setdefault(key, []).extend((w.name, s) for s in pick(group, w))
    if match and len(windows) > 1:
        names = {w.name for w in windows}
        chosen = {k: v for k, v in chosen.items() if {name for name, _ in v} == names}
    seen, out = set(), []
    for name, scene in sorted(
        (c for v in chosen.values() for c in v), key=lambda c: (when(c[1]), c[1][0].id)
    ):
        if scene[0].id not in seen:  # a scene in two overlapping windows counts once
            seen.add(scene[0].id)
            out.append((name, scene))
    return out
