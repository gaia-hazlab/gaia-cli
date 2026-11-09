"""NASA Earthdata: ASF search and Earthdata Login reads.

Granules are found with asf_search and turned into pystac Items with the
properties ``stac.group`` and ``stac.select`` use. Their COGs are read lazily by
GDAL from the URLs ASF signs for them, with the Earthdata Login token used only
to get those (see dev/opera-rtc.md).
"""

from __future__ import annotations

import datetime as dt
import logging
import os
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import asf_search
import geopandas as gpd
import pandas as pd
import pystac
import rasterio
import requests
from shapely.geometry import shape

log = logging.getLogger(__name__)
# asf_search logs each search's options, the whole AOI WKT, at INFO
logging.getLogger("asf_search").setLevel(logging.WARNING)

TOKEN = "EARTHDATA_TOKEN"
TOKENS = {
    TOKEN: "an Earthdata Login user token (they expire after 60 days), "
    "https://urs.earthdata.nasa.gov/documentation/for_users/user_token"
}
# ASF's Earthdata Login host for the files in each S3 bucket, at /<path>/<key>
HOST = "https://cumulus.asf.earthdatacloud.nasa.gov"
BUCKETS = {"asf-cumulus-prod-opera-products": "OPERA"}
RASTER_EXT = "https://stac-extensions.github.io/raster/v1.1.0/schema.json"


def token() -> str:
    """The Earthdata Login token from ``EARTHDATA_TOKEN``."""
    if not (value := os.environ.get(TOKEN)):
        raise RuntimeError(f"{TOKEN} is not set: {TOKENS[TOKEN]}")
    return value


def sign(items: Sequence[pystac.Item], assets: Sequence[str]) -> None:
    """Point *items*' *assets* at the CloudFront URLs ASF signs for them.

    ASF answers an authenticated request with a redirect to a URL signed for an
    hour, so GDAL then reads anonymously, like Planetary Computer's SAS-signed
    hrefs (``stac.sign``). The store must be written within that hour.
    """
    headers = {"Authorization": f"Bearer {token()}"}

    def resolve(href: str) -> str:
        resp = requests.get(href, headers=headers, allow_redirects=False, timeout=60)
        location = resp.headers.get("location", "")
        if resp.status_code not in (302, 303) or "urs.earthdata" in location:
            raise RuntimeError(
                f"could not sign {href} (HTTP {resp.status_code}); "
                f"is {TOKEN} current? {TOKENS[TOKEN]}"
            )
        return location

    hrefs = [it.assets[a].href for it in items for a in assets]
    with ThreadPoolExecutor(16) as pool:
        signed = dict(zip(hrefs, pool.map(resolve, hrefs), strict=True))
    for it in items:
        for a in assets:
            it.assets[a].href = signed[it.assets[a].href]


def _utc(t: str) -> str:
    ts = pd.Timestamp(t)
    ts = ts.tz_localize("UTC") if ts.tzinfo is None else ts.tz_convert("UTC")
    return ts.isoformat()


def search(
    aoi: gpd.GeoDataFrame,
    start: str | None = None,
    end: str | None = None,
    **filters: Any,
) -> list[dict]:
    """asf_search results intersecting *aoi* between *start* and *end*
    (inclusive, UTC unless given), as GeoJSON features; *filters* are
    asf_search keywords (``dataset``, ``processingLevel``, ``relativeOrbit``)."""
    kwargs: dict[str, Any] = {
        "intersectsWith": aoi.to_crs("EPSG:4326").geometry.union_all().wkt
    }
    if start:
        kwargs["start"] = _utc(start)
    if end:
        kwargs["end"] = _utc(end)
    features = asf_search.search(**kwargs, **filters).geojson()["features"]
    log.info("ASF %s: %d granule(s)", filters.get("processingLevel", ""), len(features))
    return features


def to_item(feature: dict, bands: dict[str, dict] | None = None) -> pystac.Item:
    """An asf_search *feature* as a pystac Item, one asset per COG.

    Assets are named by the file's layer, lowercased (``vv``, ``mask``,
    ``local_incidence_angle``), with hrefs on :data:`HOST`. *bands* gives
    assets a ``raster:bands`` entry (``data_type``, ``nodata``) by name.
    Static layers have no acquisition time: their time is their processing.
    """
    p = feature["properties"]
    item = pystac.Item(
        id=p["sceneName"],
        geometry=feature["geometry"],
        bbox=list(shape(feature["geometry"]).bounds),
        datetime=dt.datetime.fromisoformat(p["startTime"] or p["processingDate"]),
        properties={
            "platform": p["platform"],
            "sat:orbit_state": p["flightDirection"].lower(),
            "sat:relative_orbit": int(p["pathNumber"]),
            "sat:absolute_orbit": int(p["orbit"]),
            "opera:burst_id": p["operaBurstID"],
            "processing:datetime": p["processingDate"],
        },
    )
    if bands:
        item.stac_extensions.append(RASTER_EXT)
    for url in p["s3Urls"]:
        bucket, _, key = url.removeprefix("s3://").partition("/")
        if bucket not in BUCKETS or not key.endswith(".tif"):
            continue
        name = key.rsplit("/", 1)[-1].removesuffix(".tif")
        name = name.removeprefix(f"{p['sceneName']}_").lower()
        asset = pystac.Asset(
            f"{HOST}/{BUCKETS[bucket]}/{key}",
            media_type=pystac.MediaType.COG,
            roles=["data"],
        )
        if bands and name in bands:
            asset.extra_fields["raster:bands"] = [bands[name]]
        item.add_asset(name, asset)
    return item


def newest(items: Sequence[pystac.Item]) -> list[pystac.Item]:
    """*items* with only the newest processing of each burst acquisition."""
    best: dict[tuple, pystac.Item] = {}
    for it in sorted(items, key=lambda it: it.properties["processing:datetime"]):
        best[it.properties["opera:burst_id"], it.datetime] = it
    return [
        it for it in items if best[it.properties["opera:burst_id"], it.datetime] is it
    ]


def set_proj(items: Sequence[pystac.Item], asset: str) -> None:
    """Set ``proj:code``, ``proj:shape`` and ``proj:transform`` on *items* from
    one *asset* header per burst (a burst ID keeps its grid on every date)."""
    first: dict[str, pystac.Item] = {}
    for it in items:
        first.setdefault(it.properties["opera:burst_id"], it)

    def read(item: pystac.Item) -> dict:
        href = item.assets[asset].href
        try:
            with (
                rasterio.Env(
                    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                    CPL_VSIL_CURL_USE_HEAD="NO",  # see io.py
                ),
                rasterio.open(href) as src,
            ):
                return {
                    "proj:code": src.crs.to_string(),
                    "proj:shape": list(src.shape),
                    "proj:transform": list(src.transform),
                }
        except rasterio.errors.RasterioIOError as e:
            raise RuntimeError(f"could not read {href.split('?')[0]} ({e})") from None

    with ThreadPoolExecutor(8) as pool:
        grids = dict(zip(first, pool.map(read, first.values()), strict=True))
    for it in items:
        it.properties.update(grids[it.properties["opera:burst_id"]])
