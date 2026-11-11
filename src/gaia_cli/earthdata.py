"""NASA Earthdata: ASF and CMR search and Earthdata Login reads.

ASF granules are found with asf_search, other DAACs' (NSIDC's SMAP) with a CMR
granule search, and turned into pystac Items with the properties
``stac.group`` and ``stac.select`` use. Their files are read lazily by GDAL from
the URLs the DAAC signs for them, with the Earthdata Login token used only to
get those (see dev/opera-rtc.md).
"""

from __future__ import annotations

import datetime as dt
import logging
import os
from collections.abc import Callable, Hashable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import asf_search
import geopandas as gpd
import pandas as pd
import pystac
import rasterio
from shapely.geometry import box, mapping, shape

from . import fetch, stac

log = logging.getLogger(__name__)
# asf_search logs each search's options, the whole AOI WKT, at INFO
logging.getLogger("asf_search").setLevel(logging.WARNING)

# CMR's granule search, with results in UMM-G
CMR = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
CMR_PAGE_SIZE = 2000  # the most CMR returns per page

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
    """Point *items*' *assets* at the CloudFront URLs their DAAC signs for them.

    ASF and NSIDC answer an authenticated request with a redirect to a URL
    signed for an hour, so GDAL then reads anonymously, like Planetary
    Computer's SAS-signed hrefs (``stac.sign``). The store must be written
    within that hour.
    """
    headers = {"Authorization": f"Bearer {token()}"}

    def resolve(href: str) -> str:
        resp = fetch.get(href, headers=headers, allow_redirects=False)
        location = resp.headers.get("location", "")
        if resp.status_code not in (302, 303) or "urs.earthdata" in location:
            raise RuntimeError(
                f"could not sign {href} (HTTP {resp.status_code}); "
                f"is {TOKEN} current? {TOKENS[TOKEN]}"
            )
        return location

    # Assets may share a file (NISAR's HDF5): one request each
    hrefs = list(dict.fromkeys(it.assets[a].href for it in items for a in assets))
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
    """asf_search results intersecting *aoi*'s bounds between *start* and *end*
    (inclusive, UTC unless given), as GeoJSON features; *filters* are
    asf_search keywords (``dataset``, ``processingLevel``, ``relativeOrbit``).

    The bounds are what odc-stac reads before clipping to *aoi*, and spare
    asf_search simplifying and reorienting the polygon (with a warning)."""
    bounds = box(*aoi.to_crs("EPSG:4326").total_bounds)
    kwargs: dict[str, Any] = {"intersectsWith": bounds.wkt}
    if start:
        kwargs["start"] = _utc(start)
    if end:
        kwargs["end"] = _utc(end)
    features = asf_search.search(**kwargs, **filters).geojson()["features"]
    log.info("ASF %s: %d granule(s)", filters.get("processingLevel", ""), len(features))
    return features


def cmr_search(
    concept_id: str,
    aoi: gpd.GeoDataFrame,
    start: str | None = None,
    end: str | None = None,
) -> list[dict]:
    """CMR granules of the collection *concept_id* intersecting *aoi*'s bounds
    between *start* and *end* (inclusive, UTC unless given), as UMM-G records
    (each ``{"meta": ..., "umm": ...}``)."""
    west, south, east, north = aoi.to_crs("EPSG:4326").total_bounds
    params: dict[str, Any] = {
        "concept_id": concept_id,
        "bounding_box": f"{west},{south},{east},{north}",
        "page_size": CMR_PAGE_SIZE,
    }
    if start or end:
        params["temporal"] = ",".join(_cmr_time(t) for t in (start, end))
    records: list[dict] = []
    headers: dict[str, str] = {}
    while True:
        resp = fetch.get(CMR, params=params, headers=headers)
        resp.raise_for_status()
        page = resp.json()["items"]
        records += page
        # Later pages follow the search-after token CMR sends with each page
        after = resp.headers.get("CMR-Search-After")
        if len(page) < CMR_PAGE_SIZE or not after:
            break
        headers = {"CMR-Search-After": after}
    log.info("CMR %s: %d granule(s)", concept_id, len(records))
    return records


def _cmr_time(t: str | None) -> str:
    """*t* as CMR's temporal bound: UTC, to the second; empty is open."""
    return pd.Timestamp(_utc(t)).strftime("%Y-%m-%dT%H:%M:%SZ") if t else ""


def cmr_item(record: dict) -> pystac.Item:
    """A CMR granule *record* (UMM-G) as a pystac Item, with one ``data`` asset:
    the file's HTTPS URL.

    Its time is the granule's temporal range (``start_datetime`` to
    ``end_datetime``, ``datetime`` its start), its footprint the first
    bounding rectangle, and ``processing:datetime`` the granule's production
    time, by which :func:`newest` ranks processings.
    """
    umm = record["umm"]
    times = umm["TemporalExtent"]["RangeDateTime"]
    rect = umm["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"][
        "BoundingRectangles"
    ][0]
    bbox = [
        rect["WestBoundingCoordinate"],
        rect["SouthBoundingCoordinate"],
        rect["EastBoundingCoordinate"],
        rect["NorthBoundingCoordinate"],
    ]
    # GET DATA is the HTTPS download; GET DATA VIA DIRECT ACCESS the S3 one
    href = next(
        u["URL"]
        for u in umm["RelatedUrls"]
        if u["Type"] == "GET DATA" and u["URL"].startswith("https://")
    )
    item = pystac.Item(
        id=umm["GranuleUR"],
        geometry=mapping(box(*bbox)),
        bbox=bbox,
        datetime=dt.datetime.fromisoformat(times["BeginningDateTime"]),
        properties={
            "start_datetime": times["BeginningDateTime"],
            "end_datetime": times["EndingDateTime"],
            "processing:datetime": umm["DataGranule"]["ProductionDateTime"],
        },
    )
    item.add_asset("data", pystac.Asset(href, roles=["data"]))
    return item


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


def acquisition(item: pystac.Item) -> tuple:
    """The burst and time an item was acquired, whatever its processing."""
    return item.properties["opera:burst_id"], item.datetime


def newest(
    items: Sequence[pystac.Item], key: Callable[[pystac.Item], Hashable] = acquisition
) -> list[pystac.Item]:
    """*items* with only the newest processing of each *key* (by default, each
    burst acquisition; a static layer's is its burst ID)."""
    return stac.newest(items, key, rank=lambda it: it.properties["processing:datetime"])


def burst(item: pystac.Item) -> str:
    return item.properties["opera:burst_id"]


def set_proj(
    items: Sequence[pystac.Item],
    asset: str,
    key: Callable[[pystac.Item], Hashable] = burst,
) -> None:
    """Set ``proj:code``, ``proj:shape`` and ``proj:transform`` on *items* from
    one *asset* header per *key*: by default per burst (a burst ID keeps its
    grid on every date)."""
    first: dict[Hashable, pystac.Item] = {}
    for it in items:
        first.setdefault(key(it), it)

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
        it.properties.update(grids[key(it)])
