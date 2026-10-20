"""STAC metadata for cached PRISM grids.

Each cached COG gets its own STAC Item (one variable, one day), written next
to it when it is cached. Items are grouped into one Catalog per variable and
year, all children of a root Collection::

    collection.json
    us/an/{res}/{var}/daily/{yyyy}/catalog.json
    us/an/{res}/{var}/daily/{yyyy}/prism_{var}_us_{code}_{yyyymmdd}.json

For loading, the per-variable items for a day are merged into one item with
an asset per variable, which ``odc.stac`` expects.

Items are stamped at 00:00 UTC on the PRISM day label (see dev/prism.md);
``start_datetime``/``end_datetime`` give the actual 1200–1200 UTC period.
"""

from __future__ import annotations

import datetime as dt
import posixpath

import pystac
from rasterio.io import MemoryFile

COLLECTION_ID = "prism-daily"
COLLECTION_KEY = "collection.json"
JSON = "application/json"
GEOJSON = "application/geo+json"
COG = "image/tiff; application=geotiff; profile=cloud-optimized"
EXTENSIONS = [
    "https://stac-extensions.github.io/projection/v2.0.0/schema.json",
    "https://stac-extensions.github.io/raster/v1.1.0/schema.json",
]
UNITS = {
    "ppt": "mm",
    "tmin": "degC",
    "tmax": "degC",
    "tmean": "degC",
    "tdmean": "degC",
    "vpdmin": "hPa",
    "vpdmax": "hPa",
}

# {variable: {date: STAC item dict}}
Entries = dict[str, dict[dt.date, dict]]


def _iso(t: dt.datetime) -> str:
    return t.astimezone(dt.UTC).isoformat().replace("+00:00", "Z")


def catalog_key(key: str) -> str:
    """Key of the variable-year catalog for the grid at *key*."""
    return f"{posixpath.dirname(key)}/catalog.json"


def grid_item(
    tif: bytes,
    meta: dict,
    base: str,
    key: str,
    updated: dt.datetime,
) -> dict:
    """Build the STAC item for one cached grid from its COG bytes.

    Parameters
    ----------
    tif : bytes
        The COG.
    meta : dict
        Metadata from ``upstream.download``.
    base : str
        Public URL of the cache root.
    key : str
        Cache key of the grid, without extension.
    updated : datetime
        When the grid was fetched into the cache.
    """
    href = f"{base}/{key}.tif"
    with MemoryFile(tif) as mem, mem.open() as src:
        left, bottom, right, top = src.bounds
        epsg = src.crs.to_epsg()
        shape = [src.height, src.width]
        transform = list(src.transform)[:6]
        band = {"data_type": src.dtypes[0], "nodata": src.nodata}

    var = meta["variable"]
    date = dt.date.fromisoformat(meta["date"])
    midnight = dt.datetime(date.year, date.month, date.day, tzinfo=dt.UTC)
    ring = [[left, bottom], [right, bottom], [right, top], [left, top], [left, bottom]]

    properties = {
        "datetime": _iso(midnight),
        "start_datetime": _iso(midnight - dt.timedelta(hours=12)),
        "end_datetime": _iso(midnight + dt.timedelta(hours=12)),
        "updated": _iso(updated),
        "proj:code": f"EPSG:{epsg}",
        "proj:shape": shape,
        "proj:transform": transform,
        "prism:variable": var,
        "prism:resolution": meta.get("resolution"),
        "prism:release": meta.get("release", 0),
        "prism:upstream_etag": meta.get("upstream_etag"),
        "prism:upstream_last_modified": meta.get("upstream_last_modified"),
    }
    # The rest of PRISM's info.txt, e.g. PRISM_DATASET_VERSION -> prism:dataset_version
    for field, value in meta.get("info", {}).items():
        name = f"prism:{field.lower().removeprefix('prism_')}"
        properties.setdefault(name, value)

    links = _item_links(base, key)
    if meta.get("upstream_url"):
        links.append(
            {"rel": "via", "href": meta["upstream_url"], "type": "application/zip"}
        )

    return {
        "type": "Feature",
        "stac_version": "1.1.0",
        "stac_extensions": EXTENSIONS,
        "id": posixpath.basename(key),
        "collection": COLLECTION_ID,
        "bbox": [left, bottom, right, top],
        "geometry": {"type": "Polygon", "coordinates": [ring]},
        "properties": properties,
        "assets": {
            var: {
                "href": href,
                "type": COG,
                "roles": ["data"],
                "raster:bands": [{**band, "unit": UNITS.get(var)}],
            }
        },
        "links": links,
    }


def _item_links(base: str, key: str) -> list[dict]:
    return [
        {"rel": "self", "href": f"{base}/{key}.json", "type": GEOJSON},
        {"rel": "root", "href": f"{base}/{COLLECTION_KEY}", "type": JSON},
        {"rel": "parent", "href": f"{base}/{catalog_key(key)}", "type": JSON},
        {"rel": "collection", "href": f"{base}/{COLLECTION_KEY}", "type": JSON},
    ]


def relink_item(item: dict, base: str, key: str) -> bool:
    """Reset an item's catalog links (keeping others, e.g. ``via``).
    Returns True if *item* changed."""
    before = repr(item)
    managed = {"self", "root", "parent", "collection"}
    others = [lk for lk in item.get("links", []) if lk["rel"] not in managed]
    item["links"] = _item_links(base, key) + others
    item["collection"] = COLLECTION_ID
    return repr(item) != before


def new_catalog(base: str, key: str, item: dict) -> dict:
    """Empty variable-year catalog at *key* for grids like *item*."""
    p = item["properties"]
    year = p["datetime"][:4]
    return {
        "type": "Catalog",
        "stac_version": "1.1.0",
        "id": f"prism-{p['prism:resolution']}-{p['prism:variable']}-{year}",
        "title": f"PRISM daily {p['prism:variable']} {p['prism:resolution']} {year}",
        "description": (
            f"Cached PRISM daily {p['prism:variable']} grids at "
            f"{p['prism:resolution']} for {year}."
        ),
        "links": [
            {"rel": "self", "href": f"{base}/{key}", "type": JSON},
            {"rel": "root", "href": f"{base}/{COLLECTION_KEY}", "type": JSON},
            {"rel": "parent", "href": f"{base}/{COLLECTION_KEY}", "type": JSON},
        ],
    }


def new_collection(base: str, item: dict) -> dict:
    """Root collection, with its extent seeded from *item*."""
    p = item["properties"]
    return {
        "type": "Collection",
        "stac_version": "1.1.0",
        "id": COLLECTION_ID,
        "title": "PRISM daily climate grids (gaia cache)",
        "description": (
            "Unzipped PRISM AN daily time series COGs for the conterminous US, "
            "copied byte-for-byte from data.prism.oregonstate.edu as they are "
            "requested. Each item records its PRISM release number "
            "(`prism:release`, 8 is final) and when it was fetched (`updated`)."
        ),
        "keywords": ["PRISM", "climate", "precipitation", "temperature"],
        "license": "other",
        "providers": [
            {
                "name": "PRISM Group, Oregon State University",
                "roles": ["producer", "licensor"],
                "url": "https://prism.oregonstate.edu",
            },
            {"name": "GAIA hazlab", "roles": ["host"]},
        ],
        "extent": {
            "spatial": {"bbox": [list(item["bbox"])]},
            "temporal": {"interval": [[p["start_datetime"], p["end_datetime"]]]},
        },
        "summaries": {"prism:variable": [], "prism:resolution": []},
        "links": [
            {"rel": "self", "href": f"{base}/{COLLECTION_KEY}", "type": JSON},
            {"rel": "root", "href": f"{base}/{COLLECTION_KEY}", "type": JSON},
            {
                "rel": "license",
                "href": "https://prism.oregonstate.edu/terms/",
                "type": "text/html",
                "title": "PRISM terms of use",
            },
        ],
    }


def add_link(doc: dict, rel: str, href: str, media_type: str = JSON) -> bool:
    """Add a link unless one with the same rel and href exists.
    Returns True if *doc* changed."""
    links = doc["links"]
    if any(lk["rel"] == rel and lk["href"] == href for lk in links):
        return False
    links.append({"rel": rel, "href": href, "type": media_type})
    # Keep catalog listings in a stable order
    links.sort(key=lambda lk: (lk["rel"] in ("child", "item"), lk["href"]))
    return True


def extend_collection(collection: dict, item: dict) -> bool:
    """Grow the collection's extent and summaries to cover *item*.
    Returns True if *collection* changed."""
    before = repr(collection)
    p = item["properties"]
    bbox = collection["extent"]["spatial"]["bbox"][0]
    bbox[:] = [
        min(bbox[0], item["bbox"][0]),
        min(bbox[1], item["bbox"][1]),
        max(bbox[2], item["bbox"][2]),
        max(bbox[3], item["bbox"][3]),
    ]
    interval = collection["extent"]["temporal"]["interval"][0]
    interval[:] = [
        min(interval[0], p["start_datetime"]),
        max(interval[1], p["end_datetime"]),
    ]
    for field in ("prism:variable", "prism:resolution"):
        values = collection["summaries"].setdefault(field, [])
        if p[field] not in values:
            values.append(p[field])
            values.sort()
    return repr(collection) != before


def release(item: dict) -> int:
    """PRISM release number of a cached grid (8 is final)."""
    return int(item["properties"].get("prism:release", 0))


def updated(item: dict) -> dt.datetime:
    """When a cached grid was last fetched."""
    return dt.datetime.fromisoformat(item["properties"]["updated"])


def items(entries: Entries) -> list[pystac.Item]:
    """Merge per-variable items into one item per day, for ``odc.stac.load``."""
    by_date: dict[dt.date, list[dict]] = {}
    for per_date in entries.values():
        for date, item in per_date.items():
            by_date.setdefault(date, []).append(item)

    out = []
    for date in sorted(by_date):
        parts = by_date[date]
        first = parts[0]
        merged = {
            **{k: v for k, v in first.items() if k != "collection"},
            "id": first["id"].replace(
                f"_{first['properties']['prism:variable']}_", "_"
            ),
            "properties": {
                k: v
                for k, v in first["properties"].items()
                if not k.startswith("prism:") and k != "updated"
            },
            "assets": {k: v for p in parts for k, v in p["assets"].items()},
            "links": [],
        }
        out.append(pystac.Item.from_dict(merged))
    return out
