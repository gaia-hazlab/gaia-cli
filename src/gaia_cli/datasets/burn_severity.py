"""USFS BAER soil burn severity map, clipped to an AOI.

See dev/burn-severity.md.
"""

from __future__ import annotations

import json
import logging
import math

import geopandas as gpd
import pyproj
import requests
import rioxarray
import shapely
import xarray as xr
from rasterio.io import MemoryFile

from .. import io, provenance

log = logging.getLogger(__name__)

# Written by hand (no STAC Collection): the image service's extent in EPSG:4326
# and its years. Static: the fire is chosen by the AOI and options, not a window.
# See dev/burn-severity.md.
EXTENT = {
    "spatial": {"bbox": [[-124.47, 31.44, -76.91, 50.03]]},
    "temporal": {"interval": [["2017-01-01T00:00:00Z", None]]},
}
KIND = "raster"
# Categorical: `gaia align` uses mode for flag variables whatever this says
REGRID_METHOD = "nearest"

SERVICE = (
    "https://imagery.geoplatform.gov/iipp/rest/services/Fire_Aviation/"
    "USFS_EDW_BAER_SoilBurnSeverityClassification/ImageServer"
)
TIMEOUT = 120
# Soil burn severity classes: code -> (name, color), the BAER legend's colors
CLASSES = {
    1: ("unburned_to_very_low", "#006400"),
    2: ("low", "#7fffd4"),
    3: ("moderate", "#ffff00"),
    4: ("high", "#ff0000"),
}
# Every other value (15 from the service, 0 in some source files) is outside the
# mapped area; stored as the fill value
FILL = 15


def _get(url: str, **params) -> requests.Response:
    r = requests.get(url, params=params, timeout=TIMEOUT)
    r.raise_for_status()
    return r


def _json(url: str, **params) -> dict:
    d = _get(url, f="json", **params).json()
    if "error" in d:
        raise RuntimeError(f"{url}: {d['error'].get('message', d['error'])}")
    return d


def find_fires(aoi: gpd.GeoDataFrame) -> list[dict]:
    """The service's fires whose footprint intersects *aoi*: ``objectid``,
    ``name`` and ``year``, newest first."""
    geom = aoi.to_crs("EPSG:4326").union_all()
    d = _json(
        f"{SERVICE}/query",
        geometry=",".join(str(v) for v in geom.bounds),
        geometryType="esriGeometryEnvelope",
        inSR=4326,
        spatialRel="esriSpatialRelIntersects",
        where="category = 1",  # primary rasters, not the service's overviews
        outFields="objectid,name,beginyear",
        returnGeometry="true",
        outSR=4326,
    )
    fires = [
        {
            "objectid": f["attributes"]["objectid"],
            "name": f["attributes"]["name"],
            "year": f["attributes"]["beginyear"],
        }
        for f in d["features"]
        # the query's bbox is coarse; keep footprints that meet the AOI itself
        if shapely.union_all(
            [shapely.Polygon(r) for r in f["geometry"]["rings"]]
        ).intersects(geom)
    ]
    return sorted(fires, key=lambda f: (-f["year"], f["name"]))


def _describe(fires: list[dict]) -> str:
    return ", ".join(f"{f['name']} ({f['year']})" for f in fires)


def select_fire(
    fires: list[dict], fire: str | None = None, year: int | None = None
) -> dict:
    """The one fire in *fires* matching *fire* (part of its name, case and
    spaces ignored) and *year*; raise if none or several do."""
    matches = [
        f
        for f in fires
        if (fire is None or fire.lower().replace(" ", "_") in f["name"].lower())
        and (year is None or f["year"] == year)
    ]
    if not fires:
        raise ValueError(
            "No BAER soil burn severity map intersects the AOI (the service "
            "holds Forest Service assessments from 2017 on)"
        )
    if not matches:
        raise ValueError(
            f"No fire matching fire={fire!r}, year={year} in the AOI; "
            f"it has {_describe(fires)}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"{len(matches)} fires intersect the AOI: {_describe(matches)}. "
            "Choose one with fire (part of its name) and/or year"
        )
    return matches[0]


def _crs(sr: dict) -> pyproj.CRS:
    """An ArcGIS spatial reference: EPSG by its latest code, else an ESRI code."""
    code = sr.get("latestWkid") or sr["wkid"]
    try:
        return pyproj.CRS.from_epsg(code)
    except pyproj.exceptions.CRSError:
        return pyproj.CRS.from_user_input(f"ESRI:{code}")


def _window(info: dict, aoi: gpd.GeoDataFrame) -> tuple[list[float], int, int]:
    """The bbox and size (columns, rows) of the cells of the raster in *info*
    that cover *aoi*, snapped to its grid."""
    e = info["extent"]
    rx, ry = info["pixelSizeX"], info["pixelSizeY"]
    minx, miny, maxx, maxy = aoi.to_crs(_crs(e["spatialReference"])).total_bounds
    ncols = round((e["xmax"] - e["xmin"]) / rx)
    nrows = round((e["ymax"] - e["ymin"]) / ry)
    c0 = min(max(math.floor((minx - e["xmin"]) / rx), 0), ncols)
    c1 = min(max(math.ceil((maxx - e["xmin"]) / rx), 0), ncols)
    r0 = min(max(math.floor((e["ymax"] - maxy) / ry), 0), nrows)
    r1 = min(max(math.ceil((e["ymax"] - miny) / ry), 0), nrows)
    bbox = [
        e["xmin"] + c0 * rx,
        e["ymax"] - r1 * ry,
        e["xmin"] + c1 * rx,
        e["ymax"] - r0 * ry,
    ]
    return bbox, c1 - c0, r1 - r0


def open_fire(objectid: int, aoi: gpd.GeoDataFrame) -> xr.DataArray:
    """The cells of the service's raster *objectid* that cover *aoi*, on the
    raster's own grid and in its CRS (the export is exact, not resampled)."""
    info = _json(f"{SERVICE}/{objectid}/info")
    bbox, width, height = _window(info, aoi)
    if not (width and height):
        raise ValueError(f"The AOI misses raster {objectid}'s extent")
    wkid = info["extent"]["spatialReference"]["wkid"]
    log.info("BAER: exporting raster %s (%d x %d cells)", objectid, width, height)
    r = _get(
        f"{SERVICE}/exportImage",
        bbox=",".join(str(v) for v in bbox),
        bboxSR=wkid,  # a WKID: WKT in the query string is refused (403)
        imageSR=wkid,
        size=f"{width},{height}",
        format="tiff",
        pixelType="U8",
        interpolation="RSP_NearestNeighbor",
        mosaicRule=json.dumps(
            {"mosaicMethod": "esriMosaicLockRaster", "lockRasterIds": [objectid]}
        ),
        f="image",
    )
    if not r.headers.get("content-type", "").startswith("image/tiff"):
        raise RuntimeError(f"exportImage returned {r.headers.get('content-type')}")
    with MemoryFile(r.content) as mem, mem.open() as src:
        da = rioxarray.open_rasterio(src).load()  # ty: ignore[unresolved-attribute]
    return da.squeeze("band", drop=True)  # ty: ignore[invalid-return-type]


def source_attrs(fire: dict) -> dict:
    """Provenance attributes for the soil burn severity map of *fire*."""
    item = f"{SERVICE}/{fire['objectid']}"
    return {
        "title": f"BAER soil burn severity: {fire['name']}",
        "description": "Soil burn severity (SBS) map from a USDA Forest Service "
        "Burned Area Emergency Response (BAER) assessment: a satellite burn "
        "severity classification (BARC) checked in the field by the BAER team.",
        "license": "other",  # U.S. public domain has no SPDX identifier
        "providers": [
            {
                "name": "USDA Forest Service, BAER Imagery Support Program",
                "roles": ["producer", "licensor"],
                "url": "https://burnseverity.cr.usgs.gov/baer/",
            },
            {
                "name": "GeoPlatform imagery services",
                "roles": ["host"],
                "url": SERVICE,
            },
        ],
        "fire_name": fire["name"],
        "start_datetime": f"{fire['year']}-01-01T00:00:00Z",
        "end_datetime": f"{fire['year']}-12-31T23:59:59Z",
        "processing:lineage": "exported from the image service on the fire "
        "raster's own grid (nearest, no resampling) and clipped to the AOI",
        "links": [
            provenance.link("via", item, type="application/json"),
            provenance.link(
                "license",
                "https://www.usa.gov/government-copyright",
                type="text/html",
            ),
        ],
    }


def load(
    aoi: gpd.GeoDataFrame, fire: str | None = None, year: int | None = None
) -> xr.Dataset:
    """BAER soil burn severity for the one fire in *aoi* matching *fire* and
    *year*, clipped to *aoi* on its native grid, with provenance attributes.

    Parameters
    ----------
    fire : str
        Part of the fire's BAER name (e.g. ``bolt creek``), to choose when
        several fires intersect the AOI. Case and spaces are ignored.
    year : int
        The fire's year, to choose when several fires intersect the AOI.
    """
    chosen = select_fire(find_fires(aoi), fire, year)
    log.info("BAER: %s (%s)", chosen["name"], chosen["year"])
    raw = open_fire(chosen["objectid"], aoi)
    # Classes 1-4; everything else is unmapped. Cells outside the AOI also get FILL
    raw = raw.where(raw.isin(list(CLASSES)), FILL).rio.write_nodata(FILL)
    raw = io.clip_to_aoi(raw, aoi)
    da = raw.where(raw != FILL).astype("float32").rename("burn_severity")
    da.attrs = {
        "long_name": "Soil burn severity",
        # CF flags mark it categorical (gaia align uses mode); flag_colors for legends
        "flag_values": list(CLASSES),
        "flag_meanings": " ".join(name for name, _ in CLASSES.values()),
        "flag_colors": " ".join(color for _, color in CLASSES.values()),
    }
    # Stored as the source's bytes; read back as float with NaN outside the map
    da.encoding = {"dtype": "uint8", "_FillValue": FILL}
    ds = da.to_dataset()
    ds.attrs = source_attrs(chosen)
    return ds
