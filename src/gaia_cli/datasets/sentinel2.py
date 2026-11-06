"""Sentinel-2 L2A surface reflectance from Planetary Computer, clipped to an AOI.

Every pass over the AOI in the window is a time step of surface reflectance (0-1)
and the scene classification (SCL) on the scenes' UTM grid at 10 m, with its
platform, orbit and cloud cover as coordinates. See dev/sentinel2.md.

Examples
--------
    $ gaia stage sentinel2 -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-s2
    $ gaia stage sentinel2 -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-s2 -m 50 -b red green blue
"""

from __future__ import annotations

import logging
import re
from typing import Literal

import geopandas as gpd
import numpy as np
import pystac
import xarray as xr

from .. import provenance, stac
from ..window import Window

log = logging.getLogger(__name__)

COLLECTION = "sentinel-2-l2a"
COLLECTION_URL = f"{stac.PC_STAC}/collections/{COLLECTION}"
# Copied from COLLECTION_URL
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [["2015-06-27T10:25:31Z", None]]},
}
KIND = "raster"
# Takes pre/post windows with select and count (see stac.select)
SCENES = True
# `gaia align` default for reflectance; SCL (flag_values) always uses mode
REGRID_METHOD = "average"

# Band -> asset, by the bands' eo:bands common names
ASSETS = {
    "coastal": "B01",
    "blue": "B02",
    "green": "B03",
    "red": "B04",
    "rededge1": "B05",
    "rededge2": "B06",
    "rededge3": "B07",
    "nir": "B08",
    "nir08": "B8A",
    "nir09": "B09",
    "swir16": "B11",
    "swir22": "B12",
    "scl": "SCL",
}
Band = Literal[
    "coastal",
    "blue",
    "green",
    "red",
    "rededge1",
    "rededge2",
    "rededge3",
    "nir",
    "nir08",
    "nir09",
    "swir16",
    "swir22",
    "scl",
]
# Asset pixel sizes (m)
GSD = {"B02": 10, "B03": 10, "B04": 10, "B08": 10, "B01": 60, "B09": 60}
# Scene classification classes (ESA L2A product specification)
SCL_CLASSES = {
    0: "no_data",
    1: "saturated_or_defective",
    2: "topographic_or_cast_shadow",
    3: "cloud_shadow",
    4: "vegetation",
    5: "not_vegetated",
    6: "water",
    7: "unclassified",
    8: "cloud_medium_probability",
    9: "cloud_high_probability",
    10: "thin_cirrus",
    11: "snow_or_ice",
}
# Digital numbers are reflectance x 10000, plus 1000 from processing baseline
# 04.00 (January 2022) on (BOA_ADD_OFFSET); 0 is nodata
SCALE = 1e-4
OFFSET_FROM_BASELINE_4 = -0.1
# The largest reflectance the uint16 packing holds
PACKED_MAX = 65534 * SCALE + OFFSET_FROM_BASELINE_4


def offset(item: pystac.Item) -> float:
    """The reflectance offset of *item*'s processing baseline."""
    baseline = float(item.properties.get("s2:processing_baseline", "0"))
    return OFFSET_FROM_BASELINE_4 if baseline >= 4 else 0.0


def newest(items: list[pystac.Item]) -> list[pystac.Item]:
    """*items* with only the newest processing of each tile of each datatake.

    Planetary Computer lists reprocessed products (e.g. baselines 05.09 and
    05.10) next to the originals.
    """

    def tile(it):
        return pass_key(it), it.properties.get("s2:mgrs_tile")

    def rank(it):
        p = it.properties
        return float(p.get("s2:processing_baseline", "0")), p.get(
            "s2:generation_time", ""
        )

    best: dict[tuple, pystac.Item] = {}
    for it in sorted(items, key=rank):
        best[tile(it)] = it  # later (newer) processings replace earlier ones
    return [it for it in items if best[tile(it)] is it]


def pass_key(item: pystac.Item) -> str:
    """The pass (datatake) an item's tile belongs to, whatever its processing.

    ``s2:datatake_id`` ends in the processing baseline (``_N05.10``), so two
    processings of one pass would otherwise look like two passes.
    """
    datatake = item.properties.get("s2:datatake_id") or item.id
    return re.sub(r"_N\d{2}\.\d{2}$", "", datatake)


def source_attrs(items: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every item."""
    attrs = provenance.from_stac_collection(COLLECTION_URL)
    attrs["description"] = (
        "Sentinel-2 MSI Level-2A surface reflectance (Sen2Cor) and scene "
        "classification, as COGs on Microsoft Planetary Computer."
    )
    attrs["processing:lineage"] = (
        "tiles of each pass mosaicked, the newest processing of each, on a 10 m "
        "grid (20 and 60 m bands by nearest neighbour), clipped to the AOI; "
        "digital numbers scaled to reflectance by processing baseline, 0 set to NaN"
    )
    attrs["links"] += [
        provenance.link(
            "derived_from",
            f"{COLLECTION_URL}/items/{it.id}",
            type="application/geo+json",
        )
        for it in items
    ]
    return attrs


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    bands: tuple[Band, ...] = ("blue", "green", "red", "nir", "swir16", "scl"),
    max_cloud_cover: int = 100,
    windows: list[Window] | None = None,
) -> xr.Dataset:
    """Sentinel-2 L2A reflectance for every pass over *aoi*, lazily, at 10 m.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start (inclusive).
    end_date : str
        End (inclusive).
    bands : list of str
        Bands by common name (coastal, blue, green, red, rededge1-3, nir,
        nir08, nir09, swir16, swir22) and scl, the scene classification.
    max_cloud_cover : int
        Skip tiles with more cloud than this percentage (STAC eo:cloud_cover,
        over the whole 110 km tile, not the AOI); 100 keeps every tile.
    windows : list of Window
        Set by the runner from a stage file's windows: the passes each window
        selects are kept (see ``stac.select``). Default: every pass.

    Returns
    -------
    xarray.Dataset
        One variable per band along ``(time, y, x)``: float32 reflectance,
        stored as uint16 with a scale and offset, and ``scl`` as uint8 CF
        flags; with the coordinates ``platform``, ``relative_orbit``,
        ``cloud_cover``, ``mgrs_tiles`` and ``scene_ids`` along ``time``, and
        ``window`` when the windows are named.
    """
    if not 0 <= max_cloud_cover <= 100:
        raise ValueError(f"max_cloud_cover is a percentage, not {max_cloud_cover}")
    query = (
        {"eo:cloud_cover": {"lte": max_cloud_cover}} if max_cloud_cover < 100 else None
    )
    found = stac.search(COLLECTION, aoi, start_date, end_date, query)
    items = newest(found)
    if skipped := len(found) - len(items):
        log.info(
            "Sentinel-2: skipped %d older processing(s) of the same tiles", skipped
        )
    if not items:
        raise ValueError(
            f"no Sentinel-2 L2A scenes over the AOI from {start_date} to {end_date}"
            + (f" with at most {max_cloud_cover}% cloud" if query else "")
        )
    groups = stac.group(items, pass_key)
    labels: list[str | None] = [None] * len(groups)
    if windows:
        picked = stac.select(groups, windows)
        if not picked:
            raise ValueError("no Sentinel-2 L2A passes match the windows")
        labels = [name for name, _ in picked]
        groups = [list(scene) for _, scene in picked]
        items = [it for g in groups for it in g]
        for name, scene in picked:
            log.info("Sentinel-2: %s %s", name or "selected", scene[0].datetime)

    assets = [ASSETS[b] for b in bands]
    ds = stac.load_scenes(
        groups, assets, aoi, resolution=min(GSD.get(a, 20) for a in assets)
    )
    ds = ds.rename({ASSETS[b]: b for b in bands})
    log.info(
        "Sentinel-2: %d pass(es) of %d tile(s), %d x %d cells",
        len(groups),
        len(items),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    offsets = xr.DataArray([offset(g[0]) for g in groups], dims="time")
    for band in bands:
        if band == "scl":
            # odc-stac reads it as float32 (the items declare no data type);
            # outside the AOI and missing tiles are class 0, no_data
            ds[band] = ds[band].fillna(0).astype("uint8")
            ds[band].attrs = {
                "long_name": "Sentinel-2 L2A scene classification",
                "flag_values": list(SCL_CLASSES),
                "flag_meanings": " ".join(SCL_CLASSES.values()),
            }
            continue
        dn = ds[band]
        # Stored as uint16 (half of float32) with 0 for NaN, so the largest values
        # of pre-2022 baselines (no offset) clip at the top of that range
        reflectance = (dn.where(dn != 0) * SCALE + offsets).clip(max=PACKED_MAX)
        ds[band] = reflectance.astype("float32").assign_attrs(
            long_name=f"{band} surface reflectance ({ASSETS[band]})", units="1"
        )
        ds[band].encoding = {
            "dtype": "uint16",
            "scale_factor": SCALE,
            "add_offset": OFFSET_FROM_BASELINE_4,
            "_FillValue": 0,
        }
    first = [g[0].properties for g in groups]
    ds = ds.assign_coords(
        platform=("time", [p["platform"] for p in first]),
        relative_orbit=("time", np.array([p["sat:relative_orbit"] for p in first])),
        cloud_cover=(
            "time",
            np.array(
                [np.mean([it.properties["eo:cloud_cover"] for it in g]) for g in groups]
            ),
            {"long_name": "mean eo:cloud_cover of the pass's tiles", "units": "%"},
        ),
        mgrs_tiles=(
            "time",
            [" ".join(it.properties["s2:mgrs_tile"] for it in g) for g in groups],
        ),
        scene_ids=("time", [" ".join(it.id for it in g) for g in groups]),
    )
    if any(labels):
        ds = ds.assign_coords(window=("time", [name or "" for name in labels]))
    ds.attrs = {**source_attrs(items), **provenance.time_extent(ds)}
    return ds
