"""NISAR L-band soil moisture (SME2, 200 m) frames from ASF, clipped to an AOI.

Every frame over the AOI in the window is a time step of JPL's Level-3 soil
moisture, on the global 200 m EASE-Grid 2.0, with its uncertainty and quality
flags. Frames of one pass overlap, and are retrieved separately: they are kept
apart, not mosaicked. PROVISIONAL products, the only ones released so far.
Read remotely from each frame's HDF5 file. Needs EARTHDATA_TOKEN. See
dev/nisar-sme2.md.

Examples
--------
    $ gaia stage nisar-sme2 -i odessa.geojson -s 2026-09-01 -e 2026-10-31 -o odessa-sme2
    $ gaia stage nisar-sme2 -i odessa.geojson -s 2026-09-01 -e 2026-10-31 -o odessa-sme2 --orbit-state ascending -v soilMoisture -v waterbodyFraction
"""

from __future__ import annotations

import datetime as dt
import logging
from typing import Any, Literal

import geopandas as gpd
import numpy as np
import pystac
import xarray as xr
from shapely.geometry import shape

from .. import earthdata, provenance, regrid, stac
from ..stac import Match
from ..window import Window
from .nisar_gcov import STAC, frame
from .sentinel1_rtc import OrbitState

log = logging.getLogger(__name__)

DATASET = "NISAR"
COLLECTION = "NISAR_L3_SME2_PROVISIONAL_V1"
# CMR's
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [["2025-10-28T00:00:00Z", None]]},
}
KIND = "raster"
# Takes pre/post windows with select and count (see stac.select)
SCENES = True
# `gaia align` default for the measurements; the flags and land cover use nearest
REGRID_METHOD = "average"
TOKENS = earthdata.TOKENS
# Read in worker processes: libhdf5 makes a process's reads one at a time. Each
# fetches 2 MiB at a time, the files' HDF5 page size: an open is 6 requests,
# not 19 (GDAL fixes this when a process first reads, hence processes of its own)
PROCESS_ENV = {"CPL_VSIL_CURL_CHUNK_SIZE": str(2 * 2**20)}

GRIDS = "//science/LSAR/SME2/grids"
# Grids of the files: the soil moisture retrieval's, then ancillaryData's
Variable = Literal[
    "soilMoisture",
    "soilMoistureUncertainty",
    "retrievalQualityFlag",
    "surfaceQualityFlag",
    "landCover",
    "localIncidenceAngle",
    "waterbodyFraction",
]
VARIABLES: tuple[Variable, ...] = (
    "soilMoisture",
    "soilMoistureUncertainty",
    "retrievalQualityFlag",
    "surfaceQualityFlag",
)
ANCILLARY = {"landCover", "localIncidenceAngle", "waterbodyFraction"}

FILL = -9999
FLOAT = {"data_type": "float32", "nodata": FILL}
BANDS: dict[str, dict] = {
    "soilMoisture": FLOAT,
    "soilMoistureUncertainty": FLOAT,
    "retrievalQualityFlag": {"data_type": "int16", "nodata": FILL},
    "surfaceQualityFlag": {"data_type": "int16", "nodata": FILL},
    "landCover": {"data_type": "int8", "nodata": -127},
    "localIncidenceAngle": FLOAT,
    "waterbodyFraction": FLOAT,
}
SOIL_MOISTURE = "volume_fraction_of_condensed_water_in_soil"
ATTRS = {
    "soilMoisture": {
        "long_name": "soil moisture: the mean of the recommended retrievals "
        "(DSG, PMI, TSR), else of all of them",
        "units": "m3 m-3",
        "standard_name": SOIL_MOISTURE,
    },
    "soilMoistureUncertainty": {
        "long_name": "soil moisture uncertainty",
        "units": "m3 m-3",
    },
    # Bits as the files' own description gives them (the 2024 product
    # specification's differ)
    "retrievalQualityFlag": {
        "long_name": "retrieval quality flag",
        "flag_masks": [1, 2, 4, 8],
        "flag_meanings": "retrieval_not_recommended pmi_applied tsr_applied "
        "dsg_applied",
    },
    # Bits undocumented in the files; see dev/nisar-sme2.md
    "surfaceQualityFlag": {"long_name": "surface quality flag"},
    "landCover": {"long_name": "land cover class (ancillary)"},
    "localIncidenceAngle": {
        "long_name": "local incidence angle (ancillary)",
        "units": "degree",
    },
    "waterbodyFraction": {
        "long_name": "water body fraction (ancillary)",
        "units": "1",
    },
}


def path(name: str) -> str:
    """*name*'s grid in the files."""
    return f"{GRIDS}/ancillaryData/{name}" if name in ANCILLARY else f"{GRIDS}/{name}"


def subdataset(href: str, name: str) -> str:
    """GDAL's name for the *name* grid in the HDF5 file at *href*."""
    return f'HDF5:"/vsicurl/{href}":{path(name)}'


def to_item(feature: dict) -> pystac.Item:
    """An asf_search SME2 *feature* as a pystac Item with an asset per grid,
    all at the HDF5 file's URL until :func:`load` points them at their
    subdatasets."""
    p = feature["properties"]
    item = pystac.Item(
        id=p["sceneName"],
        geometry=feature["geometry"],
        bbox=list(shape(feature["geometry"]).bounds),
        datetime=dt.datetime.fromisoformat(p["startTime"]),
        properties={
            "platform": p["platform"],
            "sat:orbit_state": p["flightDirection"].lower(),
            "sat:relative_orbit": int(p["pathNumber"]),
            "sat:absolute_orbit": int(p["orbit"]),
            "nisar:frame": int(p["frameNumber"]),
            "nisar:collection": p["collectionName"],
            "processing:datetime": p["processingDate"],
        },
        stac_extensions=[earthdata.RASTER_EXT],
    )
    for name, band in BANDS.items():
        item.add_asset(
            name,
            pystac.Asset(
                p["url"],
                media_type="application/x-hdf5",
                roles=["data"],
                extra_fields={"raster:bands": [band]},
            ),
        )
    return item


def grid(item: pystac.Item) -> tuple:
    """What an item's grid depends on: its frame (a frame keeps its window
    of the EASE-Grid on every date)."""
    p = item.properties
    return p["sat:relative_orbit"], p["nisar:frame"]


def source_attrs(items: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every frame."""
    url = f"{STAC}/{COLLECTION}_1"
    attrs = provenance.from_stac_collection(url)
    attrs["description"] = (
        "NISAR L-band Level-3 soil moisture (SME2, provisional), 200 m on "
        "EASE-Grid 2.0, from ASF DAAC."
    )
    attrs["version"] = COLLECTION
    attrs["processing:lineage"] = (
        "the newest processing of each frame; each frame its own time step, "
        "clipped to the AOI on the native 200 m EASE-Grid 2.0; -9999 → NaN in "
        "the floats"
    )
    attrs["links"] += [
        provenance.link(
            "derived_from", f"{url}/items/{it.id}", type="application/geo+json"
        )
        for it in items
    ]
    return attrs


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: tuple[Variable, ...] = VARIABLES,
    orbit_state: OrbitState = OrbitState.any,
    relative_orbit: int = 0,
    match: Match = Match.none,
    windows: list[Window] | None = None,
) -> xr.Dataset:
    """NISAR SME2 soil moisture for every frame over *aoi*, lazily, on the
    native 200 m EASE-Grid 2.0.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start (inclusive).
    end_date : str
        End (inclusive).
    variables : list of str
        Grids to read: soilMoisture, soilMoistureUncertainty and the bit flags
        retrievalQualityFlag and surfaceQualityFlag (default), and the
        ancillary landCover, localIncidenceAngle and waterbodyFraction.
    orbit_state : OrbitState
        Pass direction: any, ascending or descending.
    relative_orbit : int
        Relative orbit (track) number; 0 keeps every track.
    match : Match
        With pre and post windows in a stage file, ``relative_orbit`` chooses
        each window's frames on each track separately and keeps the tracks
        that have frames in both.
    windows : list of Window
        Set by the runner from a stage file's windows: the frames each window
        selects are kept (see ``stac.select``). Default: every frame.

    Returns
    -------
    xarray.Dataset
        One variable per grid along ``(time, y, x)``, one time step per frame
        (the frames of a pass are retrieved separately, and overlap): float32
        with -9999 as NaN, and integer flags and classes with their fill
        value. Coordinates ``platform``, ``orbit_state``, ``relative_orbit``,
        ``absolute_orbit`` and ``frame`` along ``time``, and ``window`` (e.g.
        ``pre``, ``post``) when the windows are named.
    """
    filters: dict[str, Any] = {"dataset": DATASET, "processingLevel": "SME2"}
    if orbit_state is not OrbitState.any:
        filters["flightDirection"] = OrbitState(orbit_state).value.upper()
    if relative_orbit:
        filters["relativeOrbit"] = relative_orbit
    features = earthdata.search(aoi, start_date, end_date, **filters)
    found = [
        to_item(f)
        for f in features
        if f["properties"].get("collectionName") == COLLECTION
    ]
    # Frames are reprocessed (CRIDs P05023, then P05036) into the same collection
    items = earthdata.newest(found, key=frame)
    if older := len(features) - len(items):
        log.info("NISAR SME2: skipped %d older or other processing(s)", older)
    if not items:
        raise ValueError(
            f"no NISAR SME2 frames over the AOI from {start_date} to {end_date}"
        )
    scenes, labels = stac.choose(stac.group(items, frame), windows, "NISAR SME2", match)
    items = [it for scene in scenes for it in scene]
    bands = list(variables)
    earthdata.sign(items, bands)
    for it in items:
        for name in bands:
            it.assets[name].href = subdataset(it.assets[name].href, name)
    # One header per frame: each date has the same grid
    earthdata.set_proj(items, bands[0], key=grid)
    ds = stac.load_scenes(scenes, bands, aoi)
    log.info(
        "NISAR SME2: %d frame(s), %d x %d cells of 200 m",
        len(items),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    for name in bands:
        da = ds[name]
        if np.issubdtype(da.dtype, np.floating):
            da = da.where(da != FILL)
        else:
            # The fill value, kept as the store's _FillValue
            da.encoding["_FillValue"] = BANDS[name]["nodata"]
        # Replaced, not merged: odc-stac's nodata attribute
        da.attrs = dict(ATTRS[name])
        if not np.issubdtype(da.dtype, np.floating):
            # TODO: teach `gaia align` CF flag_masks (bit flags), as it does flag_values
            da.attrs[regrid.METHOD_ATTR] = "nearest"
        ds[name] = da
    p = [it.properties for it in items]
    ds = ds.assign_coords(
        platform=("time", [x["platform"] for x in p]),
        orbit_state=("time", [x["sat:orbit_state"] for x in p]),
        relative_orbit=("time", np.array([x["sat:relative_orbit"] for x in p])),
        absolute_orbit=("time", np.array([x["sat:absolute_orbit"] for x in p])),
        frame=("time", np.array([x["nisar:frame"] for x in p])),
    )
    ds = stac.label_windows(ds, labels)
    ds.attrs = {**source_attrs(items), **provenance.time_extent(ds)}
    return ds
