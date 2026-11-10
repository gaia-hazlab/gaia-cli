"""NISAR L-band GCOV backscatter (gamma0) frames from ASF, clipped to an AOI.

Every pass over the AOI in the window is a time step of linear-power gamma0, the
diagonal covariance terms (HHHH, HVHV) of frequency A, the pass's frames
mosaicked on their native UTM grid (10 m in 40 MHz mode), with the frames'
validity mask. VALIDATED products by default; PROVISIONAL ones (June-October
2026, before VALIDATED processing began) only on request, never mixed. Read
remotely from each frame's HDF5 file. Needs EARTHDATA_TOKEN. See
dev/nisar-gcov.md.

Examples
--------
    $ gaia stage nisar-gcov -i paradise.geojson -s 2026-10-21 -e 2026-11-08 -o paradise-nisar
    $ gaia stage nisar-gcov -i paradise.geojson -s 2026-07-01 -e 2026-10-20 -o paradise-nisar -m provisional -r 172
"""

from __future__ import annotations

import datetime as dt
import logging
from enum import Enum
from typing import Any, Literal

import geopandas as gpd
import numpy as np
import pystac
import xarray as xr
from shapely.geometry import shape

from .. import earthdata, provenance, stac
from ..stac import Match
from ..window import Window
from .sentinel1_rtc import OrbitState, pass_key

log = logging.getLogger(__name__)

DATASET = "NISAR"
STAC = "https://cmr.earthdata.nasa.gov/stac/ASF/collections"
# CMR's, for both maturities (the PROVISIONAL collection holds acquisitions
# reprocessed back to 2025-11 as well as June-October 2026)
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [["2025-10-28T00:00:00Z", None]]},
}
KIND = "raster"
# Takes pre/post windows with select and count (see stac.select)
SCENES = True
# `gaia align` default: the mean of linear power; the mask (flag_values) uses mode
REGRID_METHOD = "average"
TOKENS = earthdata.TOKENS
# Read in worker processes: libhdf5 makes a process's reads one at a time. Each
# fetches 4 MiB at a time, the files' HDF5 page size, so an open reads the
# metadata at the front of the file in a request or two (GDAL fixes this when a
# process first reads, hence processes of its own)
PROCESS_ENV = {"CPL_VSIL_CURL_CHUNK_SIZE": str(4 * 2**20)}
# The same quantities as sentinel1-rtc's and opera-rtc's vv, vh, hh, hv
COMMON_NAMES = {
    "HHHH": "gamma0_hh",
    "HVHV": "gamma0_hv",
    "VVVV": "gamma0_vv",
    "VHVH": "gamma0_vh",
}

# Diagonal covariance terms: real backscatter power (the off-diagonal ones are
# complex, and not read)
Term = Literal["HHHH", "HVHV", "VVVV", "VHVH"]
# Grids of frequency A besides the terms
Layer = Literal["mask", "numberOfLooks", "rtcGammaToSigmaFactor"]
GRIDS = "//science/LSAR/GCOV/grids/frequencyA"

FLOAT = {"data_type": "float32", "nodata": "nan"}
BANDS: dict[str, dict] = {
    **{t: FLOAT for t in COMMON_NAMES},
    "mask": {"data_type": "uint8", "nodata": 255},
    "numberOfLooks": FLOAT,
    "rtcGammaToSigmaFactor": FLOAT,
}
LAYER_ATTRS = {
    # 1 to the number of subswaths: a valid sample, from that subswath
    "mask": {
        "long_name": "valid samples subswath mask: 0 where any radar sample "
        "averaged was partly focused or invalid, else the subswath number",
        "flag_values": [0, 255],
        "flag_meanings": "invalid outside_swath",
    },
    "numberOfLooks": {"long_name": "number of radar looks", "units": "1"},
    "rtcGammaToSigmaFactor": {
        "long_name": "area normalization factor, gamma0 to sigma0",
        "units": "1",
    },
}


class Maturity(str, Enum):
    """Product maturity (asf_search ``dataMaturity``); each its own collection."""

    validated = "validated"
    provisional = "provisional"


COLLECTIONS = {
    Maturity.validated: "NISAR_L2_GCOV_V1",
    Maturity.provisional: "NISAR_L2_GCOV_PROVISIONAL_V1",
}


def subdataset(href: str, name: str) -> str:
    """GDAL's name for frequency A's *name* grid in the HDF5 file at *href*."""
    return f'HDF5:"/vsicurl/{href}":{GRIDS}/{name}'


def to_item(feature: dict) -> pystac.Item:
    """An asf_search GCOV *feature* as a pystac Item with an asset per grid of
    frequency A, all at the HDF5 file's URL until :func:`load` points them at
    their subdatasets."""
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
            "nisar:range_bandwidth": " ".join(p.get("rangeBandwidth") or []),
            "nisar:collection": p["collectionName"],
            "processing:datetime": p["processingDate"],
        },
        stac_extensions=[earthdata.RASTER_EXT],
    )
    terms = [f"{pol}{pol}" for pol in p.get("mainBandPolarization") or []]
    for name in [t for t in terms if t in COMMON_NAMES] + list(LAYER_ATTRS):
        item.add_asset(
            name,
            pystac.Asset(
                p["url"],
                media_type="application/x-hdf5",
                roles=["data"],
                extra_fields={"raster:bands": [BANDS[name]]},
            ),
        )
    return item


def frame(item: pystac.Item) -> tuple:
    """The frame and time an item was acquired, whatever its processing."""
    p = item.properties
    return p["sat:relative_orbit"], p["nisar:frame"], item.datetime


def grid(item: pystac.Item) -> tuple:
    """What an item's grid depends on: its frame and bandwidth mode."""
    p = item.properties
    return p["sat:relative_orbit"], p["nisar:frame"], p["nisar:range_bandwidth"]


def source_attrs(items: list[pystac.Item], maturity: Maturity) -> dict:
    """Provenance attributes: the collection's, plus a link to every frame."""
    url = f"{STAC}/{COLLECTIONS[maturity]}_1"
    attrs = provenance.from_stac_collection(url)
    attrs["description"] = (
        f"NISAR L-band geocoded polarimetric covariance (GCOV, {maturity.value}), "
        "radiometrically terrain corrected (gamma0), from ASF DAAC."
    )
    attrs["version"] = COLLECTIONS[maturity]
    attrs["processing:lineage"] = (
        "the newest processing of each frame; frequency A's grids of each pass's "
        "frames mosaicked, clipped to the AOI on the finest native grid; "
        "non-positive gamma0 set to NaN"
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
    terms: tuple[Term, ...] = ("HHHH", "HVHV"),
    maturity: Maturity = Maturity.validated,
    orbit_state: OrbitState = OrbitState.any,
    relative_orbit: int = 0,
    match: Match = Match.none,
    layers: tuple[Layer, ...] = ("mask",),
    windows: list[Window] | None = None,
) -> xr.Dataset:
    """NISAR GCOV gamma0 for every pass over *aoi*, lazily, on the native grid.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start (inclusive).
    end_date : str
        End (inclusive).
    terms : list of str
        Diagonal covariance terms of frequency A (HHHH, HVHV, VVVV, VHVH).
        Frames without all of them are skipped.
    maturity : Maturity
        Product maturity: validated, or provisional (June-October 2026,
        superseded by reprocessing in 2027). Never mixed.
    orbit_state : OrbitState
        Pass direction: any, ascending or descending.
    relative_orbit : int
        Relative orbit (track) number; 0 keeps every track.
    match : Match
        With pre and post windows in a stage file, ``relative_orbit`` chooses
        each window's passes on each track separately and keeps the tracks
        that have passes in both.
    layers : list of str
        Other grids of frequency A (mask, numberOfLooks,
        rtcGammaToSigmaFactor): the validity mask by default.
    windows : list of Window
        Set by the runner from a stage file's windows: the passes each window
        selects are kept (see ``stac.select``). Default: every pass.

    Returns
    -------
    xarray.Dataset
        One float32 variable per term and per float layer, and the uint8
        ``mask`` (CF flags), along ``(time, y, x)``, with the coordinates
        ``platform``, ``orbit_state``, ``relative_orbit``, ``absolute_orbit``,
        ``frames`` and ``range_bandwidth`` along ``time``, and ``window``
        (e.g. ``pre``, ``post``) when the windows are named.
    """
    maturity = Maturity(maturity)
    filters: dict[str, Any] = {
        "dataset": DATASET,
        "processingLevel": "GCOV",
        "dataMaturity": maturity.value.upper(),
    }
    if orbit_state is not OrbitState.any:
        filters["flightDirection"] = orbit_state.value.upper()
    if relative_orbit:
        filters["relativeOrbit"] = relative_orbit
    features = earthdata.search(aoi, start_date, end_date, **filters)
    # Only the one collection, whatever the search returns
    collection = COLLECTIONS[maturity]
    found = [
        to_item(f)
        for f in features
        if f["properties"].get("collectionName") == collection
    ]
    found = earthdata.newest(found, key=frame)
    if older := len(features) - len(found):
        log.info("NISAR GCOV: skipped %d older or other processing(s)", older)
    items = [it for it in found if set(terms) <= set(it.assets)]
    if skipped := len(found) - len(items):
        log.info("NISAR GCOV: skipped %d frame(s) without %s", skipped, "+".join(terms))
    if not items:
        raise ValueError(
            f"no {maturity.value} NISAR GCOV frames with {'+'.join(terms)} over the "
            f"AOI from {start_date} to {end_date}"
        )
    groups, labels = stac.choose(
        stac.group(items, pass_key), windows, "NISAR GCOV", match
    )
    items = [it for g in groups for it in g]
    bands = [*terms, *layers]
    earthdata.sign(items, bands)
    for it in items:
        for name in bands:
            it.assets[name].href = subdataset(it.assets[name].href, name)
    # One header per frame and mode: each date has the same grid (an HDF5 open
    # takes seconds, and libhdf5 makes them one at a time)
    earthdata.set_proj(items, terms[0], key=grid)
    finest = min(abs(it.properties["proj:transform"][0]) for it in items)
    ds = stac.load_scenes(groups, bands, aoi, resolution=finest)
    log.info(
        "NISAR GCOV: %d pass(es) of %d frame(s), %d x %d cells of %g m",
        len(groups),
        len(items),
        ds.sizes["y"],
        ds.sizes["x"],
        finest,
    )
    for term in terms:
        ds[term] = ds[term].where(ds[term] > 0)  # zero gamma0 isn't valid either
        ds[term].attrs = {
            "long_name": f"gamma0 {term[:2]} backscatter (covariance term {term}), "
            "linear power",
            "units": "1",
        }
    for layer in layers:
        # Replaced, not merged: odc-stac's nodata would make the mask's 255 (a
        # flag value) a _FillValue, read back as NaN
        ds[layer].attrs = LAYER_ATTRS[layer]
    first = [g[0].properties for g in groups]
    ds = ds.assign_coords(
        platform=("time", [p["platform"] for p in first]),
        orbit_state=("time", [p["sat:orbit_state"] for p in first]),
        relative_orbit=("time", np.array([p["sat:relative_orbit"] for p in first])),
        absolute_orbit=("time", np.array([p["sat:absolute_orbit"] for p in first])),
        frames=(
            "time",
            [" ".join(str(it.properties["nisar:frame"]) for it in g) for g in groups],
        ),
        range_bandwidth=("time", [p["nisar:range_bandwidth"] for p in first]),
    )
    ds = stac.label_windows(ds, labels)
    ds.attrs = {**source_attrs(items, maturity), **provenance.time_extent(ds)}
    return ds
