"""OPERA Sentinel-1 RTC backscatter (gamma0) bursts from ASF, clipped to an AOI.

Every pass over the AOI in the window is a time step of linear-power gamma0 and
its layover/shadow mask, the pass's bursts mosaicked on their native 30 m UTM
grid, with its platform and orbit as coordinates. Needs EARTHDATA_TOKEN. See
dev/opera-rtc.md.

Examples
--------
    $ gaia stage opera-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-opera
    $ gaia stage opera-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-opera -r 131
"""

from __future__ import annotations

import logging

import geopandas as gpd
import numpy as np
import pystac
import xarray as xr

from .. import earthdata, provenance, stac
from ..window import Window
from .sentinel1_rtc import Match, OrbitState, Polarization, pass_key

log = logging.getLogger(__name__)

DATASET = "OPERA-S1"
COLLECTION_URL = (
    "https://cmr.earthdata.nasa.gov/stac/ASF/collections/OPERA_L2_RTC-S1_V1_1"
)
# Copied from COLLECTION_URL
EXTENT = {
    "spatial": {"bbox": [[-180.0, -60.0, 180.0, 90.0]]},
    "temporal": {"interval": [["2016-04-14T00:00:00Z", None]]},
}
KIND = "raster"
# Takes pre/post windows with select and count (see stac.select)
SCENES = True
# `gaia align` default: the mean of linear power; the mask (flag_values) uses mode
REGRID_METHOD = "average"
TOKENS = earthdata.TOKENS

# The burst mask (the product guide: use it rather than the static one)
MASK = {0: "valid", 1: "shadow", 2: "layover", 3: "layover_and_shadow", 255: "no_data"}
BANDS = {
    **{
        pol: {"data_type": "float32", "nodata": "nan"}
        for pol in ("vv", "vh", "hh", "hv")
    },
    "mask": {"data_type": "uint8", "nodata": 255},
}


def source_attrs(items: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every burst."""
    attrs = provenance.from_stac_collection(COLLECTION_URL)
    attrs["description"] = (
        "OPERA Sentinel-1 IW burst backscatter, radiometrically terrain corrected "
        "(gamma0) by JPL, as 30 m COGs from ASF DAAC."
    )
    attrs["processing:lineage"] = (
        "the newest processing of each burst; bursts of each pass mosaicked, "
        "clipped to the AOI on the native grid"
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
    polarizations: tuple[Polarization, ...] = ("vv", "vh"),
    orbit_state: OrbitState = OrbitState.any,
    relative_orbit: int = 0,
    match: Match = Match.none,
    windows: list[Window] | None = None,
) -> xr.Dataset:
    """OPERA RTC-S1 gamma0 for every pass over *aoi*, lazily, on the native grid.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start (inclusive).
    end_date : str
        End (inclusive).
    polarizations : list of str
        Polarizations (vv, vh, hh, hv). Bursts without all of them are skipped.
    orbit_state : OrbitState
        Pass direction: any, ascending or descending.
    relative_orbit : int
        Relative orbit (track) number; 0 keeps every track.
    match : Match
        With pre and post windows in a stage file, ``relative_orbit`` chooses
        each window's passes on each track separately and keeps the tracks
        that have passes in both.
    windows : list of Window
        Set by the runner from a stage file's windows: the passes each window
        selects are kept (see ``stac.select``). Default: every pass.

    Returns
    -------
    xarray.Dataset
        One float32 variable per polarization and the uint8 ``mask`` (CF
        flags) along ``(time, y, x)``, with the coordinates ``platform``,
        ``orbit_state``, ``relative_orbit``, ``absolute_orbit`` and
        ``burst_ids`` along ``time``, and ``window`` (e.g. ``pre``, ``post``)
        when the windows are named.
    """
    filters = {"dataset": DATASET, "processingLevel": "RTC"}
    if orbit_state is not OrbitState.any:
        filters["flightDirection"] = orbit_state.value.upper()
    if relative_orbit:
        filters["relativeOrbit"] = relative_orbit
    features = earthdata.search(aoi, start_date, end_date, **filters)
    found = earthdata.newest([earthdata.to_item(f, BANDS) for f in features])
    if older := len(features) - len(found):
        log.info("OPERA RTC: skipped %d older processing(s) of the same bursts", older)
    items = [it for it in found if set(polarizations) <= set(it.assets)]
    if skipped := len(found) - len(items):
        log.info(
            "OPERA RTC: skipped %d burst(s) without %s",
            skipped,
            "+".join(polarizations),
        )
    if not items:
        raise ValueError(
            f"no OPERA RTC bursts with {'+'.join(polarizations)} over the AOI "
            f"from {start_date} to {end_date}"
        )
    groups = stac.group(items, pass_key)
    labels = [None] * len(groups)
    if windows:
        key = "sat:relative_orbit" if match is Match.relative_orbit else None
        picked = stac.select(groups, windows, key)
        if not picked:
            raise ValueError(
                "no OPERA RTC passes match the windows"
                + (" on any track with passes in every window" if key else "")
            )
        labels = [name for name, _ in picked]
        groups = [list(scene) for _, scene in picked]
        items = [it for g in groups for it in g]
        for name, scene in picked:
            log.info(
                "OPERA RTC: %s %s, relative orbit %s",
                name or "selected",
                scene[0].datetime,
                scene[0].properties["sat:relative_orbit"],
            )
    bands = [*polarizations, "mask"]
    earthdata.sign(items, bands)
    earthdata.set_proj(items, polarizations[0])
    ds = stac.load_scenes(groups, bands, aoi)
    log.info(
        "OPERA RTC: %d pass(es) of %d burst(s), %d x %d cells",
        len(groups),
        len(items),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    for pol in polarizations:
        ds[pol] = ds[pol].where(ds[pol] > 0)  # zero gamma0 isn't valid either
        ds[pol].attrs = {
            "long_name": f"gamma0 {pol.upper()} backscatter, linear power",
            "units": "1",
        }
    ds["mask"].attrs = {
        "long_name": "layover and shadow mask",
        "flag_values": list(MASK),
        "flag_meanings": " ".join(MASK.values()),
    }
    first = [g[0].properties for g in groups]
    ds = ds.assign_coords(
        platform=("time", [p["platform"] for p in first]),
        orbit_state=("time", [p["sat:orbit_state"] for p in first]),
        relative_orbit=("time", np.array([p["sat:relative_orbit"] for p in first])),
        absolute_orbit=("time", np.array([p["sat:absolute_orbit"] for p in first])),
        burst_ids=(
            "time",
            [" ".join(it.properties["opera:burst_id"] for it in g) for g in groups],
        ),
    )
    if any(labels):
        ds = ds.assign_coords(window=("time", [name or "" for name in labels]))
    ds.attrs = {**source_attrs(items), **provenance.time_extent(ds)}
    return ds
