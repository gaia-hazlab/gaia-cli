"""Sentinel-1 RTC backscatter (gamma0) from Planetary Computer, clipped to an AOI.

Every pass over the AOI in the window is a time step of linear-power gamma0 on
the scenes' native 10 m UTM grid, with its platform and orbit as coordinates.
See dev/sentinel1-rtc.md.

Examples
--------
    $ gaia stage sentinel1-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-s1
    $ gaia stage sentinel1-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-s1 -r 131
"""

from __future__ import annotations

import logging
from enum import Enum
from typing import Literal

import geopandas as gpd
import numpy as np
import pystac
import xarray as xr

from .. import provenance, stac
from ..stac import Match
from ..window import Window

log = logging.getLogger(__name__)

COLLECTION = "sentinel-1-rtc"
COLLECTION_URL = f"{stac.PC_STAC}/collections/{COLLECTION}"
# Copied from COLLECTION_URL
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [["2014-10-10T00:28:21Z", None]]},
}
KIND = "raster"
# Takes pre/post windows with select and count (see stac.select)
SCENES = True
# `gaia align` default: the mean of linear power (average dB values after aligning)
REGRID_METHOD = "average"

Polarization = Literal["vv", "vh", "hh", "hv"]


class OrbitState(str, Enum):
    """Pass direction (STAC ``sat:orbit_state``); ``any`` keeps both."""

    any = "any"
    ascending = "ascending"
    descending = "descending"


def pass_key(item: pystac.Item) -> tuple[str, int]:
    """The pass an item's frame belongs to: its platform and absolute orbit."""
    return item.properties["platform"], int(item.properties["sat:absolute_orbit"])


def passes(items: list[pystac.Item]) -> list[list[pystac.Item]]:
    """*items* grouped into passes, the passes in time order."""
    return stac.group(items, pass_key)


def source_attrs(items: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every item."""
    attrs = provenance.from_stac_collection(COLLECTION_URL)
    # The collection's own description is ~3.7k characters of Markdown; the
    # derived_from link leads to it
    attrs["description"] = (
        "Sentinel-1 IW GRD backscatter, radiometrically terrain corrected (gamma0) "
        "by Catalyst, as 10 m COGs on Microsoft Planetary Computer."
    )
    attrs["processing:lineage"] = (
        "frames of each pass mosaicked, clipped to the AOI on the native grid; "
        "nodata and non-positive values set to NaN"
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
    """Sentinel-1 RTC gamma0 for every pass over *aoi*, lazily, on the native grid.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start (inclusive).
    end_date : str
        End (inclusive).
    polarizations : list of str
        Polarizations (vv, vh, hh, hv). Scenes without all of them are skipped.
    orbit_state : OrbitState
        Pass direction: any, ascending or descending.
    relative_orbit : int
        Relative orbit (track) number; 0 keeps every track.
    match : Match
        With pre and post windows in a stage file, ``relative_orbit`` chooses
        each window's scenes on each track separately and keeps the tracks
        that have scenes in both.
    windows : list of Window
        Set by the runner from a stage file's windows: the passes each window
        selects are kept (see ``stac.select``). Default: every pass.

    Returns
    -------
    xarray.Dataset
        One float32 variable per polarization along ``(time, y, x)``, with the
        coordinates ``platform``, ``orbit_state``, ``relative_orbit``,
        ``absolute_orbit`` and ``scene_ids`` along ``time``, and ``window``
        (e.g. ``pre``, ``post``) when the windows are named.
    """
    query = {}
    if orbit_state is not OrbitState.any:
        query["sat:orbit_state"] = {"eq": orbit_state.value}
    if relative_orbit:
        query["sat:relative_orbit"] = {"eq": relative_orbit}
    found = stac.search(COLLECTION, aoi, start_date, end_date, query or None)
    items = [it for it in found if set(polarizations) <= set(it.assets)]
    if skipped := len(found) - len(items):
        log.info(
            "Sentinel-1 RTC: skipped %d scene(s) without %s",
            skipped,
            "+".join(polarizations),
        )
    if not items:
        raise ValueError(
            f"no Sentinel-1 RTC scenes with {'+'.join(polarizations)} over the AOI "
            f"from {start_date} to {end_date}"
        )
    groups, labels = stac.choose(passes(items), windows, "Sentinel-1 RTC", match)
    items = [it for g in groups for it in g]
    ds = stac.load_scenes(groups, list(polarizations), aoi)
    log.info(
        "Sentinel-1 RTC: %d pass(es) of %d scene(s), %d x %d cells",
        len(groups),
        len(items),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    # -32768 is the assets' nodata; zero and negative gamma0 aren't valid either
    ds = ds.where(ds > 0)
    for pol in polarizations:
        ds[pol].attrs.pop("nodata", None)  # now NaN
        ds[pol].attrs.update(
            long_name=f"gamma0 {pol.upper()} backscatter, linear power", units="1"
        )
    first = [g[0].properties for g in groups]
    ds = ds.assign_coords(
        platform=("time", [p["platform"] for p in first]),
        orbit_state=("time", [p["sat:orbit_state"] for p in first]),
        relative_orbit=("time", np.array([p["sat:relative_orbit"] for p in first])),
        absolute_orbit=("time", np.array([p["sat:absolute_orbit"] for p in first])),
        scene_ids=("time", [" ".join(it.id for it in g) for g in groups]),
    )
    ds = stac.label_windows(ds, labels)
    ds.attrs = {**source_attrs(items), **provenance.time_extent(ds)}
    return ds
