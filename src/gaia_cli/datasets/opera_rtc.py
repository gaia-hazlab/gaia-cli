"""OPERA Sentinel-1 RTC backscatter (gamma0) bursts from ASF, clipped to an AOI.

Every pass over the AOI in the window is a time step of linear-power gamma0, the
pass's bursts mosaicked on their native 30 m UTM grid, with its platform and
orbit as coordinates. Each pass also gets its track's static layers: the
layover/shadow mask by default, and the incidence angles on request. Needs
EARTHDATA_TOKEN. See dev/opera-rtc.md.

Examples
--------
    $ gaia stage opera-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-opera
    $ gaia stage opera-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-opera -r 131
    $ gaia stage opera-rtc -i peters-dome.geojson -s 2023-09-01 -e 2023-09-30 -o peters-opera --static mask --static local_incidence_angle
"""

from __future__ import annotations

import logging
from typing import Any, Literal

import geopandas as gpd
import numpy as np
import pystac
import xarray as xr

from .. import earthdata, provenance, stac
from ..stac import Match
from ..window import Window
from .sentinel1_rtc import OrbitState, Polarization, pass_key

log = logging.getLogger(__name__)

DATASET = "OPERA-S1"
COLLECTION_URL = (
    "https://cmr.earthdata.nasa.gov/stac/ASF/collections/OPERA_L2_RTC-S1_V1_1"
)
STATIC_URL = (
    "https://cmr.earthdata.nasa.gov/stac/ASF/collections/OPERA_L2_RTC-S1-STATIC_V1_1"
)
# Copied from COLLECTION_URL
EXTENT = {
    "spatial": {"bbox": [[-180.0, -60.0, 180.0, 90.0]]},
    "temporal": {"interval": [["2016-04-14T00:00:00Z", None]]},
}
KIND = "raster"
# Takes pre/post windows with select and count (see stac.select)
SCENES = True
# `gaia align` default: the mean of linear power (and of the angles); the mask
# (flag_values) uses mode
REGRID_METHOD = "average"
TOKENS = earthdata.TOKENS

FLOAT = {"data_type": "float32", "nodata": "nan"}
BANDS = {pol: FLOAT for pol in ("vv", "vh", "hh", "hv")}

# Static layers: one granule per burst ID, the same on every date
StaticLayer = Literal[
    "mask",
    "local_incidence_angle",
    "incidence_angle",
    "number_of_looks",
    "rtc_anf_gamma0_to_beta0",
    "rtc_anf_gamma0_to_sigma0",
]
MASK = {0: "valid", 1: "shadow", 2: "layover", 3: "layover_and_shadow", 255: "no_data"}
STATIC_BANDS = {
    "mask": {"data_type": "uint8", "nodata": 255},
    "local_incidence_angle": FLOAT,
    "incidence_angle": FLOAT,
    "number_of_looks": FLOAT,
    "rtc_anf_gamma0_to_beta0": FLOAT,
    "rtc_anf_gamma0_to_sigma0": FLOAT,
}
STATIC_ATTRS = {
    "mask": {
        "long_name": "layover and shadow mask of the pass's track",
        "flag_values": list(MASK),
        "flag_meanings": " ".join(MASK.values()),
    },
    "local_incidence_angle": {
        "long_name": "local incidence angle, from the terrain normal",
        "units": "degree",
    },
    "incidence_angle": {
        "long_name": "incidence angle, from the ellipsoid normal",
        "units": "degree",
    },
    "number_of_looks": {"long_name": "number of looks", "units": "1"},
    "rtc_anf_gamma0_to_beta0": {
        "long_name": "area normalization factor, gamma0 to beta0",
        "units": "1",
    },
    "rtc_anf_gamma0_to_sigma0": {
        "long_name": "area normalization factor, gamma0 to sigma0",
        "units": "1",
    },
}


def burst(item: pystac.Item) -> str:
    return item.properties["opera:burst_id"]


def track(item: pystac.Item) -> int:
    return item.properties["sat:relative_orbit"]


def load_static(
    aoi: gpd.GeoDataFrame,
    tracks: list[int],
    layers: list[str],
    like: xr.Dataset,
) -> tuple[xr.Dataset, list[pystac.Item]]:
    """*layers* of each of *tracks* over *aoi*, along ``relative_orbit``, on
    *like*'s grid, and the static granules read."""
    features = earthdata.search(
        aoi, dataset=DATASET, processingLevel="RTC-STATIC", relativeOrbit=tracks
    )
    found = [earthdata.to_item(f, STATIC_BANDS) for f in features]
    items = [
        it for it in earthdata.newest(found, key=burst) if set(layers) <= set(it.assets)
    ]
    missing = set(tracks) - {track(it) for it in items}
    if missing:
        raise ValueError(
            f"no OPERA RTC static layers ({', '.join(layers)}) over the AOI for "
            f"relative orbit(s) {', '.join(map(str, sorted(missing)))}"
        )
    groups = sorted(stac.group(items, track), key=lambda g: track(g[0]))
    earthdata.sign(items, layers)
    earthdata.set_proj(items, layers[0])
    static = stac.load_scenes(groups, layers, aoi, geobox=like.odc.geobox)
    static = static.drop_vars("time").rename_dims(time="relative_orbit")
    static = static.assign_coords(relative_orbit=[track(g[0]) for g in groups])
    log.info(
        "OPERA RTC: static %s of relative orbit(s) %s, %d burst(s)",
        ", ".join(layers),
        ", ".join(map(str, static.relative_orbit.values)),
        len(items),
    )
    return static, items


def source_attrs(items: list[pystac.Item], static: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every burst and
    static granule."""
    attrs = provenance.from_stac_collection(COLLECTION_URL)
    attrs["description"] = (
        "OPERA Sentinel-1 IW burst backscatter, radiometrically terrain corrected "
        "(gamma0) by JPL, as 30 m COGs from ASF DAAC."
    )
    attrs["processing:lineage"] = (
        "the newest processing of each burst; bursts of each pass mosaicked, "
        "clipped to the AOI on the native grid; static layers mosaicked per "
        "track on the same grid and given to each pass of the track"
    )
    attrs["links"] += [
        provenance.link(
            "derived_from", f"{url}/items/{it.id}", type="application/geo+json"
        )
        for url, granules in ((COLLECTION_URL, items), (STATIC_URL, static))
        for it in granules
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
    static: tuple[StaticLayer, ...] = ("mask",),
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
    static : list of str
        Static layers of each pass's track (mask, local_incidence_angle,
        incidence_angle, number_of_looks, rtc_anf_gamma0_to_beta0,
        rtc_anf_gamma0_to_sigma0): the layover/shadow mask by default.
    windows : list of Window
        Set by the runner from a stage file's windows: the passes each window
        selects are kept (see ``stac.select``). Default: every pass.

    Returns
    -------
    xarray.Dataset
        One float32 variable per polarization and one per static layer (the
        uint8 ``mask`` as CF flags) along ``(time, y, x)``, with the
        coordinates ``platform``,
        ``orbit_state``, ``relative_orbit``, ``absolute_orbit`` and
        ``burst_ids`` along ``time``, and ``window`` (e.g. ``pre``, ``post``)
        when the windows are named.
    """
    filters: dict[str, Any] = {"dataset": DATASET, "processingLevel": "RTC"}
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
    groups, labels = stac.choose(
        stac.group(items, pass_key), windows, "OPERA RTC", match
    )
    items = [it for g in groups for it in g]
    earthdata.sign(items, polarizations)
    earthdata.set_proj(items, polarizations[0])
    ds = stac.load_scenes(groups, polarizations, aoi)
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
    tracks = [track(g[0]) for g in groups]
    static_items: list[pystac.Item] = []
    if static:
        layers, static_items = load_static(aoi, sorted(set(tracks)), list(static), ds)
        # Each pass gets its track's layers
        per_pass = layers.sel(relative_orbit=xr.DataArray(tracks, dims="time"))
        for layer in static:
            ds[layer] = per_pass[layer].drop_vars("relative_orbit")
            # Replaced, not merged: odc-stac's nodata would make the mask's 255
            # (a flag value) a _FillValue, read back as NaN
            ds[layer].attrs = STATIC_ATTRS[layer]
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
    ds = stac.label_windows(ds, labels)
    ds.attrs = {**source_attrs(items, static_items), **provenance.time_extent(ds)}
    return ds
