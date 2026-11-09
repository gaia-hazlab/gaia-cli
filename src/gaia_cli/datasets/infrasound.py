"""Infrasound waveforms from EarthScope's FDSN services, as raw counts at native rate.

Pressure channels inside the AOI, one per sensor position (the first of
--channel it has, e.g. BDF before HDF), so each element of an array is
kept, whether coded by station (IM.I53H1 to H8) or by location (AV.KENI
01 to 06). Laid out like `seismic` (waveforms/infrasound/40hz, ..., next to
the sites and channels tables); pressure in Pa is counts / sensitivity. See
dev/infrasound.md.

Examples
--------
Peters Dome, 2023-09-14: every BDF or HDF sensor within 400 km:
    $ gaia stage --config examples/peters-dome-waveforms.yaml -o peters-dome-waveforms

The IMS array I53US (Fairbanks), ten minutes:
    $ gaia stage infrasound -i aoi.geojson -s 2023-09-14T02:10 -e 2023-09-14T02:20 -n IM --station "I53H?" -o i53
"""

import geopandas as gpd
import pandas as pd
import xarray as xr

from .. import waveforms

KIND = "waveforms"
# Read exactly the window: waveforms have no periods straddling its edges
PAD = pd.Timedelta(0)
# Written by hand (FDSN services have no STAC Collection): EarthScope serves
# networks worldwide
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [[None, None]]},
}
PROVENANCE = {
    "title": "Infrasound waveforms from EarthScope FDSN web services",
    "description": "Raw digitizer counts of pressure channels (microbarometers "
    "and infrasound sensors) in the AOI, at each channel's native sample rate.",
    "license": "other",
    "gaia:license-note": "Waveforms come from many networks, each with its own "
    "terms and citation (most are open; see https://www.fdsn.org/networks/). "
    "Cite the networks in `network`.",
    "providers": [
        {
            "name": "EarthScope Consortium",
            "roles": ["host"],
            "url": "https://www.earthscope.org",
        }
    ],
    "processing": "raw counts, not demeaned, filtered or resampled; samples "
    "snapped to the nearest step of a grid starting at the window's start "
    "(the shift is time_offset); pressure in Pa = counts / sensitivity "
    "(gaia_cli.waveforms.to_physical), flat-response approximation",
}


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    network: str = "*",
    station: str = "*",
    location: str = "*",
    channel: str = "BDF,HDF",
    origin: dict | None = None,
) -> xr.DataTree:
    """Raw counts of the infrasound channels in *aoi*, one child per sample rate.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon; sensors inside it are kept.
    start_date : str
        Start time (inclusive), parseable by pandas.Timestamp.
    end_date : str
        End time (inclusive).
    network : str
        FDSN network codes, comma-separated, with ``?`` and ``*`` wildcards,
        e.g. ``AK,AV,IM``. Default: every network.
    station : str
        FDSN station codes, e.g. ``I53H?``. Default: every station.
    location : str
        FDSN location codes; ``--`` is the empty code. Default: every one;
        each sensor position keeps its lowest (an array's elements at their
        own positions are all kept).
    channel : str
        FDSN channel codes in order of preference: each sensor position keeps
        the first it has, e.g. BDF rather than HDF for ``BDF,HDF`` (the
        default).
    origin : dict, optional
        ``{latitude, longitude}`` for ``distance_km`` and ``back_azimuth``;
        the stage file's ``origin``, else the AOI centroid.

    Returns
    -------
    xarray.DataTree
        Children ``sites`` and ``channels`` (tables), and one per sample rate
        (``20hz``, ``40hz``, …) with ``counts`` on ``(channel, time)``; see
        :mod:`gaia_cli.waveforms`. Provenance attributes are on the tree.
    """
    return waveforms.load(
        aoi,
        start_date,
        end_date,
        network,
        station,
        location,
        channel,
        origin,
        PROVENANCE,
        arrays=True,
    )


def to_dataset(ds: xr.Dataset) -> xr.Dataset:
    """A group from :func:`load`, xvec CF-encoded for the store."""
    return waveforms.to_dataset(ds)
