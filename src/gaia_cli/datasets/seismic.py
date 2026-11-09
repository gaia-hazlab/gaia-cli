"""Seismic waveforms from EarthScope's FDSN services, as raw counts at native rate.

Channels at stations inside the AOI, one band per station (the first of
--channel it has, e.g. HH before BH). Each sample rate is its own group
(waveforms/seismic/50hz, ...), next to the sites and channels tables;
velocity is counts / sensitivity. See dev/seismic.md.

Examples
--------
Peters Dome, 2023-09-14: every HH or BH station within 250 km:
    $ gaia stage --config examples/peters-dome-waveforms.yaml -o peters-dome-waveforms

One station, two minutes:
    $ gaia stage seismic -i aoi.geojson -s 2023-09-14T02:06 -e 2023-09-14T02:08 -n AK --station BPAW -o bpaw
"""

import geopandas as gpd
import pandas as pd
import xarray as xr

from .. import waveforms

KIND = "waveforms"
# Read exactly the window: waveforms have no periods straddling its edges
PAD = pd.Timedelta(0)
# Written by hand (FDSN services have no STAC Collection): EarthScope serves
# networks worldwide, from the 1970s on
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [[None, None]]},
}
PROVENANCE = {
    "title": "Seismic waveforms from EarthScope FDSN web services",
    "description": "Raw digitizer counts of seismic channels at stations in "
    "the AOI, at each channel's native sample rate.",
    "license": "other",
    "gaia:license-note": "Waveforms come from many seismic networks, each "
    "with its own terms and citation (most are open; see "
    "https://www.fdsn.org/networks/). Cite the networks in `network`.",
    "providers": [
        {
            "name": "EarthScope Consortium",
            "roles": ["host"],
            "url": "https://www.earthscope.org",
        }
    ],
    "processing": "raw counts, not demeaned, filtered or resampled; samples "
    "snapped to the nearest step of a grid starting at the window's start "
    "(the shift is time_offset); physical units = counts / sensitivity "
    "(gaia_cli.waveforms.to_physical), flat-response approximation",
}


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    network: str = "*",
    station: str = "*",
    location: str = "*",
    channel: str = "HH?,BH?",
    origin: dict | None = None,
) -> xr.DataTree:
    """Raw counts of the seismic channels in *aoi*, one child per sample rate.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon; stations inside it are kept.
    start_date : str
        Start time (inclusive), parseable by pandas.Timestamp.
    end_date : str
        End time (inclusive).
    network : str
        FDSN network codes, comma-separated, with ``?`` and ``*`` wildcards,
        e.g. ``AK,AT``. Default: every network.
    station : str
        FDSN station codes, e.g. ``BPAW``. Default: every station.
    location : str
        FDSN location codes; ``--`` is the empty code. Default: every one;
        each station keeps its lowest.
    channel : str
        FDSN channel codes in order of preference: each station keeps every
        channel of the first it has, e.g. HHZ, HHN and HHE but not its BH
        channels for ``HH?,BH?`` (the default).
    origin : dict, optional
        ``{latitude, longitude}`` for ``distance_km`` and ``back_azimuth``;
        the stage file's ``origin``, else the AOI centroid.

    Returns
    -------
    xarray.DataTree
        Children ``sites`` and ``channels`` (tables), and one per sample rate
        (``50hz``, ``100hz``, …) with ``counts`` on ``(channel, time)``; see
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
    )


def to_dataset(ds: xr.Dataset) -> xr.Dataset:
    """A group from :func:`load`, xvec CF-encoded for the store."""
    return waveforms.to_dataset(ds)
