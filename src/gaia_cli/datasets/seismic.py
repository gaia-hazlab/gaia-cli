"""Seismic waveforms from EarthScope's FDSN services, as raw counts at native rate.

Channels at stations inside the AOI, one band per station (the first of
--channel it has, e.g. HH before BH). Each sample rate is its own group
(waveforms/seismic/50hz, ...), next to the sites and channels tables;
velocity is counts / sensitivity. See dev/seismic.md.

Examples
--------
Peters Dome, 2023-09-14: every HH or BH station within 250 km:
    $ gaia stage --config examples/peters-dome-seismic.yaml -o peters-dome-seismic

One station, two minutes:
    $ gaia stage seismic -i aoi.geojson -s 2023-09-14T02:06 -e 2023-09-14T02:08 -n AK --station BPAW -o bpaw
"""

import datetime as dt
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from .. import provenance, waveforms

log = logging.getLogger(__name__)

KIND = "waveforms"
# Read exactly the window: waveforms have no periods straddling its edges
PAD = pd.Timedelta(0)
# Written by hand (FDSN services have no STAC Collection): EarthScope serves
# networks worldwide, from the 1970s on
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [[None, None]]},
}
# Fetched with a second's margin, so samples near the window's edges have a
# grid step to snap to
MARGIN = dt.timedelta(seconds=1)
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


def _time(value: str) -> pd.Timestamp:
    t = pd.Timestamp(value)
    if isinstance(t, pd.Timestamp):  # not NaT
        return t
    raise ValueError(f"{value!r} is not a time")


def _origin(aoi: gpd.GeoDataFrame, origin: dict | None) -> dict:
    """The origin used for distances: the stage file's, else the AOI centroid."""
    if origin is not None:
        return {**origin, "source": "stage file origin"}
    centroid = aoi.to_crs(aoi.estimate_utm_crs()).centroid.to_crs(4326).iloc[0]
    return {
        "latitude": round(centroid.y, 6),
        "longitude": round(centroid.x, 6),
        "source": "AOI centroid",
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
    start, end = _time(start_date), _time(end_date)
    url = waveforms.FDSN_URL
    log.info("querying channels from %s", url)
    table = waveforms.query_channels(
        aoi, start, end, network, station, location, channel, url
    )
    table = waveforms.prefer_bands(waveforms.select_epochs(table, start, end), channel)
    if table.empty:
        raise ValueError(f"no {channel} channels at stations in the AOI")
    waveforms.check_size(table, start, end)
    log.info(
        "requesting %d channel(s) at %d station(s), about %.0f MB of counts",
        len(table),
        table.groupby(["network", "station"]).ngroups,
        waveforms.estimated_bytes(table, start, end) / 1e6,
    )
    buffers = waveforms.dataselect(table, start - MARGIN, end + MARGIN, url)
    segments = waveforms.decode(buffers)
    if not segments:
        raise RuntimeError("no waveform data returned")
    used = _origin(aoi, origin)
    groups = waveforms.to_tree(
        table, segments, start, end, (used["longitude"], used["latitude"])
    )
    query = {
        "url": url,
        "network": network,
        "station": station,
        "location": location,
        "channel": channel,
    }
    links = [
        provenance.link("via", f"{url}/station/1/", title="FDSN station service"),
        provenance.link("via", f"{url}/dataselect/1/", title="FDSN dataselect service"),
    ]
    rates = {k: g for k, g in groups.items() if k not in waveforms.TABLES}
    if not rates:
        raise RuntimeError("no waveform data returned for the channels requested")
    for name, ds in rates.items():
        log.info("%s: %d channel(s), %d samples each", name, *ds.counts.shape)
    times = np.concatenate([g.time.values[[0, -1]] for g in rates.values()])
    # The networks to cite: those with data
    networks = {s.split(".")[0] for g in rates.values() for s in g.site.values}
    tree = xr.DataTree.from_dict(groups)
    tree.attrs = {
        **PROVENANCE,
        "links": links,
        **provenance.time_extent(xr.Dataset(coords={"time": times})),
        "networks": sorted(networks),
        "gaia:origin": used,
        "gaia:query": query,
    }
    return tree


def to_dataset(ds: xr.Dataset) -> xr.Dataset:
    """A group from :func:`load`, xvec CF-encoded for the store."""
    return waveforms.to_dataset(ds)
