"""Waveforms from FDSN web services, stored as raw counts at their native rate.

The FDSN station service gives the channels in the AOI with their metadata
(``format=text``), and dataselect their miniSEED, which pymseed decodes. A
``waveforms/*`` group has these children (see dev/seismic.md):

- ``sites``: one row per position (``NET.STA.LOC``) along an xvec
  ``geometry`` dimension, like a ``stations/*`` group, with ``site`` ids,
  ``elevation``, ``distance_km`` and ``back_azimuth``;
- ``channels``: one row per channel requested (``NET.STA.LOC.CHA``), with its
  ``site``, orientation, ``sensitivity`` and the ``group`` holding its counts
  (empty for a channel that returned no data);
- one group per sample rate (``50hz``, ``100hz``, …): ``counts`` on
  ``(channel, time)``, with ``site`` and ``time_offset`` along ``channel``.

Each rate's time axis starts at the window's start, one step per sample. A
channel's samples go at the nearest step, and ``time_offset`` records the
shift: sample *i* was recorded at ``time[i] + time_offset``.
"""

from __future__ import annotations

import fnmatch
import logging
import time as clock
from collections import defaultdict
from typing import NamedTuple

import geopandas as gpd
import numpy as np
import pandas as pd
import pymseed
import pyproj
import requests
import shapely
import xarray as xr
from scipy import signal
from shapely.geometry import Point

log = logging.getLogger(__name__)

FDSN_URL = "https://service.earthscope.org/fdsnws"
# Gaps hold this value; xarray masks it to NaN on read
FILL = np.int32(np.iinfo(np.int32).min)
# Raw int32 counts allowed in memory for one load(), estimated before downloading
MAX_BYTES = 2 * 1024**3
STATIONS_PER_REQUEST = 50  # stations per dataselect POST
RETRIES = 3  # attempts per request, on 5xx errors and timeouts
TIMEOUT = 300  # seconds
CHUNK_SECONDS = 3600  # Zarr chunks: one channel by up to an hour of samples
SEED_ID = "SEED id NET.STA.LOC.CHA"
SITE_ID = "site id NET.STA.LOC (NET.STA for an empty location)"
# Text station columns (level=channel) -> channel table columns
COLUMNS = {
    "Network": "network",
    "Station": "station",
    "Location": "location",
    "Channel": "code",
    "Latitude": "latitude",
    "Longitude": "longitude",
    "Elevation": "elevation",
    "Depth": "depth",
    "Azimuth": "azimuth",
    "Dip": "dip",
    "SensorDescription": "sensor",
    "Scale": "sensitivity",
    "ScaleFreq": "sensitivity_frequency",
    "ScaleUnits": "sensitivity_units",
    "SampleRate": "sample_rate",
    "StartTime": "start",
    "EndTime": "end",
}
NUMERIC = (
    "latitude",
    "longitude",
    "elevation",
    "depth",
    "azimuth",
    "dip",
    "sensitivity",
    "sensitivity_frequency",
    "sample_rate",
)
# Variables of the `sites` group, from the channel table
SITE_VARS = {
    "network": {"long_name": "FDSN network code"},
    "station": {"long_name": "station code"},
    "location": {"long_name": "location code"},
    "elevation": {"long_name": "station elevation", "units": "m"},
}
# Variables of the `channels` group, from the channel table
CHANNEL_VARS = {
    "code": {"long_name": "channel code (band, instrument, orientation)"},
    "depth": {"long_name": "sensor depth below the station", "units": "m"},
    "azimuth": {"long_name": "sensor azimuth, clockwise from north", "units": "degree"},
    "dip": {"long_name": "sensor dip, down from horizontal", "units": "degree"},
    "sensitivity": {
        "long_name": "overall sensitivity: counts per sensitivity_units",
    },
    "sensitivity_frequency": {
        "long_name": "frequency of the sensitivity",
        "units": "Hz",
    },
    "sensitivity_units": {"long_name": "input units of the sensitivity"},
}


class Segment(NamedTuple):
    """Contiguous samples of one channel, as decoded."""

    rate: float  # Hz
    start: int  # first sample, ns since the epoch
    data: np.ndarray  # int32


def fdsn_time(t: pd.Timestamp) -> str:
    """*t* as an FDSN time, to the microsecond."""
    return t.strftime("%Y-%m-%dT%H:%M:%S.%f")


# ---- channels ----------------------------------------------------------------


def _units(units: str) -> str:
    # Services mix the case of the same units (m/s, M/S, M/S**2, Pa)
    u = units.strip().lower()
    return "Pa" if u == "pa" else u


def parse_channel_text(text: str) -> pd.DataFrame:
    """The channel table of an FDSN ``format=text&level=channel`` response.

    One row per channel epoch, with the columns of :data:`COLUMNS`. Empty
    location codes stay ``""``, empty numbers are NaN, an open ``end`` is
    NaT, and ``sensitivity_units`` are lowercased (``Pa`` aside).
    """
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return pd.DataFrame(columns=list(COLUMNS.values()))
    header = [c.strip().lstrip("#").strip() for c in lines[0].split("|")]
    rows = [[v.strip() for v in line.split("|")] for line in lines[1:]]
    table = pd.DataFrame(rows, columns=header).rename(columns=COLUMNS)
    for col in NUMERIC:
        table[col] = pd.to_numeric(table[col], errors="coerce")
    for col in ("start", "end"):
        table[col] = pd.to_datetime(table[col].replace("", None), format="ISO8601")
    table["sensitivity_units"] = table.sensitivity_units.map(_units)
    return table[list(COLUMNS.values())]


def seed_ids(table: pd.DataFrame) -> pd.Series:
    """``NET.STA.LOC.CHA`` for each row of a channel table."""
    return table.network + "." + table.station + "." + table.location + "." + table.code


def site_ids(table: pd.DataFrame) -> pd.Series:
    """``NET.STA.LOC`` for each row of a channel table, ``NET.STA`` for an
    empty location."""
    loc = ("." + table.location).where(table.location != "", "")
    return table.network + "." + table.station + loc


def query_channels(
    aoi: gpd.GeoDataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
    network: str = "*",
    station: str = "*",
    location: str = "*",
    channel: str = "*",
    url: str = FDSN_URL,
) -> pd.DataFrame:
    """Channel epochs active in ``[start, end]`` at stations inside *aoi*.

    The station service is queried with the AOI's WGS84 bounds, and stations
    outside the polygon are dropped. The other arguments are FDSN codes, with
    ``?`` and ``*`` wildcards and commas between alternatives; an empty
    location is ``--``.
    """
    aoi = aoi.to_crs("EPSG:4326")
    west, south, east, north = aoi.total_bounds
    params = {
        "network": network,
        "station": station,
        "location": location or "--",
        "channel": channel,
        "starttime": fdsn_time(start),
        "endtime": fdsn_time(end),
        "minlatitude": south,
        "maxlatitude": north,
        "minlongitude": west,
        "maxlongitude": east,
        "level": "channel",
        "format": "text",
        "includerestricted": "false",
        "nodata": 404,
    }
    r = requests.get(f"{url}/station/1/query", params=params, timeout=TIMEOUT)
    if r.status_code == 404:
        return parse_channel_text("")
    r.raise_for_status()
    table = parse_channel_text(r.text)
    inside = shapely.intersects_xy(
        aoi.geometry.union_all(), table.longitude, table.latitude
    )
    return table[inside].reset_index(drop=True)


def select_epochs(
    table: pd.DataFrame, start: pd.Timestamp, end: pd.Timestamp
) -> pd.DataFrame:
    """One epoch per channel: the one covering *start*, else the first
    overlapping ``[start, end]``.

    A sensitivity that changes inside the window is logged: the data keeps
    the one chosen.
    """
    active = table[
        (table.start <= end) & (table.end.isna() | (table.end > start))
    ].copy()
    active["covers"] = (active.start <= start) & (
        active.end.isna() | (active.end > start)
    )
    rows = []
    for sid, epochs in active.groupby(seed_ids(active), sort=False):
        if epochs.sensitivity.nunique() > 1:
            log.warning(
                "%s: sensitivity changes in the window (%s); using the window start's",
                sid,
                ", ".join(f"{s:g}" for s in epochs.sensitivity),
            )
        rows.append(epochs.sort_values(["covers", "start"], ascending=[False, True]))
    if not rows:
        return active.drop(columns="covers")
    return pd.concat([r.iloc[:1] for r in rows]).drop(columns="covers")


def prefer_bands(table: pd.DataFrame, channel: str) -> pd.DataFrame:
    """Per station, the channels of its preferred band and location.

    *channel* is the comma-separated FDSN list in order of preference
    (``"HH?,BH?"``). Each station keeps every channel of the first pattern
    it has, of one band code (the first two letters) and of its lowest
    location code, like peters-dome's ``_pick_channels``. So a station with
    both HH and BH channels keeps HHZ, HHN and HHE only, unlike FDSN's own
    union of every pattern.
    """
    patterns = [p.strip() for p in channel.split(",") if p.strip()]

    def rank(code: str) -> int:
        return next(
            (i for i, p in enumerate(patterns) if fnmatch.fnmatchcase(code, p)),
            len(patterns),
        )

    table = table.assign(_rank=table.code.map(rank), _band=table.code.str[:2])
    keep = []
    for _, chans in table.groupby(["network", "station"], sort=False):
        best = chans.sort_values(["_rank", "location", "_band"]).iloc[0]
        keep.append(
            chans[
                (chans._rank == best._rank)
                & (chans.location == best.location)
                & (chans._band == best._band)
            ]
        )
    if not keep:
        return table.drop(columns=["_rank", "_band"])
    return pd.concat(keep).drop(columns=["_rank", "_band"]).reset_index(drop=True)


def estimated_bytes(
    table: pd.DataFrame, start: pd.Timestamp, end: pd.Timestamp
) -> float:
    """In-memory size of the int32 counts of every channel in *table*."""
    return float(table.sample_rate.sum()) * (end - start).total_seconds() * 4


def check_size(
    table: pd.DataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
    limit: float = MAX_BYTES,
):
    """Raise ValueError if the channels' counts would exceed *limit* bytes."""
    size = estimated_bytes(table, start, end)
    if size > limit:
        raise ValueError(
            f"{len(table)} channels over {end - start} would be about "
            f"{size / 1e9:.1f} GB of counts, over the {limit / 1e9:.1f} GB limit; "
            "narrow the AOI, the window or the channels"
        )


# ---- waveforms -----------------------------------------------------------------


def _post(url: str, body: str) -> bytes:
    """POST *body* to *url*, retrying 5xx errors and timeouts with backoff.

    No data (204 or 404) is ``b""``; other errors raise.
    """
    error: Exception = RuntimeError("no attempts")
    for attempt in range(RETRIES):
        try:
            r = requests.post(url, data=body, timeout=TIMEOUT)
        except (requests.ConnectionError, requests.Timeout) as e:
            error = e
        else:
            if r.status_code in (204, 404):
                return b""
            if r.status_code < 500:
                r.raise_for_status()
                return r.content
            error = requests.HTTPError(f"{r.status_code} {r.reason}", response=r)
        if attempt + 1 < RETRIES:
            clock.sleep(2**attempt)
    raise error


def dataselect(
    table: pd.DataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
    url: str = FDSN_URL,
) -> list[bytes]:
    """miniSEED for every channel in *table* over ``[start, end]``.

    One POST per :data:`STATIONS_PER_REQUEST` stations. A request that
    still fails after its retries is logged and skipped, so one bad chunk
    doesn't lose the rest.
    """
    stations = list(dict.fromkeys(zip(table.network, table.station, strict=True)))
    t0, t1 = fdsn_time(start), fdsn_time(end)
    out = []
    for k in range(0, len(stations), STATIONS_PER_REQUEST):
        chunk = set(stations[k : k + STATIONS_PER_REQUEST])
        rows = table[[s in chunk for s in zip(table.network, table.station)]]
        body = "\n".join(
            f"{r.network} {r.station} {r.location or '--'} {r.code} {t0} {t1}"
            for r in rows.itertuples()
        )
        try:
            data = _post(f"{url}/dataselect/1/query", body + "\n")
        except requests.RequestException as e:
            log.warning(
                "dataselect request %d (%d stations) failed, skipped: %s",
                k // STATIONS_PER_REQUEST + 1,
                len(chunk),
                e,
            )
            continue
        if not data:
            log.info("dataselect request %d: no data", k // STATIONS_PER_REQUEST + 1)
        out.append(data)
    return out


def decode(buffers: list[bytes]) -> dict[str, list[Segment]]:
    """``{NET.STA.LOC.CHA: [Segment, …]}`` from miniSEED *buffers*.

    pymseed joins contiguous records into segments. Channels of non-integer
    samples (rare for seismic data) are logged and skipped, since counts are
    stored as int32.
    """
    segments: dict[str, list[Segment]] = defaultdict(list)
    for buf in buffers:
        if not buf:
            continue
        for tid in pymseed.MS3TraceList.from_buffer(buf, unpack_data=True):
            sid = ".".join(pymseed.sourceid2nslc(tid.sourceid))
            for seg in tid:
                if seg.sampletype != "i":
                    log.warning("%s: %s samples skipped", sid, seg.sampletype)
                    continue
                data = np.array(seg.datasamples, dtype=np.int32)
                segments[sid].append(Segment(seg.samprate, seg.starttime, data))
    return dict(segments)


def nominal_rate(rate: float) -> float:
    """*rate* without float jitter: a whole number of Hz if within 1e-4 of
    one, else rounded to 1e-6 Hz."""
    nearest = round(rate)
    if nearest and abs(rate - nearest) <= 1e-4 * rate:
        return float(nearest)
    return round(rate, 6)


def group_name(rate: float) -> str:
    """The group of a sample rate: ``50hz``, ``0.1hz``."""
    return f"{rate:g}hz"


def place(
    segments: list[Segment], start: int, n: int, rate: float, sid: str = ""
) -> tuple[np.ndarray, int]:
    """One channel's *segments* on the grid of *n* steps from *start* (ns).

    Each segment goes at the step nearest its first sample; steps without a
    sample hold :data:`FILL`, and where segments overlap the later one wins.
    Returns the samples and the channel's ``time_offset`` in ns: the shift of
    its longest segment from the grid, under half a step.
    """
    out = np.full(n, FILL, dtype=np.int32)
    filled = np.zeros(n, dtype=bool)
    longest = max(segments, key=lambda s: len(s.data))
    offset = 0
    overlaps = 0
    for seg in sorted(segments, key=lambda s: s.start):
        pos = (seg.start - start) * rate / 1e9
        i0 = int(np.floor(pos + 0.5))
        shift = round((pos - i0) / rate * 1e9)
        if seg is longest:
            offset = shift
        lo, hi = max(i0, 0), min(i0 + len(seg.data), n)
        if lo >= hi:
            continue
        overlaps += bool(filled[lo:hi].any())
        out[lo:hi] = seg.data[lo - i0 : hi - i0]
        filled[lo:hi] = True
    if overlaps:
        log.info("%s: %d overlapping segment(s); later data kept", sid, overlaps)
    return out, offset


def distances(
    lon: np.ndarray, lat: np.ndarray, origin: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """Distance (km) and back azimuth (degrees, station to *origin*, clockwise
    from north) of each point from *origin* ``(lon, lat)``, on WGS84."""
    lon0 = np.full(len(lon), origin[0])
    lat0 = np.full(len(lat), origin[1])
    _, back, dist = pyproj.Geod(ellps="WGS84").inv(lon0, lat0, lon, lat)
    return np.asarray(dist) / 1e3, np.asarray(back) % 360


def _sites(table: pd.DataFrame, origin: tuple[float, float]) -> xr.Dataset:
    """The ``sites`` group: one row per ``NET.STA.LOC`` of *table*, along an
    xvec ``geometry`` dimension (a site's first channel gives its position)."""
    table = table.assign(site=site_ids(table).to_numpy())
    spread = table.groupby("site")[["latitude", "longitude"]].agg(
        lambda s: s.max() - s.min()
    )
    for site in spread.index[(spread > 1e-4).any(axis=1)]:
        log.info("%s: channels at different positions; kept the first", site)
    rows = table.drop_duplicates("site").sort_values("site")
    dist, back = distances(rows.longitude.to_numpy(), rows.latitude.to_numpy(), origin)
    points = [Point(x, y) for x, y in zip(rows.longitude, rows.latitude, strict=True)]
    ds = xr.Dataset(
        {
            **{
                name: ("geometry", rows[name].to_numpy(), attrs)
                for name, attrs in SITE_VARS.items()
            },
            "distance_km": (
                "geometry",
                dist,
                {"long_name": "distance from the origin", "units": "km"},
            ),
            "back_azimuth": (
                "geometry",
                back,
                {"long_name": "azimuth from the site to the origin", "units": "degree"},
            ),
        },
        coords={
            "geometry": np.array(points, dtype=object),
            "site": ("geometry", rows.site.to_numpy(), {"long_name": SITE_ID}),
        },
    )
    return ds.xvec.set_geom_indexes("geometry", crs="EPSG:4326")


def _channels(table: pd.DataFrame, rate: pd.Series, group: pd.Series) -> xr.Dataset:
    """The ``channels`` group: every channel of *table* (indexed by SEED id),
    with the sample *rate* of its data and the *group* holding its counts."""
    return xr.Dataset(
        {
            "site": ("channel", site_ids(table).to_numpy(), {"long_name": SITE_ID}),
            **{
                name: ("channel", table[name].to_numpy(), attrs)
                for name, attrs in CHANNEL_VARS.items()
            },
            "sample_rate": (
                "channel",
                rate.to_numpy(),
                {
                    "long_name": "sample rate of the data, else of the metadata",
                    "units": "Hz",
                },
            ),
            "group": (
                "channel",
                group.to_numpy(),
                {"long_name": "group of the counts; empty if no data was returned"},
            ),
        },
        coords={"channel": ("channel", table.index.to_numpy(), {"long_name": SEED_ID})},
    )


def to_tree(
    table: pd.DataFrame,
    segments: dict[str, list[Segment]],
    start: pd.Timestamp,
    end: pd.Timestamp,
    origin: tuple[float, float],
) -> dict[str, xr.Dataset]:
    """The groups of a waveform dataset: ``sites``, ``channels``, and one per
    sample rate, keyed by :func:`group_name`.

    Each rate's group covers ``[start, end]`` (inclusive) with ``counts`` on
    ``(channel, time)``, and ``site`` and ``time_offset`` along ``channel``.
    ``sites`` and ``channels`` cover every channel in *table*, including
    those without data; distances are from *origin* ``(lon, lat)``.
    """
    table = table.set_index(seed_ids(table)).sort_index()
    rate = table.sample_rate.map(nominal_rate)
    group = pd.Series("", index=table.index, dtype=object)
    by_rate: dict[float, dict[str, list[Segment]]] = defaultdict(dict)
    for sid, segs in sorted(segments.items()):
        if sid not in table.index:
            log.debug("%s: data for a channel not requested, skipped", sid)
            continue
        # Trust the data's rate over the metadata's, keeping the main one
        rates: dict[float, int] = defaultdict(int)
        for s in segs:
            rates[nominal_rate(s.rate)] += len(s.data)
        main = max(rates, key=lambda r: rates[r])
        if len(rates) > 1:
            log.warning("%s: rates %s; kept %g Hz", sid, sorted(rates), main)
        if main != rate[sid]:
            log.info("%s: data at %g Hz, metadata %g Hz", sid, main, rate[sid])
        by_rate[main][sid] = [s for s in segs if nominal_rate(s.rate) == main]
        rate[sid], group[sid] = main, group_name(main)
    if missing := int((group == "").sum()):
        log.info("%d channel(s) returned no data", missing)

    t0 = start.value
    groups = {
        "sites": _sites(table, origin),
        "channels": _channels(table, rate, group),
    }
    for r in sorted(by_rate):
        chans = by_rate[r]
        # Every step from start up to and including end
        n = int(np.floor((end.value - t0) * r / 1e9 + 1e-9)) + 1
        steps = np.round(np.arange(n) * (1e9 / r)).astype("int64")
        counts = np.empty((len(chans), n), dtype=np.int32)
        offsets = np.zeros(len(chans), dtype="int64")
        for i, (sid, segs) in enumerate(chans.items()):
            counts[i], offsets[i] = place(segs, t0, n, r, sid)
        ds = xr.Dataset(
            {
                "counts": (
                    ("channel", "time"),
                    counts,
                    {"long_name": "raw digitizer counts", "units": "count"},
                )
            },
            coords={
                "channel": ("channel", list(chans), {"long_name": SEED_ID}),
                "time": ("time", (start + pd.to_timedelta(steps)).to_numpy()),
                "site": (
                    "channel",
                    site_ids(table.loc[list(chans)]).to_numpy(),
                    {"long_name": SITE_ID},
                ),
                "time_offset": (
                    "channel",
                    offsets.astype("timedelta64[ns]"),
                    {"long_name": "sample i was recorded at time[i] + time_offset"},
                ),
            },
        )
        ds.counts.encoding = {
            "_FillValue": FILL,
            "chunks": (1, min(n, round(CHUNK_SECONDS * r))),
        }
        groups[group_name(r)] = ds
    return groups


def to_dataset(ds: xr.Dataset) -> xr.Dataset:
    """A group of a waveform dataset, xvec CF-encoded for the store (``sites``;
    the others have no geometry), keeping ``attrs``."""
    encoded = ds.xvec.encode_cf()
    encoded.attrs = dict(ds.attrs)
    return encoded


# ---- derived views ---------------------------------------------------------------

# Physical quantity of each sensitivity unit
QUANTITIES = {
    "m/s": "velocity",
    "m/s**2": "acceleration",
    "m": "displacement",
    "Pa": "pressure",
}
COMPONENTS = ("Z", "N", "E", "1", "2")
TABLES = ("sites", "channels")


def _decoded(ds: xr.Dataset) -> xr.Dataset:
    return ds.xvec.decode_cf() if "geometry_container" in ds else ds


def rate_groups(node: xr.DataTree) -> list[str]:
    """The names of *node*'s sample-rate groups (``40hz``, ``50hz``, …), by rate."""
    names = [name for name in node.children if name not in TABLES]
    return sorted(names, key=lambda name: float(name.removesuffix("hz")))


def _counts(ds: xr.Dataset) -> xr.DataArray:
    """``counts`` as floats, gaps NaN, as read from the store or not."""
    da = ds.counts
    if da.dtype.kind in "iu":
        da = da.where(da != da.encoding.get("_FillValue", FILL))
    return da.astype("float64")


def to_physical(node: xr.DataTree, group: str) -> xr.DataArray:
    """The counts of *node*'s *group* / sensitivity, in each channel's
    ``sensitivity_units``.

    Named for the quantity (``velocity`` in m/s for seismometers,
    ``acceleration``, ``pressure`` in Pa) when every channel has the same
    units, else ``physical``. Dividing by the sensitivity alone assumes a
    flat response, true of broadband sensors well inside their passband
    (see dev/seismic.md).

    Parameters
    ----------
    node : xarray.DataTree
        A waveform dataset's node, e.g. ``tree["waveforms/seismic"]``, as read
        from the store or as returned by ``load()``.
    group : str
        One of its sample-rate groups, e.g. ``"50hz"``.
    """
    ds = node[group].to_dataset()
    meta = node["channels"].to_dataset().sel(channel=ds.channel)
    units = sorted(set(meta.sensitivity_units.values.tolist()))
    da = _counts(ds) / meta.sensitivity
    da.attrs = {"long_name": "counts / sensitivity"}
    if len(units) == 1:
        da.name = QUANTITIES.get(units[0], "physical")
        da.attrs["units"] = units[0]
    else:
        da.name = "physical"
        da.attrs["comment"] = "units differ by channel: see sensitivity_units"
    return da


def _resample(
    t: np.ndarray, x: np.ndarray, grid: np.ndarray, sos: np.ndarray | None
) -> np.ndarray:
    """*x* at times *t* (s), low-passed with *sos* per contiguous run of
    samples and linearly interpolated at *grid*; NaN outside every run."""
    out = np.full(len(grid), np.nan)
    good = ~np.isnan(x)
    edges = np.flatnonzero(np.diff(np.concatenate([[0], good.astype(int), [0]])))
    for a, b in zip(edges[::2].tolist(), edges[1::2].tolist(), strict=True):
        run = x[a:b]
        if sos is not None and b - a > 1:
            run = signal.sosfiltfilt(
                sos, run, padlen=min(3 * (2 * len(sos) + 1), b - a - 1)
            )
        inside = (grid >= t[a]) & (grid <= t[b - 1])
        out[inside] = np.interp(grid[inside], t[a:b], run)
    return out


def to_grid(node: xr.DataTree, rate: float = 20.0, physical: bool = True) -> xr.Dataset:
    """Every channel with data of a waveform dataset on one common time grid.

    Channels sampled faster than *rate* are first low-passed at 0.4 × *rate*
    (4th-order Butterworth, zero phase, per contiguous run of samples), then
    every channel is linearly interpolated from its true sample times
    (``time + time_offset``) onto a grid at *rate* spanning every group.
    Gaps, and times outside a channel's data, are NaN.

    Parameters
    ----------
    node : xarray.DataTree
        A waveform dataset's node, e.g. ``tree["waveforms/seismic"]``, as read
        from the store or as returned by ``load()``.
    rate : float
        Target sample rate in Hz. Default 20.
    physical : bool
        Grid counts / sensitivity (:func:`to_physical`) rather than counts.

    Returns
    -------
    xarray.Dataset
        One variable on ``(channel, time)``. The ``channels`` table's columns
        (``sample_rate`` the native one) and the ``sites`` table's, but not
        ``geometry``, are coordinates along ``channel``.
    """
    names = rate_groups(node)
    if not names:
        raise ValueError("no sample-rate groups found")
    meta = node["channels"].to_dataset()
    groups = [node[name].to_dataset() for name in names]
    start = min(g.time.values[0] for g in groups)
    end = max(g.time.values[-1] for g in groups)
    span = (end - start) / np.timedelta64(1, "ns") / 1e9
    n = int(np.floor(span * rate + 1e-9)) + 1
    grid = np.arange(n) / rate
    parts = []
    for name, g in zip(names, groups, strict=True):
        values = to_physical(node, name) if physical else _counts(g)
        native = float(meta.sample_rate.sel(channel=g.channel.values[0]))
        sos = (
            signal.butter(4, 0.4 * rate, fs=native, output="sos")
            if native > rate
            else None
        )
        out = np.empty((g.sizes["channel"], n))
        for i in range(g.sizes["channel"]):
            t = g.time.values + g.time_offset.values[i]
            seconds = (t - start) / np.timedelta64(1, "ns") / 1e9
            out[i] = _resample(seconds, values.values[i], grid, sos)
        parts.append(
            xr.Dataset(
                {values.name or "counts": (("channel", "time"), out, values.attrs)},
                coords={"channel": g.channel.values},
            )
        )
    ds = xr.concat(parts, dim="channel")
    # In the channels table's order (by SEED id), not by group
    ds = ds.sel(channel=[c for c in meta.channel.values if c in ds.indexes["channel"]])
    ds = ds.assign_coords(
        time=start + np.round(grid * 1e9).astype("int64").astype("timedelta64[ns]")
    )
    channels = meta.sel(channel=ds.channel)
    sites = _decoded(node["sites"].to_dataset())
    sites = sites.drop_vars("geometry").swap_dims(geometry="site")
    sites = sites.sel(site=channels.site)
    ds = ds.assign_coords(
        {k: channels[k] for k in channels.data_vars}
        | {k: sites[k] for k in sites.data_vars}
    )
    ds.attrs = {
        **node.attrs,
        "gaia:to_grid": f"{rate:g} Hz: anti-alias low-pass at {0.4 * rate:g} Hz "
        "(4th-order Butterworth, zero phase) where faster, then linear "
        "interpolation from time + time_offset",
    }
    return ds


def component(code: str, azimuth: float) -> str:
    """A channel's component: its last letter, with ``1``/``2`` read as ``N``/``E``
    when its azimuth is within 5° of 0°/90° (as peters-dome's ``_component``)."""
    c = code[-1]
    if c in "ZNE" or np.isnan(azimuth):
        return c
    az = azimuth % 180
    if min(az, 180 - az) <= 5:
        return "N"
    if abs(az - 90) <= 5:
        return "E"
    return c


def by_component(ds: xr.Dataset, sites: xr.Dataset | xr.DataTree) -> xr.Dataset:
    """A :func:`to_grid` result on ``(geometry, component, time)``.

    ``geometry`` holds the *sites* with data, in their order, with their variables (``site``,
    ``distance_km``, …) as coordinates, as a ``stations/*`` group's
    :func:`gaia_cli.stations.to_grid` does. ``component`` holds the present
    ones of Z, N, E, 1, 2 (see :func:`component`), and ``channel`` the SEED
    id of each (site, component), ``""`` where there is none.
    """
    if isinstance(sites, xr.DataTree):
        sites = sites.to_dataset()
    sites = _decoded(sites)
    (name,) = ds.data_vars
    comps = [
        component(c, a) for c, a in zip(ds.code.values, ds.azimuth.values, strict=True)
    ]
    present = set(ds.site.values.tolist())
    ids = [site for site in sites.site.values.tolist() if site in present]
    order = [c for c in COMPONENTS if c in comps]
    order += sorted(set(comps) - set(COMPONENTS))
    values = np.full((len(ids), len(order), ds.sizes["time"]), np.nan)
    channel = np.full((len(ids), len(order)), "", dtype=object)
    for i, (site, comp) in enumerate(zip(ds.site.values, comps, strict=True)):
        s, c = ids.index(site), order.index(comp)
        if channel[s, c]:
            log.info(
                "%s: %s already has a %s channel", ds.channel.values[i], site, comp
            )
            continue
        values[s, c] = ds[name].values[i]
        channel[s, c] = ds.channel.values[i]
    rows = sites.isel(geometry=[sites.site.values.tolist().index(s) for s in ids])
    out = rows.set_coords(list(rows.data_vars)).assign_coords(
        component=order,
        time=ds.time.values,
        channel=(("geometry", "component"), channel),
    )
    out[name] = (("geometry", "component", "time"), values, ds[name].attrs)
    out.attrs = dict(ds.attrs)
    return out
