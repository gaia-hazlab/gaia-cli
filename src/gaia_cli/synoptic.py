import logging
import os
import re
from collections.abc import Sequence
from typing import Annotated

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import xarray as xr
from cyclopts import Parameter
from shapely.geometry import Point

from . import provenance, stations

log = logging.getLogger(__name__)

TOKEN = os.environ.get("SYNOPTIC_TOKEN")
SYNOPTIC_API = "https://api.synopticdata.com/v2"
RESET_MM = 10  # drops in raw gauge accumulation larger than this are treated as resets
PRECIP = "precip_accum"  # staged as "precipitation", a running total
# Written by hand (no STAC Collection): stations are mostly in the US, Alaska
# included, but the API serves networks worldwide
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [[None, None]]},
}
KIND = "stations"
TIMED = True
PROVENANCE = {
    "title": "Synoptic Data station observations",
    "description": "Every report of precipitation and other requested "
    "variables at weather stations, from the Synoptic Data API.",
    "license": "other",
    "gaia:license-note": "Synoptic redistributes data from many station "
    "networks with differing terms; check https://synopticdata.com and the "
    "networks' terms before sharing.",
    "providers": [
        {
            "name": "Synoptic Data PBC",
            "roles": ["host", "processor", "licensor"],
            "url": "https://synopticdata.com",
        }
    ],
    "version": "v2",  # API version
    "processing:lineage": "every report kept as a CF contiguous ragged array; "
    "precipitation is Synoptic's derived running total, or precip_accum "
    "reports de-reset when that is unavailable",
}


def source_attrs(ds: xr.Dataset) -> dict:
    """Provenance attributes for station data *ds* (from :func:`load`)."""
    # The API URL never includes the request parameters, which hold the token
    via = provenance.link("via", f"{SYNOPTIC_API}/stations/timeseries")
    return {**PROVENANCE, "links": [via], **provenance.time_extent(ds)}


def _check_token():
    if TOKEN is None:
        raise ValueError(
            "SYNOPTIC_TOKEN environment variable not set. "
            "Please set it to your Synoptic API token."
        )


def _check_response(data: dict):
    # Synoptic returns HTTP 200 for API errors, with details in SUMMARY
    summary = data.get("SUMMARY", {})
    if summary.get("RESPONSE_CODE") != 1:
        raise RuntimeError(
            f"Synoptic API error {summary.get('RESPONSE_CODE')}: "
            f"{summary.get('RESPONSE_MESSAGE')}"
        )


def accumulate_raw_precip(series: pd.Series) -> pd.Series:
    """Convert raw gauge accumulation to accumulation since the first timestamp.

    Fallback for when Synoptic derived precip (``precip=1``) is unavailable.
    Small decreases are treated as sensor jitter and large decreases
    (> ``RESET_MM``) as gauge resets. Precipitation falling between the last
    observation before a reset and the first one after it is not counted.
    """
    s = pd.to_numeric(series, errors="coerce").dropna()
    segment = (s.diff() < -RESET_MM).cumsum()
    s = s.groupby(segment).cummax()
    increments = s.diff().clip(lower=0).fillna(0)
    return increments.cumsum().reindex(series.index)


def get_station_metadata(aoi: gpd.GeoDataFrame, variables: Sequence[str] = (PRECIP,)):
    """Retrieve station metadata from the Synoptic API within an area of interest.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        A GeoDataFrame containing a single polygon geometry that defines the
        area of interest (AOI). The function uses the AOI bounds to query the
        Synoptic API.
    variables : sequence of str
        Synoptic variables; stations reporting any of them are returned.

    Returns
    -------
    geopandas.GeoDataFrame
        A GeoDataFrame of station metadata returned by the Synoptic API. The
        GeoDataFrame uses EPSG:4326 and contains point geometries and station
        attributes such as ``STID``, ``NAME``, ``LATITUDE``, ``LONGITUDE``,
        and ``ELEVATION`` when available.

    Raises
    ------
    requests.exceptions.RequestException
        If the HTTP request to the Synoptic API fails.

    Notes
    -----
    The API requires a bounding-box style query, so the returned stations are
    then filtered to those inside the AOI polygon.
    """
    aoi = aoi.to_crs("EPSG:4326")
    endpoint = "/stations/metadata"
    url = SYNOPTIC_API + endpoint

    # API requires search by Bounding Box
    params = {  # state="wa",
        # NOTE: may be best to do point in polygon first to get list of stations to avoid larger bbox!
        # stid Single or comma separated list of SynopticLabs station IDs. Use a ! before any value to exclude matching stations
        "bbox": ",".join(map(str, aoi.total_bounds)),  # minx,miny,maxx,maxy
        "token": TOKEN,
        "status": "active",
        "sensorvars": True,  # return sensor variable info
        "hfmetars": 0,  # exclude 5min NOAA METARs sampling
        "vars": ",".join(variables),  # stations reporting any of these
        "output": "geojson",
    }
    response = requests.get(url, params=params)
    data = response.json()
    if "features" not in data:  # geojson output only includes SUMMARY on error
        _check_response(data)
    gfp = gpd.GeoDataFrame.from_features(data["features"], crs="EPSG:4326")

    # Reduce to stations only within AOI Polygon
    gfp = gfp[gfp.geometry.within(aoi.geometry.iloc[0])]

    log.info("found %d stations with %s in the AOI", len(gfp), ", ".join(variables))

    return gfp


def get_station_timeseries(
    station_list: list,
    start_date: str,
    end_date: str,
    variables: Sequence[str] = (PRECIP,),
):
    """Fetch timeseries observations for a list of stations from Synoptic API.

    Parameters
    ----------
    station_list : list
        A list of Synoptic station identifiers (strings or ints).
    start_date : str
        Start timestamp for the query. Parseable by pandas.Timestamp (e.g.
        '2022-01-01 00:00').
    end_date : str
        End timestamp for the query. Parseable by pandas.Timestamp.
    variables : sequence of str
        Synoptic variables to fetch. With ``precip_accum``, Synoptic's derived
        running total is requested too.

    Returns
    -------
    dict
        The parsed JSON response from the Synoptic API as a Python dict. The
        dict has a top-level ``STATION`` key containing a list of station
        time series records.

    Raises
    ------
    requests.exceptions.RequestException
        If the HTTP request to the Synoptic API fails.
    """
    start = pd.Timestamp(start_date).strftime("%Y%m%d%H%M")
    end = pd.Timestamp(end_date).strftime("%Y%m%d%H%M")

    endpoint = "/stations/timeseries"
    url = SYNOPTIC_API + endpoint

    params = {
        "token": TOKEN,
        "stid": ",".join(map(str, station_list)),
        "start": start,
        "end": end,
        # vars = 'precip' # 'Querying too many station hours. 195840.0 hours were requested.' ...
        "vars": ",".join(variables),
        "obtimezone": "UTC",  # ensure UTC returns
    }
    if PRECIP in variables:
        # Enable derived precip (precip_accumulated_set_1d, precip_intervals_set_1d)
        # precip =0, (default) 'precip_accum_set_1' (Millimeters) cumulative sum of the interval values for the requested time period, in this case hourly timestamps
        params["precip"] = 1

    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()

    # Derived precip requires a paid Synoptic tier, fall back to raw accumulation
    if "precip" in params and data.get("SUMMARY", {}).get("RESPONSE_CODE") == 403:
        log.info("derived precip unavailable, computing from raw accumulation")
        params.pop("precip")
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()

    _check_response(data)
    return data


def _column(obs: pd.DataFrame, var: str) -> str | None:
    """The first sensor set reporting *var*, e.g. ``air_temp_set_1``."""
    sets = sorted(c for c in obs if re.fullmatch(rf"{re.escape(var)}_set_\d+d?", c))
    return sets[0] if sets else None


def _precip_source(station: dict) -> str:
    """The reported variable Synoptic derived the running total from."""
    sets = station.get("SENSOR_VARIABLES", {}).get("precip_accumulated", {})
    return next((v.get("derived_from", "") for v in sets.values()), PRECIP)


def construct_ragged(data: dict, variables: Sequence[str] = (PRECIP,)) -> xr.Dataset:
    """Every report in *data* as a CF contiguous ragged array.

    Parameters
    ----------
    data : dict
        The JSON-like dictionary returned by :func:`get_station_timeseries`.
    variables : sequence of str
        Synoptic variables to keep. ``precip_accum`` becomes
        ``precipitation``, the running total since the first report: Synoptic's
        derived ``precip_accumulated_set_1d`` when present, otherwise the raw
        accumulation de-reset with :func:`accumulate_raw_precip`. Other
        variables keep their name and take the first sensor set (``_set_1``),
        as numbers.

    Returns
    -------
    xarray.Dataset
        Station variables (``station``, ``name``, ``elevation`` in m,
        ``network_id``, ``precip_source``, ``row_size``) along ``geometry``,
        with an xvec index; ``time`` and each variable along ``obs``. Reports
        with none of the variables, and stations with no such reports, are
        dropped. See :mod:`gaia_cli.stations`.
    """
    units = data.get("UNITS", {})
    meta, reports = [], []
    for station in data["STATION"]:
        obs = pd.DataFrame(station["OBSERVATIONS"])
        obs.index = pd.DatetimeIndex(obs.pop("date_time")).tz_localize(None)
        obs.index.name = "time"
        table = pd.DataFrame(index=obs.index)
        for var in variables:
            if var == PRECIP:
                if "precip_accumulated_set_1d" in obs:
                    table["precipitation"] = obs["precip_accumulated_set_1d"]
                elif "precip_accum_set_1" in obs:
                    table["precipitation"] = accumulate_raw_precip(
                        obs["precip_accum_set_1"]
                    )
            elif (column := _column(obs, var)) is not None:
                table[var] = obs[column]
        table = table.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        if table.empty:
            log.debug("Station %s: no reports of %s", station["STID"], variables)
            continue
        reports.append(table.sort_index())
        # NOTE: Synoptic API always returns elevation in feet
        elevation = station.get("ELEVATION")
        meta.append(
            {
                "station": station["STID"],
                "name": station.get("NAME") or "",
                "elevation": float(elevation) * 0.3048 if elevation else np.nan,
                "network_id": station.get("MNET_ID") or "",
                "precip_source": _precip_source(station),
                "geometry": Point(
                    float(station["LONGITUDE"]), float(station["LATITUDE"])
                ),
            }
        )
    if not reports:
        raise ValueError(f"no station reported any of {', '.join(variables)}")

    stations_ds = xr.Dataset.from_dataframe(pd.DataFrame(meta).set_index("station"))
    stations_ds = stations_ds.set_coords("geometry")
    if PRECIP not in variables:
        stations_ds = stations_ds.drop_vars("precip_source")
    ds = stations.ragged(stations_ds, reports)
    ds.elevation.attrs = {"long_name": "station elevation", "units": "m"}
    if "precip_source" in ds:
        ds.precip_source.attrs["long_name"] = (
            "Synoptic variable the precipitation running total is derived from"
        )
    if "precipitation" in ds:
        ds.precipitation.attrs = {
            "long_name": "precipitation accumulated since the first report",
            "units": "mm",
        }
    for var in variables:
        if var in ds and var in units:
            ds[var].attrs["units"] = units[var]
    return ds


# CLI entrypoint
def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    start_date: Annotated[str, Parameter(name=["--start", "-s"])],
    end_date: Annotated[str, Parameter(name=["--end", "-e"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    variables: Annotated[
        tuple[str, ...], Parameter(name=["--variables", "-v"], consume_multiple=True)
    ] = (PRECIP,),
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
):
    """Stage weather station reports from the Synoptic API.
    Every report is kept, as a CF ragged array (see dev/stations.md).
    Requires SYNOPTIC_TOKEN environment variable (https://docs.synopticdata.com).
    Requires AWS_PROFILE environment variable if writing to S3.

    Parameters
    ----------
    vectorPath : str
        Path to a vector file (e.g. GeoJSON, Shapefile) that contains the
        area-of-interest polygon. The first geometry in the file is used.
    start_date : str
        Start date/time for the timeseries (parseable by pandas.Timestamp).
    end_date : str
        End date/time for the timeseries (parseable by pandas.Timestamp).
    output_path : str
        Output folder for the Zarr store (named after the folder), its STAC
        Item and the resolved stage.yaml.
    variables : list of str
        Synoptic variables, e.g. ``precip_accum air_temp``. ``precip_accum``
        is staged as ``precipitation``, a running total. Default precip_accum.
    dry_run : bool
        Validate the flags, read the AOI and print the resolved configuration
        without fetching any data.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.
    skip_existing : bool
        Add to an existing store, unless it already has this dataset's group.

    Examples
    --------
    Stage precipitation data for Skagit Basin December 2025 Atmospheric River:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage synoptic -i $AOI -s 2025-12-01 -e 2025-12-31 -o skagit-precip

    Returns
    -------
    None
    """
    from .stagefile import stage_dataset  # stagefile imports this module

    stage_dataset(
        "synoptic",
        vectorPath,
        start_date,
        end_date,
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
        variables=variables,
    )


def to_dataset(ds: xr.Dataset) -> xr.Dataset:
    """*ds* from :func:`load`, xvec CF-encoded for a DataTree group."""
    encoded = ds.xvec.encode_cf()
    encoded.attrs = dict(ds.attrs)
    return encoded


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: list[str] | tuple[str, ...] = (PRECIP,),
) -> xr.Dataset:
    """Every report of *variables* at stations in *aoi*, with provenance attributes.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date/time (parseable by pandas.Timestamp).
    end_date : str
        End date/time (parseable by pandas.Timestamp).
    variables : list of str
        Synoptic variables (see :func:`construct_ragged`).

    Returns
    -------
    xarray.Dataset
        A CF contiguous ragged array; :func:`gaia_cli.stations.to_grid` and
        :func:`gaia_cli.stations.to_dataframe` give other views of it.
    """
    _check_token()
    log.info("querying station metadata from %s", SYNOPTIC_API)
    stations_gdf = get_station_metadata(aoi, variables)
    log.info("fetching timeseries for %d station(s)", len(stations_gdf))
    data = get_station_timeseries(
        stations_gdf.stid.to_list(), start_date, end_date, variables
    )
    ds = construct_ragged(data, variables)
    ds.attrs = source_attrs(ds)
    return ds
