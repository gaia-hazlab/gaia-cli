import os

import geopandas as gpd
import pandas as pd
import requests
import xarray as xr
from cyclopts import Parameter
from shapely.geometry import Point
from typing import Annotated
import xvec

from . import io

TOKEN = os.environ.get("SYNOPTIC_TOKEN")
SYNOPTIC_API = "https://api.synopticdata.com/v2"


def _check_token():
    if TOKEN is None:
        raise ValueError(
            "SYNOPTIC_TOKEN environment variable not set. "
            "Please set it to your Synoptic API token."
        )


def get_station_metadata(aoi: gpd.GeoDataFrame, bbox=False, buffer=0):
    """Retrieve station metadata from the Synoptic API within an area of interest.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        A GeoDataFrame containing a single polygon geometry that defines the
        area of interest (AOI). The function uses the AOI bounds to query the
        Synoptic API. CRS should be EPSG:4326 or will be interpreted as such.
    bbox : bool, optional
        If True, the API query will use the bounding box of ``aoi`` and the
        results will not be further clipped to the AOI polygon. Default is
        False (clip to polygon).
    buffer : float, optional
        If > 0 and ``bbox`` is False, the AOI geometry will be buffered by
        this amount (in the AOI CRS units) before testing station containment.

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
    The API requires a bounding-box style query; when ``bbox`` is False the
    function will additionally filter the returned stations to those whose
    points lie within the provided AOI polygon.
    """
    endpoint = "/stations/metadata"
    url = SYNOPTIC_API + endpoint

    precip_vars = ["precip_accum"]  # NOTE: there are a lot! e.g. hourly etc

    # API requires search by Bounding Box
    params = dict(  # state="wa",
        # NOTE: may be best to do point in polygon first to get list of stations to avoid larger bbox!
        # stid Single or comma separated list of SynopticLabs station IDs. Use a ! before any value to exclude matching stations
        bbox=",".join(map(str, aoi.total_bounds)),  # minx,miny,maxx,maxy
        token=TOKEN,
        status="active",
        sensorvars=True,  # return sensor variable info
        hfmetars=0,  # exclude 5min NOAA METARs sampling
        vars=",".join(precip_vars),  # restrict to precip only
        output="geojson",
    )
    response = requests.get(url, params=params)
    data = response.json()
    gfp = gpd.GeoDataFrame.from_features(data["features"], crs="EPSG:4326")

    # Reduce to stations only within AOI Polygon
    if not bbox:
        if buffer > 0:
            aoi.geometry = aoi.geometry.buffer(buffer)
        gfp = gfp[gfp.geometry.within(aoi.geometry.iloc[0])]

    print(f"Found {len(gfp)} stations with precip data in the area of interest.")

    return gfp


def get_station_timeseries(station_list: list, start_date: str, end_date: str):
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

    params = dict(
        token=TOKEN,
        stid=",".join(map(str, station_list)),
        start=start,
        end=end,
        # vars = 'precip' # 'Querying too many station hours. 195840.0 hours were requested.' ...
        vars="precip_accum",  # works, and gives hourly sampling
        precip=1,  # Enable derived precip (precip_accumulated_set_1d, precip_intervals_set_1d)
        # precip =0, (default) 'precip_accum_set_1' (Millimeters) cumulative sum of the interval values for the requested time period, in this case hourly timestamps
        obtimezone="UTC",  # ensure UTC returns
    )

    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()


def construct_uniform_datarray(data: dict):
    """Construct a uniform hourly xarray DataArray of accumulated precipitation.

    Parameters
    ----------
    data : dict
        The JSON-like dictionary returned by :func:`get_station_timeseries`.
        Expected to contain a top-level ``STATION`` list where each station
        contains an ``OBSERVATIONS`` mapping with a ``date_time`` index and
        precipitation fields such as ``precip_accumulated_set_1d``.

    Returns
    -------
    xarray.DataArray
        A 2D (time x station) DataArray of hourly accumulated precipitation
        values. Each station is assigned coordinates for ``station``,
        ``elevation``, and ``geometry`` (a shapely Point).

    Notes
    -----
    The function resamples station data to hourly timestamps using a
    backward-fill with a one-step limit (``bfill(limit=1)``) to avoid
    propagating values through long gaps. Stations with different time
    coverage are concatenated using ``join='outer'``.
    """

    data_arrays = []
    for ind in range(len(data["STATION"])):
        station = data["STATION"][ind]

        # Create DataFrame
        df_station = pd.DataFrame(station["OBSERVATIONS"])
        df_station["date_time"] = pd.to_datetime(
            df_station["date_time"]
        ).dt.tz_localize(None)  # Remove timezone
        df_station = df_station.set_index("date_time")
        df_station.index.name = "time"

        # Check if already hourly (median interval ~1 hour)
        time_diff = df_station.index.to_series().diff().dropna()
        median_interval = time_diff.median()

        if median_interval != pd.Timedelta(hours=1):
            print(
                f"Station {station['STID']} ({station['NAME']}): original interval = {median_interval}, resampling to hourly"
            )
            # Resample to hourly - use forward fill for accumulated precip (last known value)

        # I think this is more honest to the data, as it doesn't propagate values forward when there are multiple missing timestamps in a row, which could be the case if a station goes offline for a period of time.
        # 0:28 -> 0:00, but if multiple adjacent timesteps are missing, the rest get NaNs
        df_station = df_station.resample("h").bfill(limit=1)

        # Convert to DataArray
        varname = "precip_accumulated_set_1d"
        da_station = df_station[varname].to_xarray()

        # Assign coordinates
        # NOTE: Synoptic API always returns elevation in feet
        elevation_m = float(station["ELEVATION"]) * 0.3048
        da_station = da_station.assign_coords(
            station=station["STID"],
            elevation=elevation_m,
            geometry=Point(station["LONGITUDE"], station["LATITUDE"]),
        )
        da_station.name = "precipitation"

        data_arrays.append(da_station)

    # Concatenate all DataArrays along the station dimension
    # TODO: revisit keyword options here to make sure we are doing the right thing with coordinates and attributes
    return xr.concat(
        data_arrays, dim="station", coords="different", compat="equals", join="outer"
    )


# CLI entrypoint
def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])] = None,
    start_date: Annotated[str, Parameter(name=["--start", "-s"])] = None,
    end_date: Annotated[str, Parameter(name=["--end", "-e"])] = None,
    output_path: Annotated[str, Parameter(name=["--output", "-o"])] = None,
):
    """Stage in-situ precipitation from the Synoptic API.
    Precipitation timeseries is resampled to construct a uniform hourly DataArray.
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
        Path to write the resulting Zarr store.

    Examples
    --------
    Stage precipitation data for Skagit Basin December 2025 Atmospheric River:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage synoptic -i $AOI -s 2025-12-01 -e 2025-12-31 -o skagit_precip.zarr

    Returns
    -------
    None
    """
    _check_token()
    io.check_s3_env(output_path)

    aoi = io.load_aoi(vectorPath)
    DA = load(aoi, start_date, end_date)

    # Save as xvec-encoded Zarr (point geometries)
    io.save_xvec_zarr(DA, output_path)

    print("DONE!")


def load(aoi: gpd.GeoDataFrame, start_date: str, end_date: str) -> xr.DataArray:
    """Fetch and build an hourly precipitation DataArray for stations in *aoi*.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date/time (parseable by pandas.Timestamp).
    end_date : str
        End date/time (parseable by pandas.Timestamp).

    Returns
    -------
    xarray.DataArray
    """
    _check_token()
    stations_gdf = get_station_metadata(aoi)
    data = get_station_timeseries(stations_gdf.stid.to_list(), start_date, end_date)
    return construct_uniform_datarray(data)
