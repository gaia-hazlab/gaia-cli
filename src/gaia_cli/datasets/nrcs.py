"""Soil moisture and other hourly reports at NRCS SNOTEL and SCAN stations.

From the NRCS Air and Water Database (AWDB) REST API: SNOTEL, SNOTEL-lite, SCAN and
Tribal SCAN stations, in the western US, Alaska and (SCAN) across the US. No token
needed. Every hourly report is kept, as a CF ragged array (see dev/stations.md), in
SI units and UTC. Sensors below ground are named for their depth (``SMS_5cm``).

Examples
--------
    $ gaia stage nrcs -i wa.geojson -s 2026-10-01 -e 2026-10-31 -o wa-soil
    $ gaia stage nrcs -i wa.geojson -s 2026-01-01 -e 2026-03-31 -v SMS STO WTEQ -o wa-snow
"""

import logging
from collections.abc import Sequence

import geopandas as gpd
import pandas as pd
import xarray as xr

from .. import fetch, provenance, stations

log = logging.getLogger(__name__)

AWDB_API = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"
# SNOTEL, SNOTEL-lite, SCAN and Tribal SCAN: the networks with soil sensors
NETWORKS = ("SNTL", "SNTLT", "SCAN", "CSCAN")
SOIL_MOISTURE = "SMS"  # AWDB element code: volumetric soil moisture, percent
CHUNK = 10  # stations per data request
# AWDB's stored units to SI: value * scale + offset
SI = {
    "in": (25.4, 0.0, "mm"),
    "degF": (5 / 9, -160 / 9, "degC"),
    "mph": (0.44704, 0.0, "m s-1"),
}
# Soil moisture elements, stored in percent, are staged as a volume fraction
VOLUMETRIC = {"SMS", "SMV", "SMN", "SMX"}
DESCRIPTIONS = {
    "SMS": ("volume_fraction_of_condensed_water_in_soil", "volumetric soil moisture"),
    "STO": ("soil_temperature", "soil temperature"),
    "TOBS": ("air_temperature", "air temperature"),
    "WTEQ": ("lwe_thickness_of_surface_snow_amount", "snow water equivalent"),
    "SNWD": ("surface_snow_thickness", "snow depth"),
    "PREC": (
        "lwe_thickness_of_precipitation_amount",
        "precipitation accumulated since the start of the water year (October 1)",
    ),
}
# Written by hand from AWDB's station list (2026-11-09): SCAN reaches Puerto Rico,
# the Virgin Islands and Hawaii, SNOTEL northern Alaska. Hourly soil moisture
# starts in the 1990s at the oldest stations
EXTENT = {
    "spatial": {"bbox": [[-165.2, 17.7, -64.7, 70.3]]},
    "temporal": {"interval": [[None, None]]},
}
KIND = "stations"
PROVENANCE = {
    "title": "NRCS SNOTEL and SCAN station observations",
    "description": "Every hourly report of soil moisture and other requested "
    "elements at NRCS SNOTEL, SNOTEL-lite, SCAN and Tribal SCAN stations, from the "
    "Air and Water Database (AWDB).",
    "license": "other",  # U.S. public domain has no SPDX identifier
    "gaia:license-note": "U.S. Government work, in the public domain.",
    "providers": [
        {
            "name": "USDA Natural Resources Conservation Service, National Water "
            "and Climate Center",
            "roles": ["producer", "host", "licensor"],
            "url": "https://wcc.sc.egov.usda.gov/awdbRestApi/swagger-ui/index.html",
        }
    ],
    "version": "v1",  # API version
    "processing:lineage": "every hourly report kept as a CF contiguous ragged "
    "array; data AWDB marks suspect left out; times converted from each station's "
    "local standard time to UTC; inches to mm, degF to degC, mph to m s-1 and "
    "soil moisture percent to m3 m-3",
}


def source_attrs(ds: xr.Dataset) -> dict:
    """Provenance attributes for station data *ds* (from :func:`load`)."""
    links = [
        provenance.link("via", f"{AWDB_API}/data"),
        provenance.link(
            "describedby",
            "https://wcc.sc.egov.usda.gov/awdbRestApi/swagger-ui/index.html",
            type="text/html",
        ),
    ]
    return {**PROVENANCE, "links": links, **provenance.time_extent(ds)}


def get_stations(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: Sequence[str] = (SOIL_MOISTURE,),
) -> gpd.GeoDataFrame:
    """Stations in *aoi* with any of *variables*, open at some time in the window.

    Closed stations are included: AWDB keeps their records.
    """
    params = {
        "stationTriplets": ",".join(f"*:*:{n}" for n in NETWORKS),
        "elements": ",".join(f"{v}:*" for v in variables),
        "activeOnly": "false",
    }
    response = fetch.get(f"{AWDB_API}/stations", params=params)
    response.raise_for_status()
    table = pd.DataFrame(response.json())
    gdf = gpd.GeoDataFrame(
        table,
        geometry=gpd.points_from_xy(table.longitude, table.latitude),
        crs="EPSG:4326",
    )
    gdf = gdf[gdf.within(aoi.to_crs("EPSG:4326").union_all())]
    start, end = pd.Timestamp(start_date), pd.Timestamp(end_date)
    gdf = gdf[
        (pd.to_datetime(gdf.beginDate) <= end) & (pd.to_datetime(gdf.endDate) >= start)
    ]
    log.info("found %d stations with %s in the AOI", len(gdf), ", ".join(variables))
    return gdf


def get_data(
    triplets: Sequence[str],
    start_date: str,
    end_date: str,
    variables: Sequence[str] = (SOIL_MOISTURE,),
) -> list[dict]:
    """Hourly AWDB data for *triplets*, a list with one entry per station.

    AWDB's dates are each station's local standard time, UTC-4 to UTC-11 in the
    US, so the request starts 12 hours early; :func:`construct_ragged` converts
    the times to UTC and trims them to the window.
    """
    begin = pd.Timestamp(start_date) - pd.Timedelta(hours=12)
    params = {
        "elements": ",".join(f"{v}:*" for v in variables),
        "duration": "HOURLY",
        "beginDate": begin.strftime("%Y-%m-%d %H:%M"),
        "endDate": pd.Timestamp(end_date).strftime("%Y-%m-%d %H:%M"),
    }
    data = []
    for i in range(0, len(triplets), CHUNK):
        chunk = ",".join(triplets[i : i + CHUNK])
        response = fetch.get(
            f"{AWDB_API}/data", params={**params, "stationTriplets": chunk}
        )
        response.raise_for_status()
        data.extend(response.json())
    return data


def _name(element: dict) -> str:
    """The staged name of an AWDB station element: ``SMS_5cm`` for SMS at -2 in."""
    name = element["elementCode"]
    if (h := element.get("heightDepth")) is not None:
        name += f"_{round(abs(h) * 2.54)}cm" + ("_above" if h > 0 else "")
    if element.get("ordinal", 1) > 1:
        name += f"_{element['ordinal']}"
    return name


def _attrs(element: dict) -> tuple[float, float, dict]:
    """Scale and offset to SI, and the staged variable's attributes."""
    code, unit = element["elementCode"], element.get("storedUnitCode", "")
    if code in VOLUMETRIC:
        scale, offset, unit = 0.01, 0.0, "m3 m-3"
    else:
        scale, offset, unit = SI.get(unit, (1.0, 0.0, unit))
    standard_name, long_name = DESCRIPTIONS.get(code, ("", code))
    attrs = {"standard_name": standard_name} if standard_name else {}
    attrs |= {"long_name": long_name, "units": unit, "awdb_element": code}
    if (h := element.get("heightDepth")) is not None:
        key = "height" if h > 0 else "depth"
        attrs[key] = f"{abs(h) * 2.54:.4g} cm"
        attrs["long_name"] += f" at {abs(h) * 2.54:.4g} cm {key}"
    return scale, offset, attrs


def construct_ragged(
    sites: gpd.GeoDataFrame, data: list[dict], start_date: str, end_date: str
) -> xr.Dataset:
    """Every report in *data* within the window, as a CF contiguous ragged array.

    Parameters
    ----------
    sites : geopandas.GeoDataFrame
        Stations from :func:`get_stations`.
    data : list of dict
        AWDB data from :func:`get_data`.
    start_date, end_date : str
        The window, inclusive, in UTC.

    Returns
    -------
    xarray.Dataset
        Station variables (``station``, the AWDB triplet; ``name``, ``elevation``
        in m, ``network_id``, ``utc_offset``, ``row_size``) along ``geometry``,
        with an xvec index; ``time`` (UTC) and one variable per station element
        along ``obs``, named by :func:`_name`. Reports with no values, and stations
        with no such reports, are dropped.
    """
    start, end = pd.Timestamp(start_date), pd.Timestamp(end_date)
    by_triplet = sites.set_index("stationTriplet")
    attrs, meta, reports = {}, [], []
    for record in data:
        triplet = record["stationTriplet"]
        site = by_triplet.loc[triplet]
        # AWDB dates are local standard time all year, a fixed offset from UTC
        offset = pd.Timedelta(hours=float(site.dataTimeZone))
        columns = {}
        for series in record.get("data", []):
            element = series["stationElement"]
            name = _name(element)
            scale, add, attrs[name] = _attrs(element)
            values = {
                pd.Timestamp(v["date"]) - offset: v["value"] * scale + add
                for v in series.get("values", [])
                if v.get("value") is not None
            }
            columns[name] = pd.Series(values, dtype="float64")
        table = pd.DataFrame(columns)
        table.index = pd.DatetimeIndex(table.index, name="time")
        table = table[(table.index >= start) & (table.index <= end)]
        table = table.dropna(how="all").sort_index()
        if table.empty:
            log.debug("Station %s: no reports in the window", triplet)
            continue
        reports.append(table)
        meta.append(
            {
                "station": triplet,
                "name": site["name"],
                "elevation": float(site.elevation) * 0.3048,  # AWDB: feet
                "network_id": site.networkCode,
                "utc_offset": float(site.dataTimeZone),
                "geometry": site.geometry,
            }
        )
    if not reports:
        raise ValueError("no NRCS station reported in the window")

    stations_ds = xr.Dataset.from_dataframe(pd.DataFrame(meta).set_index("station"))
    ds = stations.ragged(stations_ds.set_coords("geometry"), reports)
    ds.elevation.attrs = {"long_name": "station elevation", "units": "m"}
    ds.network_id.attrs["long_name"] = "AWDB network code (SNTL, SNTLT, SCAN, CSCAN)"
    ds.utc_offset.attrs = {
        "long_name": "the station's local standard time offset from UTC, which "
        "AWDB reports in",
        "units": "hours",
    }
    for name, value in attrs.items():
        if name in ds:
            ds[name].attrs = value
    return ds


# CLI entrypoint


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: tuple[str, ...] = (SOIL_MOISTURE,),
) -> xr.Dataset:
    """Every hourly report of *variables* at NRCS stations in *aoi*.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date/time, UTC (parseable by pandas.Timestamp).
    end_date : str
        End date/time, UTC (parseable by pandas.Timestamp).
    variables : list of str
        AWDB element codes, e.g. ``SMS STO WTEQ``: ``SMS`` volumetric soil
        moisture, ``STO`` soil temperature, ``WTEQ`` snow water equivalent,
        ``SNWD`` snow depth, ``PREC`` water-year precipitation, ``TOBS`` air
        temperature.

    Returns
    -------
    xarray.Dataset
        A CF contiguous ragged array; :func:`gaia_cli.stations.to_grid` and
        :func:`gaia_cli.stations.to_dataframe` give other views of it.
    """
    log.info("querying NRCS stations from %s", AWDB_API)
    sites = get_stations(aoi, start_date, end_date, variables)
    if sites.empty:
        raise ValueError(f"no NRCS station with {', '.join(variables)} in the AOI")
    log.info("fetching hourly data for %d station(s)", len(sites))
    data = get_data(sites.stationTriplet.to_list(), start_date, end_date, variables)
    ds = construct_ragged(sites, data, start_date, end_date)
    ds.attrs = source_attrs(ds)
    return ds
