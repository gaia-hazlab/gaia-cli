from __future__ import annotations

from enum import Enum
from typing import Annotated

import geopandas as gpd
import icechunk
import xarray as xr
from cyclopts import Parameter

from . import io

DEFAULT_BUCKET = "dynamical-noaa-hrrr"
DEFAULT_PREFIX = "noaa-hrrr-analysis/v0.2.0.icechunk/"
DEFAULT_REGION = "us-west-2"


class HRRRVariable(str, Enum):
    """Available variables in the HRRR analysis dataset (dynamical.org)."""

    categorical_freezing_rain_surface = "categorical_freezing_rain_surface"
    categorical_ice_pellets_surface = "categorical_ice_pellets_surface"
    categorical_rain_surface = "categorical_rain_surface"
    categorical_snow_surface = "categorical_snow_surface"
    composite_reflectivity = "composite_reflectivity"
    dew_point_temperature_2m = "dew_point_temperature_2m"
    downward_long_wave_radiation_flux_surface = (
        "downward_long_wave_radiation_flux_surface"
    )
    downward_short_wave_radiation_flux_surface = (
        "downward_short_wave_radiation_flux_surface"
    )
    geopotential_height_cloud_ceiling = "geopotential_height_cloud_ceiling"
    percent_frozen_precipitation_surface = "percent_frozen_precipitation_surface"
    precipitable_water_atmosphere = "precipitable_water_atmosphere"
    precipitation_surface = "precipitation_surface"
    pressure_reduced_to_mean_sea_level = "pressure_reduced_to_mean_sea_level"
    pressure_surface = "pressure_surface"
    relative_humidity_2m = "relative_humidity_2m"
    snow_area_fraction_surface = "snow_area_fraction_surface"
    snow_thickness_surface = "snow_thickness_surface"
    snow_water_equivalent_surface = "snow_water_equivalent_surface"
    snowfall_surface = "snowfall_surface"
    temperature_2m = "temperature_2m"
    total_cloud_cover_atmosphere = "total_cloud_cover_atmosphere"
    wind_gust_surface = "wind_gust_surface"
    wind_u_10m = "wind_u_10m"
    wind_u_80m = "wind_u_80m"
    wind_v_10m = "wind_v_10m"
    wind_v_80m = "wind_v_80m"


def open_hrrr(
    bucket: str = DEFAULT_BUCKET,
    prefix: str = DEFAULT_PREFIX,
    region: str = DEFAULT_REGION,
) -> xr.Dataset:
    """Open the HRRR analysis dataset from Icechunk (dynamical.org).

    Parameters
    ----------
    bucket : str
        S3 bucket hosting the Icechunk repository.
    prefix : str
        Key prefix within the bucket.
    region : str
        AWS region.

    Returns
    -------
    xarray.Dataset
    """
    storage = icechunk.s3_storage(
        bucket=bucket, prefix=prefix, region=region, anonymous=True
    )
    repo = icechunk.Repository.open(storage)
    session = repo.readonly_session("main")
    return xr.open_zarr(session.store, chunks=None)


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])] = None,
    start_date: Annotated[str, Parameter(name=["--start", "-s"])] = None,
    end_date: Annotated[str, Parameter(name=["--end", "-e"])] = None,
    output_path: Annotated[str, Parameter(name=["--output", "-o"])] = None,
    variable: Annotated[
        HRRRVariable, Parameter(name=["--variable", "-v"])
    ] = HRRRVariable.precipitation_surface,
):
    """Stage HRRR analysis data clipped to an AOI.
    The HRRR precipitation rate (kg m-2 s-1) is converted to hourly
    accumulation (mm).  Requires AWS_PROFILE environment variable if
    writing to S3.

    Parameters
    ----------
    vectorPath : str
        Path to a vector file (e.g. GeoJSON, Shapefile) defining the AOI.
    start_date : str
        Start date (inclusive), e.g. '2025-12-01'.
    end_date : str
        End date (inclusive), e.g. '2025-12-20'.
    output_path : str
        Path to write the resulting Zarr store.
    variable : HRRRVariable
        HRRR variable name to extract. Default ``precipitation-surface``.

    Examples
    --------
    Stage HRRR precipitation for Skagit Basin December 2025:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage hrrr -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit_hrrr.zarr
    """
    io.check_s3_env(output_path)

    aoi = io.load_aoi(vectorPath)
    clipped = load(aoi, start_date, end_date, variable)
    io.save_zarr(clipped, output_path)

    print("DONE!")


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variable: HRRRVariable = HRRRVariable.precipitation_surface,
) -> xr.DataArray:
    """Load HRRR analysis data clipped to *aoi*.

    Precipitation rate (kg m⁻² s⁻¹) is converted to hourly accumulation (mm).

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    variable : HRRRVariable
        HRRR variable name to extract.

    Returns
    -------
    xarray.DataArray
    """
    ds = open_hrrr()
    da = ds[variable.value].sel(time=slice(start_date, end_date))

    # HRRR is Lambert Conformal → clip_box, reproject, then clip to polygon
    clipped = io.clip_to_aoi(da, aoi, reproject=True)

    # Convert precipitation rate (kg m⁻² s⁻¹ ≡ mm/s) → hourly accumulation (mm)
    if variable == HRRRVariable.precipitation_surface:
        clipped = clipped.rename("precipitation")
        clipped = clipped * 3600
        clipped.attrs = {
            "long_name": "Hourly precipitation",
            "units": "mm",
            "comment": "Converted from precipitation_surface (kg m-2 s-1) × 3600 s",
        }

    return clipped
