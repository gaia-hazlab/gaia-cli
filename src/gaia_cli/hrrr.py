from __future__ import annotations

import logging
from enum import Enum
from typing import Annotated

import geopandas as gpd
import icechunk
import xarray as xr
from cyclopts import Parameter

from . import io, provenance

log = logging.getLogger(__name__)

DEFAULT_BUCKET = "dynamical-noaa-hrrr"
DEFAULT_PREFIX = "noaa-hrrr-analysis/v0.2.0.icechunk/"
DEFAULT_REGION = "us-west-2"
STAC_COLLECTION = "https://stac.dynamical.org/noaa-hrrr-analysis/collection.json"
# Copied from STAC_COLLECTION
EXTENT = {
    "spatial": {"bbox": [[-134.095474, 21.138123, -60.917194, 52.615654]]},
    "temporal": {"interval": [["2014-10-01T00:00:00Z", None]]},
}
KIND = "raster"
TIMED = True
# `gaia align` default: area-weighted, so totals and means are preserved
REGRID_METHOD = "conservative"


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
        With the Icechunk snapshot read in ``attrs["gaia:icechunk-snapshot"]``
        and the repository URL in ``encoding["source"]``.
    """
    # Pin the AWS endpoint so an AWS_PROFILE with a custom endpoint_url
    # (e.g. for writing output elsewhere) doesn't redirect this public read.
    storage = icechunk.s3_storage(
        bucket=bucket,
        prefix=prefix,
        region=region,
        anonymous=True,
        endpoint_url=f"https://s3.{region}.amazonaws.com",
    )
    log.info("opening dynamical.org Icechunk repo s3://%s/%s", bucket, prefix)
    repo = icechunk.Repository.open(storage)
    session = repo.readonly_session("main")
    ds = xr.open_zarr(session.store, chunks=None)
    ds.attrs["gaia:icechunk-snapshot"] = session.snapshot_id
    ds.encoding["source"] = f"s3://{bucket}/{prefix}"
    return ds


def source_attrs(ds: xr.Dataset) -> dict:
    """Provenance attributes for data read from *ds* (from :func:`open_hrrr`).

    Copied from dynamical.org's STAC collection, except ``version`` and
    ``license``, which come from the store itself: the collection describes
    dynamical's latest version, which may be newer than ``DEFAULT_PREFIX``.
    """
    attrs = provenance.from_stac_collection(STAC_COLLECTION)
    for key, upstream in (("version", "dataset_version"), ("license", "license")):
        if upstream in ds.attrs:
            attrs[key] = ds.attrs[upstream]
    if "gaia:icechunk-snapshot" in ds.attrs:
        attrs["gaia:icechunk-snapshot"] = ds.attrs["gaia:icechunk-snapshot"]
    if "source" in ds.encoding:
        attrs["links"].append(
            provenance.link("via", ds.encoding["source"], title="Icechunk repository")
        )
    return attrs


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    start_date: Annotated[str, Parameter(name=["--start", "-s"])],
    end_date: Annotated[str, Parameter(name=["--end", "-e"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    variable: Annotated[
        HRRRVariable, Parameter(name=["--variable", "-v"])
    ] = HRRRVariable.precipitation_surface,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
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
        Output folder for the Zarr store (named after the folder), its STAC
        Item and the resolved stage.yaml.
    variable : HRRRVariable
        HRRR variable name to extract. Default ``precipitation-surface``.
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
    Stage HRRR precipitation for Skagit Basin December 2025:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage hrrr -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-hrrr
    """
    from .stagefile import stage_dataset  # stagefile imports this module

    stage_dataset(
        "hrrr",
        vectorPath,
        start_date,
        end_date,
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
        variable=variable,
    )


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variable: HRRRVariable = HRRRVariable.precipitation_surface,
) -> xr.Dataset:
    """Load HRRR analysis data clipped to *aoi* on its native Lambert Conformal grid.

    Precipitation rate (kg m⁻² s⁻¹) is converted to hourly accumulation (mm).
    Averaged variables (``step_type="avg"``) cover the hour ending at ``time``,
    i.e. ``(time - 1h, time]``; instantaneous ones are valid at ``time``
    (see dev/hrrr.md).

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
    xarray.Dataset
        The variable, plus CF ``time_bounds`` for averaged variables, with
        provenance attributes (see :func:`source_attrs`).
    """
    ds = open_hrrr()
    da = ds[variable.value].sel(time=slice(start_date, end_date))
    step_type = da.attrs.get("step_type")

    clipped = io.clip_to_aoi(da, aoi)
    log.info(
        "HRRR %s: reading %d hour(s) x %d x %d cells (%s to %s)",
        variable.value,
        clipped.sizes["time"],
        clipped.sizes["y"],
        clipped.sizes["x"],
        str(clipped.time.values[0])[:16],
        str(clipped.time.values[-1])[:16],
    )
    # Read now (rather than lazily during the Zarr write) so timings are attributable
    clipped = clipped.load()

    # Convert precipitation rate (kg m⁻² s⁻¹ ≡ mm/s) → hourly accumulation (mm)
    if variable == HRRRVariable.precipitation_surface:
        clipped = clipped.rename("precipitation")
        clipped = clipped * 3600
        clipped.attrs = {
            "long_name": "Precipitation",  # period given by time_bounds
            "units": "mm",
            "comment": "Converted from precipitation_surface (kg m-2 s-1) × 3600 s",
        }
        cell_methods = "time: sum"
    else:
        cell_methods = "time: mean"

    if step_type == "avg":
        out = io.with_time_bounds(
            clipped,
            start_hours=-1,
            end_hours=0,
            cell_methods=cell_methods,
        )
    else:
        clipped.attrs["cell_methods"] = "time: point"
        out = clipped.to_dataset(name=clipped.name or variable.value)
    out.attrs = {
        **source_attrs(ds),
        "processing:lineage": f"{variable.value} clipped to the AOI on its native grid",
        **provenance.time_extent(out),
    }
    return out
