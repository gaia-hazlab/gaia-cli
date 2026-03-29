from __future__ import annotations

import geopandas as gpd
import numpy as np
import odc.stac
import pystac
import xarray as xr
from cyclopts import Parameter
from typing import Annotated
import rioxarray

from . import io

DEFAULT_CATALOG = "https://raw.githubusercontent.com/gaia-hazlab/prism-stac/refs/heads/main/stac/catalog.json"

# GDAL settings for efficient remote reading
odc.stac.configure_rio(cloud_defaults=True)

def open_prism(
    catalog_url: str = DEFAULT_CATALOG,
    start: str | None = None,
    end: str | None = None,
) -> xr.Dataset:
    """Open PRISM daily precipitation from a STAC catalog.

    Parameters
    ----------
    catalog_url : str
        URL of the PRISM STAC catalog JSON.
    start : str, optional
        Start date for time slice (inclusive).
    end : str, optional
        End date for time slice (inclusive).

    Returns
    -------
    xarray.Dataset
    """
    cat = pystac.read_file(catalog_url)
    items = list(cat.get_all_items())

    ds = odc.stac.load(items,
                       bands=['ppt'],
                       #crs='EPSG:4326', # Revisit this, we're overwriting EPSG 4269 -> 4326 which should be fine for coarse data
                       chunks={})

    ds = ds.rename({"ppt": "precipitation", "longitude":"x", "latitude":"y"})
    ds["precipitation"].attrs["units"] = "mm"
    # Use nans instead of -9999
    ds = ds.where(ds.precipitation != ds.precipitation.rio.nodata)
    ds["precipitation"].rio.write_nodata(np.nan, inplace=True)

    if start or end:
        ds = ds.sel(time=slice(start, end))

    return ds


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])] = None,
    start_date: Annotated[str, Parameter(name=["--start", "-s"])] = None,
    end_date: Annotated[str, Parameter(name=["--end", "-e"])] = None,
    output_path: Annotated[str, Parameter(name=["--output", "-o"])] = None,
    catalog_url: Annotated[str, Parameter(name=["--catalog"])] = DEFAULT_CATALOG,
):
    """Stage PRISM daily *precipitation* clipped to an AOI.

    WARNING: currently this only uses a single variable (precipitation) for 2025-12-01 to 2025-12-31

    Requires AWS_PROFILE environment variable if writing to S3.

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
    catalog_url : str
        URL of the PRISM STAC catalog JSON.

    Examples
    --------
    Stage PRISM for Skagit Basin December 2025:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage prism -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit_prism.zarr
    """
    io.check_s3_env(output_path)

    aoi = io.load_aoi(vectorPath)
    clipped = load(aoi, start_date, end_date, catalog_url)
    io.save_zarr(clipped, output_path)

    print("DONE!")


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    catalog_url: str = DEFAULT_CATALOG,
) -> xr.DataArray:
    """Load PRISM precipitation clipped to *aoi*.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    catalog_url : str
        URL of the PRISM STAC catalog JSON.

    Returns
    -------
    xarray.DataArray
    """
    ds = open_prism(catalog_url, start=start_date, end=end_date)
    da = ds["precipitation"]
    return io.clip_to_aoi(da, aoi, reproject=True).compute()
