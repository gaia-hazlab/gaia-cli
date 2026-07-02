from __future__ import annotations

import geopandas as gpd
import rioxarray
import xarray as xr
from cyclopts import Parameter
from typing import Annotated

from . import io

def open_nlcd(year: int = 2021) -> xr.DataArray:
    url = (
        "https://usgs.osn.mghpcc.org/hytest/nlcd/annual-nlcd-cu-c1v1/mosaic/"
        f"Annual_NLCD_LndCov_{year}_CU_C1V1.tif"
    )
    da = rioxarray.open_rasterio(url, chunks="auto", cache=False)
    da = da.squeeze("band", drop=True).rename("landcover")
    da.attrs["long_name"] = "Annual NLCD land cover"
    da.attrs["year"] = year
    return da


def load(aoi: gpd.GeoDataFrame, year: int = 2021) -> xr.DataArray:
    da = open_nlcd(year)
    aoi_native = aoi.to_crs(da.rio.crs)
    minx, miny, maxx, maxy = aoi_native.total_bounds

    small = da.rio.clip_box(minx=minx, miny=miny, maxx=maxx, maxy=maxy)
    clipped = small.rio.clip(aoi_native.geometry)

    return clipped


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])] = None,
    output_path: Annotated[str, Parameter(name=["--output", "-o"])] = None,
    year: Annotated[int, Parameter(name=["--year", "-y"])] = 2021,
):
    """Stage NLCD land cover clipped to an AOI."""
    io.check_s3_env(output_path)

    aoi = io.load_aoi(vectorPath)
    clipped = load(aoi, year=year)
    io.save_zarr(clipped, output_path)

    print("DONE!")
