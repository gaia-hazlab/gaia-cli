from __future__ import annotations

from typing import Annotated

import geopandas as gpd
import rioxarray
import xarray as xr
from cyclopts import Parameter
from xarray.backends.file_manager import FILE_CACHE

from . import io


def open_nlcd(year: int = 2021) -> xr.DataArray:
    url = (
        "https://usgs.osn.mghpcc.org/hytest/nlcd/annual-nlcd-cu-c1v1/mosaic/"
        f"Annual_NLCD_LndCov_{year}_CU_C1V1.tif"
    )
    # rioxarray's stub is stricter than its runtime API: chunks="auto" is a valid dask spec,
    # and without group/variable args this always returns a single DataArray, not the full union
    da = rioxarray.open_rasterio(url, chunks="auto", cache=False)  # ty: ignore[invalid-argument-type]
    da = da.squeeze("band", drop=True).rename("landcover")  # ty: ignore[unresolved-attribute, invalid-argument-type]
    da.attrs["long_name"] = "Annual NLCD land cover"
    da.attrs["year"] = year
    return da  # ty: ignore[invalid-return-type]


def load(aoi: gpd.GeoDataFrame, year: int = 2021) -> xr.DataArray:
    da = open_nlcd(year)
    aoi_native = aoi.to_crs(da.rio.crs)
    minx, miny, maxx, maxy = aoi_native.total_bounds

    small = da.rio.clip_box(minx=minx, miny=miny, maxx=maxx, maxy=maxy)
    clipped = small.rio.clip(aoi_native.geometry)

    return clipped


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    year: Annotated[int, Parameter(name=["--year", "-y"])] = 2021,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
):
    """Stage NLCD land cover clipped to an AOI."""
    io.check_output(output_path, overwrite)

    aoi = io.load_aoi(vectorPath)
    clipped = load(aoi, year=year)
    io.save_zarr(clipped, output_path, overwrite=overwrite)

    # Close cached remote rasterio/GDAL handles opened by rioxarray so their
    # cleanup doesn't run during interpreter shutdown, where it can crash
    # sys.excepthook on Python 3.14.
    FILE_CACHE.clear()

    print("DONE!")
