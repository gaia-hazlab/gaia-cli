from __future__ import annotations

from pathlib import Path
from typing import Annotated

from bmi_topography import Topography
from cyclopts import Parameter
import geopandas as gpd
import rioxarray
import xarray as xr
from xarray.backends.file_manager import FILE_CACHE

from . import io


def open_dem(
    aoi: gpd.GeoDataFrame,
    dem_type: str = "USGS10m",
    buffer_deg: float = 0.05,
    api_key: str | None = None,
) -> xr.DataArray:
    aoi_4326 = aoi.to_crs(epsg=4326)
    west, south, east, north = aoi_4326.total_bounds

    topo = Topography(
        dem_type=dem_type,
        south=south - buffer_deg,
        north=north + buffer_deg,
        west=west - buffer_deg,
        east=east + buffer_deg,
        output_format="GTiff",
        api_key=api_key,
    )

    dem_path = Path(topo.fetch())
    da = rioxarray.open_rasterio(dem_path, chunks="auto", cache=False)
    da = da.squeeze("band", drop=True).rename("elevation")
    da.attrs["long_name"] = f"{dem_type} elevation"
    da.attrs["units"] = "m"
    da.attrs["dem_type"] = dem_type
    return da


def load(
    aoi: gpd.GeoDataFrame,
    dem_type: str = "USGS10m",
    buffer_deg: float = 0.05,
    api_key: str | None = None,
) -> xr.DataArray:
    da = open_dem(aoi, dem_type=dem_type, buffer_deg=buffer_deg, api_key=api_key)
    aoi_native = aoi.to_crs(da.rio.crs)
    return da.rio.clip(aoi_native.geometry)


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])] = None,
    output_path: Annotated[str, Parameter(name=["--output", "-o"])] = None,
    dem_type: Annotated[str, Parameter(name=["--dem-type"])] = "USGS10m",
    buffer_deg: Annotated[float, Parameter(name=["--buffer-deg"])] = 0.05,
    api_key: Annotated[str | None, Parameter(name=["--api-key"])] = None,
):
    """Stage an OpenTopography DEM clipped to an AOI."""
    io.check_s3_env(output_path)

    aoi = io.load_aoi(vectorPath)
    clipped = load(aoi, dem_type=dem_type, buffer_deg=buffer_deg, api_key=api_key)
    io.save_zarr(clipped, output_path)
    FILE_CACHE.clear()

    print("DONE!")
