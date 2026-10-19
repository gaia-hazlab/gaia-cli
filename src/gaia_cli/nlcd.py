from __future__ import annotations

import logging
from typing import Annotated

import geopandas as gpd
import rioxarray
import xarray as xr
from cyclopts import Parameter
from xarray.backends.file_manager import FILE_CACHE

from . import io, provenance

log = logging.getLogger(__name__)


def _url(year: int) -> str:
    return (
        "https://usgs.osn.mghpcc.org/hytest/nlcd/annual-nlcd-cu-c1v1/mosaic/"
        f"Annual_NLCD_LndCov_{year}_CU_C1V1.tif"
    )


def source_attrs(year: int) -> dict:
    """Provenance attributes for the *year* NLCD land cover mosaic."""
    return {
        "title": "Annual NLCD land cover",
        "description": "USGS Annual National Land Cover Database (Collection 1 "
        "v1) land cover, conterminous US, 30 m.",
        "license": "other",  # U.S. public domain has no SPDX identifier
        "providers": [
            {
                "name": "U.S. Geological Survey",
                "roles": ["producer", "licensor"],
                "url": "https://www.usgs.gov/centers/eros/science/annual-national-land-cover-database",
            },
            {
                "name": "USGS HyTEST",
                "roles": ["host"],
                "url": "https://usgs.osn.mghpcc.org",
            },
        ],
        "version": "C1V1",
        "start_datetime": f"{year}-01-01T00:00:00Z",
        "end_datetime": f"{year}-12-31T23:59:59Z",
        "processing:lineage": "clipped to the AOI on the native grid",
        "links": [
            provenance.link("via", _url(year), type="image/tiff; application=geotiff"),
            provenance.link(
                "license",
                "https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits",
                type="text/html",
            ),
        ],
    }


def open_nlcd(year: int = 2021) -> xr.DataArray:
    url = _url(year)
    log.info("NLCD: opening %s", url)
    # rioxarray's stub is stricter than its runtime API: chunks="auto" is a valid dask spec,
    # and without group/variable args this always returns a single DataArray, not the full union
    da = rioxarray.open_rasterio(url, chunks="auto", cache=False)  # ty: ignore[invalid-argument-type]
    da = da.squeeze("band", drop=True).rename("landcover")  # ty: ignore[unresolved-attribute, invalid-argument-type]
    da.attrs["long_name"] = "Annual NLCD land cover"
    da.attrs["year"] = year
    return da  # ty: ignore[invalid-return-type]


def load(aoi: gpd.GeoDataFrame, year: int = 2021) -> xr.DataArray:
    return io.clip_to_aoi(open_nlcd(year), aoi)


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
    ds = clipped.to_dataset()
    ds.attrs = source_attrs(year)
    io.save_zarr(ds, output_path, overwrite=overwrite)

    # Close cached remote rasterio/GDAL handles opened by rioxarray so their
    # cleanup doesn't run during interpreter shutdown, where it can crash
    # sys.excepthook on Python 3.14.
    FILE_CACHE.clear()

    log.info("Done")
