from __future__ import annotations

import logging
from typing import Annotated

import geopandas as gpd
import rioxarray
import xarray as xr
from cyclopts import Parameter

from . import io, provenance

log = logging.getLogger(__name__)

# Written by hand (no STAC Collection): the mosaic's bounds in EPSG:4326 and the
# years published as Collection 1 v1. Static: the year is an option, not a window.
EXTENT = {
    "spatial": {"bbox": [[-129.2773, 21.8051, -63.1184, 52.9217]]},
    "temporal": {"interval": [["1985-01-01T00:00:00Z", "2024-12-31T23:59:59Z"]]},
}
KIND = "raster"
TIMED = False

# Land cover classes: code -> (name, color), the COG's own color table
CLASSES = {
    11: ("Open Water", "#466b9f"),
    12: ("Perennial Ice/Snow", "#d1def8"),
    21: ("Developed, Open Space", "#dec5c5"),
    22: ("Developed, Low Intensity", "#d99282"),
    23: ("Developed, Medium Intensity", "#eb0000"),
    24: ("Developed, High Intensity", "#ab0000"),
    31: ("Barren Land", "#b3ac9f"),
    41: ("Deciduous Forest", "#68ab5f"),
    42: ("Evergreen Forest", "#1c5f2c"),
    43: ("Mixed Forest", "#b5c58f"),
    52: ("Shrub/Scrub", "#ccb879"),
    71: ("Grassland/Herbaceous", "#dfdfc2"),
    81: ("Pasture/Hay", "#dcd939"),
    82: ("Cultivated Crops", "#ab6c28"),
    90: ("Woody Wetlands", "#b8d9eb"),
    95: ("Emergent Herbaceous Wetlands", "#6c9fb8"),
}


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
    # Lazy but not dask: a clip then reads one window of the COG, instead of every
    # whole dask chunk it touches
    da = rioxarray.open_rasterio(url, cache=False)
    da = da.squeeze("band", drop=True).rename("landcover")  # ty: ignore[unresolved-attribute, invalid-argument-type]
    da.attrs["long_name"] = "Annual NLCD land cover"
    da.attrs["year"] = year
    # CF flags mark it categorical (gaia align uses mode); flag_colors for legends
    da.attrs["flag_values"] = list(CLASSES)
    da.attrs["flag_meanings"] = " ".join(
        name.replace(", ", "_").replace("/", "_").replace(" ", "_")
        for name, _ in CLASSES.values()
    )
    da.attrs["flag_colors"] = " ".join(color for _, color in CLASSES.values())
    return da  # ty: ignore[invalid-return-type]


def load(aoi: gpd.GeoDataFrame, year: int = 2021) -> xr.Dataset:
    """NLCD land cover for *year* clipped to *aoi*, with provenance attributes."""
    ds = io.clip_to_aoi(open_nlcd(year), aoi).to_dataset()
    ds.attrs = source_attrs(year)
    return ds


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    year: Annotated[int, Parameter(name=["--year", "-y"])] = 2021,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
):
    """Stage NLCD land cover clipped to an AOI.

    Parameters
    ----------
    vectorPath : str
        Path to a vector file (e.g. GeoJSON, Shapefile) defining the AOI.
    output_path : str
        Output folder for the Zarr store (named after the folder), its STAC
        Item and the resolved stage.yaml.
    year : int
        Land cover year (1985-2024). Default 2021.
    dry_run : bool
        Validate the flags, read the AOI and print the resolved configuration
        without fetching any data.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.
    skip_existing : bool
        Add to an existing store, unless it already has this dataset's group.
    """
    from .stagefile import stage_dataset  # stagefile imports this module

    stage_dataset(
        "nlcd",
        vectorPath,
        None,
        None,
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
        year=year,
    )
