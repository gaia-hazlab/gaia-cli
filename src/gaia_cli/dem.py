from __future__ import annotations

import logging
from enum import Enum
from typing import Annotated, Any

import geopandas as gpd
import rasterio
import rioxarray
import xarray as xr
from cyclopts import Parameter

from . import io, provenance

log = logging.getLogger(__name__)

# Written by hand: Cop30 is global, and covers every usgs10m AOI. Static: the
# product is an option, not a window. See docs/dem.md.
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [["1923-01-01T00:00:00Z", None]]},
}
KIND = "raster"
TIMED = False
# Each product is its own group, dem/<product>, keyed so in stage files
VARIANT = "product"


class Product(str, Enum):
    """DEM sources, each a public VRT mosaic of COGs."""

    cop30 = "cop30"
    usgs10m = "usgs10m"


# Heights are above a geoid, which differs by product; vertical_datum is the EPSG code
PRODUCTS: dict[Product, dict[str, Any]] = {
    Product.cop30: {
        "url": "https://opentopography.s3.sdsc.edu/raster/COP30/COP30_hh.vrt",
        "title": "Copernicus DEM GLO-30",
        "description": "Copernicus global digital surface model, 1 arc-second "
        "(~30 m), from TanDEM-X acquisitions 2011-2015.",
        "license": "other",
        "license_href": "https://spacedata.copernicus.eu/documents/20126/0/CSCDA_ESA_Mission-specific+Annex.pdf",
        "providers": [
            {
                "name": "European Space Agency",
                "roles": ["producer", "licensor"],
                "url": "https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM",
            },
            {
                "name": "OpenTopography",
                "roles": ["host"],
                "url": "https://opentopography.org",
            },
        ],
        "sci:doi": "10.5069/G9028PQB",
        "sci:citation": "European Space Agency (2024). Copernicus Global Digital "
        "Elevation Model. Distributed by OpenTopography.",
        "start_datetime": "2011-01-01T00:00:00Z",
        "end_datetime": "2015-07-01T00:00:00Z",
        "vertical_datum": "EPSG:3855",  # EGM2008 height
    },
    Product.usgs10m: {
        "url": "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/USGS_Seamless_DEM_13.vrt",
        "title": "USGS 3DEP 1/3 arc-second DEM",
        "description": "USGS 3D Elevation Program seamless bare-earth DEM, "
        "1/3 arc-second (~10 m); conterminous US, Hawaii, territories and "
        "parts of Alaska. Updated continuously from source data of many dates.",
        "license": "other",  # U.S. public domain has no SPDX identifier
        "license_href": "https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits",
        "providers": [
            {
                "name": "U.S. Geological Survey",
                "roles": ["producer", "licensor", "host"],
                "url": "https://www.usgs.gov/3d-elevation-program",
            },
        ],
        "sci:doi": "10.5069/G98K778D",
        "sci:citation": "United States Geological Survey (2021). United States "
        "Geological Survey 3D Elevation Program 1/3 arc-second Digital Elevation "
        "Model. Distributed by OpenTopography.",
        "vertical_datum": "EPSG:5703",  # NAVD88 height
    },
}


def source_attrs(product: Product) -> dict:
    """Provenance attributes for the *product* DEM."""
    p = PRODUCTS[Product(product)]
    attrs = {
        k: p[k]
        for k in (
            "title",
            "description",
            "license",
            "providers",
            "version",
            "sci:doi",
            "sci:citation",
            "start_datetime",
            "end_datetime",
        )
        if k in p
    }
    attrs["processing:lineage"] = "clipped to the AOI on the native grid"
    attrs["links"] = [
        provenance.link("via", p["url"], type="application/xml"),
        provenance.link("license", p["license_href"]),
        provenance.link("cite-as", f"https://doi.org/{p['sci:doi']}"),
    ]
    return attrs


def open_dem(product: Product = Product.cop30) -> xr.DataArray:
    p = PRODUCTS[Product(product)]
    log.info("DEM: opening %s", p["url"])
    # Lazy but not dask: a clip then reads one window of the VRT, instead of every
    # whole dask chunk it touches. masked: usgs10m's -999999 nodata becomes NaN
    # (cop30 sets none)
    da = rioxarray.open_rasterio(p["url"], cache=False, masked=True)
    da = da.squeeze("band", drop=True).rename("elevation")  # ty: ignore[unresolved-attribute, invalid-argument-type]
    da.attrs = {
        "long_name": f"{p['title']} elevation",
        "standard_name": "surface_altitude",
        "units": "m",
        "vertical_datum": p["vertical_datum"],
    }
    return da  # ty: ignore[invalid-return-type]


def load(aoi: gpd.GeoDataFrame, product: Product = Product.cop30) -> xr.Dataset:
    """*product* elevation clipped to *aoi* on its native grid, with provenance attributes.

    Heights stay relative to the product's geoid (``vertical_datum``).
    """
    # Each VRT lists its tiles, so GDAL needn't list directories or probe for
    # sidecar files. rioxarray reads lazily, so load before leaving the Env.
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        ds = io.clip_to_aoi(open_dem(product), aoi).to_dataset().load()
    ds.attrs = source_attrs(product)
    return ds


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    product: Annotated[Product, Parameter(name=["--product", "-p"])] = Product.cop30,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
):
    """Stage a digital elevation model clipped to an AOI.

    Parameters
    ----------
    vectorPath : str
        Path to a vector file (e.g. GeoJSON, Shapefile) defining the AOI.
    output_path : str
        Output folder for the Zarr store (named after the folder), its STAC
        Item and the resolved stage.yaml.
    product : Product
        ``cop30`` (Copernicus GLO-30, global, EGM2008 heights) or ``usgs10m``
        (USGS 3DEP 1/3 arc-second, US, NAVD88 heights). Default cop30.
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
        "dem",
        vectorPath,
        None,
        None,
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
        product=product,
    )
