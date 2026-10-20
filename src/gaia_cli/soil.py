from __future__ import annotations

import logging
from enum import Enum
from typing import Annotated

import geopandas as gpd
import numpy as np
import rasterio
import rioxarray
import xarray as xr
from cyclopts import Parameter

from . import io, provenance, regrid

log = logging.getLogger(__name__)

# Written by hand (no STAC Collection): the bounds of the SOLUS100 COGs in
# EPSG:4326, which hold the STATSGO COG's, and the release date of SOLUS's layer
# table. Static: the product is an option, not a window. See dev/soil.md.
EXTENT = {
    "spatial": {"bbox": [[-130.5738, 21.1114, -62.664, 52.3929]]},
    "temporal": {"interval": [["2023-12-15T00:00:00Z", "2023-12-15T00:00:00Z"]]},
}
KIND = "raster"
TIMED = False
# `gaia align` default, as in the landlab workflow contract. statsgo's kffact
# sets its own (nearest)
REGRID_METHOD = "bilinear"
# Each product is its own group, soil/<product>, keyed so in stage files
VARIANT = "product"


class Product(str, Enum):
    """Soil property maps."""

    solus100 = "solus100"
    statsgo = "statsgo"


SOLUS_URL = "https://storage.googleapis.com/solus100pub"

# The layers landlab reads (dev/soil.md). The COGs hold integers times
# `scalar`, which only Final_Layer_Table_20231215.csv records: the files' own
# scale is 1.
SOLUS_LAYERS: dict[str, dict] = {
    "anylithicdpt": {
        "file": "anylithicdpt_cm_p.tif",
        "scalar": 1,
        "units": "cm",
        "long_name": "Depth to bedrock (lithic or paralithic)",
    },
    "cec7": {
        "file": "cec7_0_cm_p.tif",
        "scalar": 10,
        "units": "meq/100g",
        "long_name": "Cation exchange capacity (pH 7)",
    },
    "claytotal": {
        "file": "claytotal_0_cm_p.tif",
        "scalar": 1,
        "units": "percent",
        "long_name": "Clay content (mass fraction)",
    },
    "dbovendry": {
        "file": "dbovendry_0_cm_p.tif",
        "scalar": 100,
        "units": "g cm-3",
        "long_name": "Bulk density (oven dry)",
    },
    "ph1to1h2o": {
        "file": "ph1to1h2o_0_cm_p.tif",
        "scalar": 100,
        "units": "1",
        "long_name": "pH (1:1 soil:water)",
    },
    "sandtotal": {
        "file": "sandtotal_0_cm_p.tif",
        "scalar": 1,
        "units": "percent",
        "long_name": "Sand content (mass fraction)",
    },
    "silttotal": {
        "file": "silttotal_0_cm_p.tif",
        "scalar": 1,
        "units": "percent",
        "long_name": "Silt content (mass fraction)",
    },
}


# USGS's 30 m rasterization of the 1995 STATSGO KFFACT polygons. The public S3
# copy: the ScienceBase file link (catalog/file/get/...?name=) is a 404 since
# the move to sciencebase.usgs.gov. NAD27 Albers (EPSG:5069), as in the source
STATSGO_ITEM = "6750c172d34ed8d3858534d8"
STATSGO_URL = (
    "https://prod-is-usgs-sb-prod-publish.s3.amazonaws.com/"
    f"{STATSGO_ITEM}/statsgo-KFFACT.tif"
)
# The source's missing-value code: map units with no K factor, mostly water
STATSGO_MISSING = -0.1


def _url(layer: str) -> str:
    return f"{SOLUS_URL}/{SOLUS_LAYERS[layer]['file']}"


def source_attrs(product: Product = Product.solus100) -> dict:
    """Provenance attributes for the *product* soil maps."""
    if Product(product) == Product.statsgo:
        return _statsgo_attrs()
    return {
        "title": "SOLUS100 soil property maps",
        "description": "Soil Landscapes of the United States (SOLUS100): "
        "predicted soil properties, conterminous US, 100 m. Surface (0 cm) "
        "values, plus depth to bedrock.",
        "license": "CC-BY-4.0",
        "providers": [
            {
                "name": "USDA Agricultural Research Service and Natural "
                "Resources Conservation Service",
                "roles": ["producer", "licensor"],
                "url": "https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_Soil_Landscapes_of_the_United_States_100-meter_SOLUS100_soil_property_maps_project_repository/25033856",
            },
            {
                "name": "Google Cloud Storage (solus100pub)",
                "roles": ["host"],
                "url": f"{SOLUS_URL}/index.html",
            },
        ],
        "sci:doi": "10.15482/USDA.ADC/25033856",
        "sci:citation": "Nauman, T. W., et al. (2024). Soil landscapes of the "
        "United States (SOLUS): Developing predictive soil property maps of the "
        "conterminous United States using hybrid training sets. Soil Science "
        "Society of America Journal. https://doi.org/10.1002/saj2.20769",
        "processing:lineage": "clipped to the AOI on the native grid; stored "
        "integers divided by the layer table's scalar",
        "links": [
            *(
                provenance.link("via", _url(k), type="image/tiff; application=geotiff")
                for k in SOLUS_LAYERS
            ),
            provenance.link(
                "describedby",
                f"{SOLUS_URL}/Final_Layer_Table_20231215.csv",
                type="text/csv",
            ),
            provenance.link("license", "https://creativecommons.org/licenses/by/4.0/"),
            provenance.link("cite-as", "https://doi.org/10.15482/USDA.ADC/25033856"),
        ],
    }


def _statsgo_attrs() -> dict:
    return {
        "title": "STATSGO K-factor (KFFACT)",
        "description": "Soil erodibility factor of the universal soil loss "
        "equation from STATSGO (1:250,000, October 1994), averaged over each map "
        "unit's layers and components (Schwarz and Alexander, 1995), rasterized "
        "at 30 m by USGS. Conterminous US.",
        "license": "other",  # U.S. public domain has no SPDX identifier
        "providers": [
            {
                "name": "U.S. Geological Survey",
                "roles": ["producer", "licensor", "processor", "host"],
                "url": f"https://www.sciencebase.gov/catalog/item/{STATSGO_ITEM}",
            },
            {
                "name": "USDA Natural Resources Conservation Service (STATSGO)",
                "roles": ["producer"],
                "url": "https://doi.org/10.5066/P94JAULO",
            },
        ],
        "sci:doi": "10.5066/P13WAPYV",
        "sci:citation": "King, J.M. (2025). STATSGO Cloud Optimized GeoTIFFs for "
        "the Continental US. U.S. Geological Survey data release. "
        "https://doi.org/10.5066/P13WAPYV",
        "sci:publications": [
            {
                "doi": "10.5066/P94JAULO",
                "citation": "Schwarz, G.E. and Alexander, R.B. (1995). Soils data "
                "for the Conterminous United States Derived from the NRCS State "
                "Soil Geographic (STATSGO) Data Base. U.S. Geological Survey data "
                "release. https://doi.org/10.5066/P94JAULO",
            }
        ],
        "processing:lineage": "clipped to the AOI on the native grid (NAD27 "
        "Albers); the source's missing-value code -0.1 set to NaN",
        "links": [
            provenance.link("via", STATSGO_URL, type="image/tiff; application=geotiff"),
            provenance.link(
                "describedby",
                f"https://www.sciencebase.gov/catalog/item/{STATSGO_ITEM}",
                type="text/html",
            ),
            provenance.link(
                "license",
                "https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits",
                type="text/html",
            ),
            provenance.link("cite-as", "https://doi.org/10.5066/P13WAPYV"),
        ],
    }


def open_layer(layer: str) -> xr.DataArray:
    url = _url(layer)
    log.info("SOLUS: opening %s", url)
    # Lazy but not dask: a clip then reads one window of the COG, instead of every
    # whole dask chunk it touches. masked: nodata (255 or 65535) becomes NaN
    da = rioxarray.open_rasterio(url, cache=False, masked=True)
    return da.squeeze("band", drop=True).rename(layer)  # ty: ignore[unresolved-attribute, invalid-argument-type, invalid-return-type]


def _decode(raw: xr.DataArray, layer: str) -> xr.DataArray:
    """*raw* in physical units, packed back to the source integers on write."""
    spec = SOLUS_LAYERS[layer]
    scale = 1 / spec["scalar"]
    # float64 times scale, as xarray decodes the store, so that load() matches it
    da = raw.astype("float64") * scale
    da.attrs = {
        "long_name": spec["long_name"],
        "units": spec["units"],
        "depth": "all" if layer == "anylithicdpt" else "0 cm",
    }
    da.encoding = {
        "dtype": raw.encoding["dtype"],
        "_FillValue": raw.encoding["_FillValue"],
        "scale_factor": scale,
    }
    return da


def _load_statsgo(aoi: gpd.GeoDataFrame) -> xr.DataArray:
    log.info("STATSGO: opening %s", STATSGO_URL)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        raw = rioxarray.open_rasterio(STATSGO_URL, cache=False, masked=True)
        raw = io.clip_to_aoi(raw.squeeze("band", drop=True), aoi).load()  # ty: ignore[unresolved-attribute, invalid-argument-type]
    # float32, as the source; the code is float32 -0.1 too
    da = raw.where(raw != raw.dtype.type(STATSGO_MISSING)).rename("kffact")
    da.attrs = {
        "long_name": "Soil erodibility factor (K-factor), STATSGO KFFACT",
        # TODO: check. USGS's 2025 metadata defines KFFACT as the saturated
        # hydraulic conductivity of the fine soil in inches per hour, but the 1995
        # source calls it the USLE K factor, in (t ac h)/(100 ac ft tonf in).
        # See dev/soil.md
        "units": "inches per hour",
        "depth": "all",
        "comment": "Thickness-weighted over a map unit's layers, then "
        "averaged over its components. The source's -0.1 (no K factor, mostly "
        "water) is NaN.",
        # A map of polygons: interpolating between them makes values that
        # aren't in it
        regrid.METHOD_ATTR: "nearest",
    }
    da.encoding = {"dtype": "float32", "_FillValue": np.float32("nan")}
    return da


def load(aoi: gpd.GeoDataFrame, product: Product = Product.solus100) -> xr.Dataset:
    """*product* soil properties clipped to *aoi* on their native grid, with
    provenance attributes.
    """
    if Product(product) == Product.statsgo:
        ds = _load_statsgo(aoi).to_dataset()
        ds.attrs = source_attrs(product)
        return ds
    # Each layer is one COG, so GDAL needn't list the bucket for sidecar files.
    # rioxarray reads lazily, so load before leaving the Env.
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        layers = [
            _decode(io.clip_to_aoi(open_layer(layer), aoi).load(), layer)
            for layer in SOLUS_LAYERS
        ]
    ds = xr.Dataset({da.name: da for da in layers})
    ds.attrs = source_attrs(product)
    return ds


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    product: Annotated[Product, Parameter(name=["--product", "-p"])] = Product.solus100,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
):
    """Stage soil property maps clipped to an AOI.

    Parameters
    ----------
    vectorPath : str
        Path to a vector file (e.g. GeoJSON, Shapefile) defining the AOI.
    output_path : str
        Output folder for the Zarr store (named after the folder), its STAC
        Item and the resolved stage.yaml.
    product : Product
        ``solus100`` (USDA SOLUS100, CONUS, 100 m): the seven layers landlab
        reads. ``statsgo`` (USGS, CONUS, 30 m, from 1:250,000 polygons): the
        K-factor, ``kffact``. Default solus100.
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
        "soil",
        vectorPath,
        None,
        None,
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
        product=product,
    )
