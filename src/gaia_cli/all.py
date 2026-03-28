"""``gaia stage all`` – load every precipitation product and save as a DataTree."""

from __future__ import annotations

import xarray as xr
import xvec
from cyclopts import Parameter
from typing import Annotated

from . import io
from .hrrr import HRRRVariable
from . import hrrr as _hrrr
from . import prism as _prism
from . import synoptic as _synoptic


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])] = None,
    start_date: Annotated[str, Parameter(name=["--start", "-s"])] = None,
    end_date: Annotated[str, Parameter(name=["--end", "-e"])] = None,
    output_path: Annotated[str, Parameter(name=["--output", "-o"])] = None,
    hrrr_variable: Annotated[
        HRRRVariable, Parameter(name=["--hrrr-variable"])
    ] = HRRRVariable.precipitation_surface,
):
    """Stage **all** precipitation products clipped to an AOI as a single DataTree.

    Loads Synoptic station data, PRISM gridded data, and HRRR analysis
    data, combines them into an ``xarray.DataTree``, and writes the result
    as a single Zarr store.

    Requires ``SYNOPTIC_TOKEN`` and (for S3 output) ``AWS_PROFILE``
    environment variables.

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
    hrrr_variable : HRRRVariable
        HRRR variable to extract. Default ``precipitation-surface``.

    Examples
    --------
    Stage all products for Skagit Basin December 2025:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage all -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit_all.zarr
    """
    io.check_s3_env(output_path)
    aoi = io.load_aoi(vectorPath)

    # ---- load each product ------------------------------------------------
    print("Loading Synoptic station data …")
    da_synoptic = _synoptic.load(aoi, start_date, end_date)

    print("Loading PRISM gridded data …")
    da_prism = _prism.load(aoi, start_date, end_date)

    print("Loading HRRR analysis data …")
    da_hrrr = _hrrr.load(aoi, start_date, end_date, variable=hrrr_variable)

    # ---- assemble DataTree ------------------------------------------------
    # Synoptic point data needs xvec CF encoding for Zarr compatibility
    da_synoptic = da_synoptic.swap_dims({"station": "geometry"})
    da_synoptic = da_synoptic.xvec.set_geom_indexes("geometry", crs="EPSG:4326")
    ds_synoptic = da_synoptic.to_dataset().xvec.encode_cf()

    ds_prism = da_prism.to_dataset(name=da_prism.name or "precipitation")
    ds_hrrr = da_hrrr.to_dataset(name=da_hrrr.name or "data")

    dt = xr.DataTree.from_dict(
        {
            "stations/synoptic": ds_synoptic,
            "rasters/prism": ds_prism,
            "rasters/hrrr": ds_hrrr,
        }
    )

    # ---- save -------------------------------------------------------------
    io.save_datatree(dt, output_path)
    print("DONE!")
