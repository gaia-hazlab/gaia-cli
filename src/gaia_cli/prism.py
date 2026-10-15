from __future__ import annotations

import logging
import sys
from typing import Annotated

import geopandas as gpd
import numpy as np
import odc.stac
import xarray as xr
from cyclopts import Parameter

from . import io, prism_cache

log = logging.getLogger(__name__)

# GDAL settings for efficient remote reading
odc.stac.configure_rio(cloud_defaults=True)

CELL_METHODS = {"ppt": "time: sum", "tmin": "time: minimum", "tmax": "time: maximum"}
# Variables renamed on load to match the other gaia products
RENAMES = {"ppt": "precipitation"}


def open_prism(
    variables: list[str] | tuple[str, ...] = ("ppt",),
    start: str | None = None,
    end: str | None = None,
    revalidate: bool = True,
) -> xr.Dataset:
    """Open PRISM daily 800m grids from the gaia S3 cache.

    Grids missing from the cache are fetched from PRISM first, which needs
    ``AWS_PROFILE`` with write access to the cache (see docs/prism.md).

    Parameters
    ----------
    variables : list of str
        PRISM variables, e.g. ``["ppt", "tmean"]``. ``ppt`` is renamed
        ``precipitation``.
    start : str
        Start date (inclusive), e.g. '2025-12-01'.
    end : str, optional
        End date (inclusive). Defaults to *start*.
    revalidate : bool
        Re-fetch provisional grids that PRISM has since re-released.

    Returns
    -------
    xarray.Dataset
    """
    if start is None:
        raise ValueError("start date is required")
    cache = prism_cache.PrismCache.from_env()
    entries = cache.ensure(variables, start, end or start, revalidate=revalidate)
    items = prism_cache.items(entries)
    if not items:
        raise ValueError(f"No PRISM data available for {start} to {end}")

    ds = odc.stac.load(items, bands=list(variables), chunks={})
    ds = ds.rename({"longitude": "x", "latitude": "y"})

    dates = [t.date() for t in ds.time.to_index().to_pydatetime()]
    for var in variables:
        # Use nans instead of -9999
        ds[var] = ds[var].where(ds[var] != ds[var].rio.nodata)
        ds[var].rio.write_nodata(np.nan, inplace=True)
        ds[var].attrs["units"] = prism_cache.UNITS[var]
        # PRISM release number per time step; < 8 means provisional
        ds[var].attrs["prism_release"] = [
            prism_cache.release(entries[var][d]) if d in entries[var] else 0
            for d in dates
        ]

    return ds.rename({v: RENAMES[v] for v in variables if v in RENAMES})


def stage(
    vectorPath: Annotated[str, Parameter(name=["--input", "-i"])],
    start_date: Annotated[str, Parameter(name=["--start", "-s"])],
    end_date: Annotated[str, Parameter(name=["--end", "-e"])],
    output_path: Annotated[str, Parameter(name=["--output", "-o"])],
    variables: Annotated[
        tuple[str, ...], Parameter(name=["--variables", "-v"], consume_multiple=True)
    ] = ("ppt",),
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
):
    """Stage PRISM daily 800m grids clipped to an AOI.

    Data come from the gaia PRISM cache (s3://gaia/public/prism). Grids not
    yet cached are fetched from PRISM, which requires AWS_PROFILE with write
    access to the cache. AWS_PROFILE is also required if writing to S3.

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
    variables : list of str
        PRISM variables (ppt, tmin, tmax, tmean, tdmean, vpdmin, vpdmax).
        Default ppt, written as ``precipitation``.
    overwrite : bool
        Replace the output store if it already exists. Default False.

    Examples
    --------
    Stage PRISM for Skagit Basin December 2025:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage prism -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit_prism.zarr
    """
    io.check_output(output_path, overwrite)

    aoi = io.load_aoi(vectorPath)
    clipped = load(aoi, start_date, end_date, variables)
    io.save_zarr(clipped, output_path, overwrite=overwrite)

    log.info("Done")


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: list[str] | tuple[str, ...] = ("ppt",),
) -> xr.Dataset:
    """Load PRISM daily grids clipped to *aoi* on their native grid.

    A PRISM day is 1200–1200 UTC with a day-ending label, so the value for
    date D covers ``(D-1 12:00, D 12:00]`` UTC (see docs/prism.md).

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    variables : list of str
        PRISM variables. Default ``("ppt",)``, returned as ``precipitation``.

    Returns
    -------
    xarray.Dataset
        One variable per PRISM variable plus CF ``time_bounds``.
    """
    ds = open_prism(variables, start=start_date, end=end_date)
    parts = []
    for var in variables:
        da = io.clip_to_aoi(ds[RENAMES.get(var, var)], aoi)
        log.info(
            "PRISM %s: reading %d day(s) x %d x %d cells from cached COGs",
            var,
            da.sizes["time"],
            da.sizes["y"],
            da.sizes["x"],
        )
        da = da.compute()
        parts.append(
            io.with_time_bounds(
                da,
                start_hours=-12,
                end_hours=12,
                cell_methods=CELL_METHODS.get(var, "time: mean"),
            )
        )
    return xr.merge(parts, combine_attrs="drop_conflicts")


def fill(
    start_date: Annotated[str, Parameter(name=["--start", "-s"])],
    end_date: Annotated[str, Parameter(name=["--end", "-e"])],
    variables: Annotated[
        tuple[str, ...], Parameter(name=["--variables", "-v"], consume_multiple=True)
    ] = ("ppt",),
    revalidate: bool = True,
):
    """Fetch PRISM grids into the gaia cache (s3://gaia/public/prism).

    Downloads grids that are missing and, unless --no-revalidate, re-fetches
    provisional grids that PRISM has since re-released. Requires AWS_PROFILE
    with write access to the cache.

    Parameters
    ----------
    start_date : str
        Start date (inclusive), e.g. '2025-12-01'.
    end_date : str
        End date (inclusive).
    variables : list of str
        PRISM variables. Default ppt.
    revalidate : bool
        Check cached provisional grids for newer releases. Default True.

    Examples
    --------
        $ gaia prism fill -s 2025-12-01 -e 2025-12-31 -v ppt tmean
    """
    cache = prism_cache.PrismCache.from_env()
    if not cache.writable:
        sys.exit(
            "Error: set PRISM_CACHE_PROFILE (or AWS_PROFILE) to a profile that "
            "can write to the cache."
        )
    entries = cache.ensure(variables, start_date, end_date, revalidate=revalidate)
    for var, cached in entries.items():
        provisional = sum(
            prism_cache.release(it) < prism_cache.FINAL_RELEASE
            for it in cached.values()
        )
        log.info("%s: %d grid(s) cached, %d provisional", var, len(cached), provisional)


def rebuild_catalog():
    """Regenerate the cache's STAC catalogs and collection from its items.

    The catalogs are normally updated as grids are cached; use this to
    repair them. Requires AWS_PROFILE with write access to the cache.
    """
    cache = prism_cache.PrismCache.from_env()
    if not cache.writable:
        sys.exit(
            "Error: set PRISM_CACHE_PROFILE (or AWS_PROFILE) to a profile that "
            "can write to the cache."
        )
    n = cache.rebuild_catalog()
    log.info("Rebuilt catalog from %d item(s): %s/collection.json", n, cache.href_base)
