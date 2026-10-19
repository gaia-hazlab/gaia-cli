from __future__ import annotations

import logging
import sys
from typing import Annotated, Literal

import geopandas as gpd
import numpy as np
import odc.stac
import xarray as xr
from cyclopts import Parameter

from . import io, prism_cache, provenance

log = logging.getLogger(__name__)

# GDAL settings for efficient remote reading
odc.stac.configure_rio(cloud_defaults=True)

# STAC Collection extent: the bbox of the gaia cache's collection (prism_cache/stac.py);
# the archive starts in 1981 (https://prism.oregonstate.edu/documents/PRISM_datasets.pdf)
EXTENT = {
    "spatial": {"bbox": [[-125.020833, 24.0625, -66.479167, 49.9375]]},
    "temporal": {"interval": [["1981-01-01T00:00:00Z", None]]},
}
KIND = "raster"
TIMED = True

# The daily PRISM variables gaia stages, of those it publishes (prism_cache.VARIABLES)
Variable = Literal["ppt", "tmin", "tmax", "tmean"]

# How each day's value combines over time; tmean is a daily mean
CELL_METHODS = {
    "ppt": "time: sum",
    "tmin": "time: minimum",
    "tmax": "time: maximum",
}
# CF names; tmin, tmax and tmean differ only in cell_methods
STANDARD_NAMES = {
    "ppt": "lwe_thickness_of_precipitation_amount",
    "tmin": "air_temperature",
    "tmax": "air_temperature",
    "tmean": "air_temperature",
}
LONG_NAMES = {
    "ppt": "total precipitation",
    "tmin": "minimum temperature",
    "tmax": "maximum temperature",
    "tmean": "mean temperature",
}
# Variables renamed on load to match the other gaia products
RENAMES = {"ppt": "precipitation"}


def open_prism(
    variables: list[Variable] | tuple[Variable, ...],
    start: str,
    end: str | None = None,
    revalidate: bool = True,
) -> xr.Dataset:
    """Open PRISM daily 800m grids from the gaia S3 cache.

    Grids missing from the cache are fetched from PRISM first, which needs
    ``AWS_PROFILE`` with write access to the cache (see docs/prism.md).

    Parameters
    ----------
    variables : list of str
        PRISM variables, e.g. ``["ppt", "tmin", "tmax"]``. ``ppt`` is renamed
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
        With provenance attributes from the cache's STAC collection.
    """
    cache = prism_cache.PrismCache.from_env()
    entries = cache.ensure(variables, start, end or start, revalidate=revalidate)
    items = prism_cache.items(entries)
    if not items:
        raise ValueError(f"No PRISM data available for {start} to {end}")

    ds = odc.stac.load(items, bands=list(variables), chunks={})
    ds = ds.rename({"longitude": "x", "latitude": "y"})
    ds.attrs = provenance.collection_attrs(
        cache.read_collection(), f"{cache.href_base}/{prism_cache.COLLECTION_KEY}"
    )

    dates = [t.date() for t in ds.time.to_index().to_pydatetime()]
    for var in variables:
        # Use nans instead of -9999
        ds[var] = ds[var].where(ds[var] != ds[var].rio.nodata)
        ds[var].rio.write_nodata(np.nan, inplace=True)
        ds[var].attrs["units"] = prism_cache.UNITS[var]
        ds[var].attrs["standard_name"] = STANDARD_NAMES[var]
        ds[var].attrs["long_name"] = f"PRISM daily {LONG_NAMES[var]}"
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
        tuple[Variable, ...],
        Parameter(name=["--variables", "-v"], consume_multiple=True),
    ] = ("ppt",),
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
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
        Output folder for the Zarr store (named after the folder), its STAC
        Item and the resolved stage.yaml.
    variables : list of str
        PRISM variables (ppt, tmin, tmax, tmean).
        Default ppt, written as ``precipitation``; the others keep their
        names.
    dry_run : bool
        Validate the flags, read the AOI and print the resolved configuration
        without fetching any data.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.
    skip_existing : bool
        Add to an existing store, unless it already has this dataset's group.

    Examples
    --------
    Stage PRISM for Skagit Basin December 2025:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage prism -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-prism

    With daily minimum and maximum temperature:
        $ gaia stage prism -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit-prism -v ppt tmin tmax
    """
    from .stagefile import stage_dataset  # stagefile imports this module

    stage_dataset(
        "prism",
        vectorPath,
        start_date,
        end_date,
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
        variables=variables,
    )


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: list[Variable] | tuple[Variable, ...] = ("ppt",),
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
        PRISM variables (ppt, tmin, tmax, tmean).
        Default ``("ppt",)``, returned as ``precipitation``.

    Returns
    -------
    xarray.Dataset
        One variable per PRISM variable plus CF ``time_bounds``, with
        provenance attributes.
    """
    # Days are selected by their date label, so timestamps keep only the date
    start_date, end_date = (str(t)[:10] for t in (start_date, end_date))
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
    out = xr.merge(parts, compat="override", combine_attrs="drop_conflicts")
    out.attrs = {
        **ds.attrs,
        "processing:lineage": "clipped to the AOI on the native grid; -9999 → NaN",
        **provenance.time_extent(out),
    }
    return out


def fill(
    start_date: Annotated[str, Parameter(name=["--start", "-s"])],
    end_date: Annotated[str, Parameter(name=["--end", "-e"])],
    variables: Annotated[
        tuple[Variable, ...],
        Parameter(name=["--variables", "-v"], consume_multiple=True),
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
