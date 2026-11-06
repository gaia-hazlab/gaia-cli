"""ERA5 hourly single-level reanalysis (Earthmover), clipped to an AOI.

Needs a free Earthmover account: an API token in ARRAYLAKE_TOKEN, or a cached
`arraylake auth login`. Variables keep their ECMWF short names; temperatures are
converted to degC and water amounts to mm. See dev/era5.md.

Examples
--------
    $ gaia stage era5 -i aoi.geojson -s 2023-09-01 -e 2023-09-30 -o peters-era5
    $ gaia stage era5 -i aoi.geojson -s 2023-09-01 -e 2023-09-30 -o peters-era5 -v tp t2m d2m u10 v10
"""

from __future__ import annotations

import contextlib
import logging
import os
from typing import Literal

import geopandas as gpd
import pandas as pd
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr

from .. import io, provenance

log = logging.getLogger(__name__)

REPO = "earthmover-public/era5"
GROUP = "single/temporal"  # chunked for long time series at few cells
# The repo's bucket config names no endpoint, so the Icechunk S3 client would take one
# from the active AWS profile (e.g. an output store elsewhere); pin AWS.
S3_ENDPOINT = "https://s3.us-east-1.amazonaws.com"
# The data run 1940 to the store's time_coverage_end, which the free tier moves
# forward quarterly (dev/era5.md)
EXTENT = {
    "spatial": {"bbox": [[-180.0, -90.0, 180.0, 90.0]]},
    "temporal": {"interval": [["1940-01-01T00:00:00Z", None]]},
}
KIND = "raster"
# `gaia align` default: area-weighted, so precipitation totals are preserved
REGRID_METHOD = "conservative"

# The single-level variables in GROUP
Variable = Literal[
    "blh", "cape", "cp", "d2m", "fdir", "fg10", "fsr", "hcc", "ie", "lcc", "lsp",
    "mcc", "msl", "sd", "sf", "skt", "slhf", "sp", "ssr", "ssrd", "sst", "stl1",
    "stl2", "stl3", "stl4", "swvl1", "swvl2", "t2m", "tcc", "tcw", "tcwv", "tisr",
    "tp", "tsr", "u10", "u100", "v10", "v100", "zust",
]  # fmt: skip
# Water amounts in metres of water equivalent, staged in mm
TO_MM = {"tp", "cp", "lsp", "sf", "sd"}
# Static fields stored as coordinates of every variable; not staged
STATIC = ("lsm", "z_sfc", "sdor", "slor")
# GRIB_stepType -> how an hour-ending value combines over time
CELL_METHODS = {"accum": "time: sum", "max": "time: maximum", "instant": "time: point"}


@contextlib.contextmanager
def _aws_endpoint(url: str = S3_ENDPOINT):
    """Set ``AWS_ENDPOINT_URL_S3`` to *url* inside the block, then restore it."""
    old = os.environ.get("AWS_ENDPOINT_URL_S3")
    os.environ["AWS_ENDPOINT_URL_S3"] = url
    try:
        yield
    finally:
        if old is None:
            os.environ.pop("AWS_ENDPOINT_URL_S3", None)
        else:
            os.environ["AWS_ENDPOINT_URL_S3"] = old


def open_era5(repo: str = REPO, group: str = GROUP) -> xr.Dataset:
    """Open an ERA5 group from Earthmover's Icechunk repository, lazily.

    Call inside ``_aws_endpoint()``: chunks are read later through the same session.

    Returns
    -------
    xarray.Dataset
        On ``x``/``y`` (longitude −180–180, latitude descending) in EPSG:4326, with
        ``time``; the Icechunk snapshot read in ``attrs["gaia:icechunk-snapshot"]``.
    """
    from arraylake import Client
    from arraylake.exceptions import ArraylakeHttpError, AuthException

    log.info("opening Earthmover Icechunk repo %s (%s)", repo, group)
    try:
        session = Client().get_repo(repo).readonly_session(branch="main")
    except (ArraylakeHttpError, AuthException) as e:
        raise RuntimeError(
            f"could not open {repo} ({e}). ERA5 needs a free Earthmover account: set "
            "ARRAYLAKE_TOKEN to an API key, or run `arraylake auth login`."
        ) from e
    ds = xr.open_zarr(session.store, group=group, chunks=None, zarr_format=3)
    ds = ds.rename(valid_time="time", longitude="x", latitude="y")
    ds = ds.assign_coords(x=(ds.x + 180) % 360 - 180).sortby("x")
    ds = ds.rio.write_crs("EPSG:4326").rio.set_spatial_dims(x_dim="x", y_dim="y")
    ds.attrs["gaia:icechunk-snapshot"] = session.snapshot_id
    return ds


def source_attrs(ds: xr.Dataset) -> dict:
    """Provenance attributes for data read from *ds* (from :func:`open_era5`)."""
    a = ds.attrs
    attrs = {
        "title": "ERA5 hourly reanalysis on single levels",
        "description": a.get("summary", ""),
        "license": "CC-BY-4.0",
        "providers": [
            {
                "name": "ECMWF / Copernicus Climate Change Service",
                "roles": ["producer", "licensor"],
                "url": "https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels",
            },
            {
                "name": "Earthmover",
                "roles": ["processor", "host"],
                "url": "https://earthmover.io",
            },
        ],
        "sci:doi": "10.1002/qj.3803",
        "sci:citation": "Hersbach, H., et al. (2020). The ERA5 global reanalysis. "
        "Quarterly Journal of the Royal Meteorological Society, 146(730), 1999-2049.",
        "acknowledgement": a.get("acknowledgement", ""),
        "links": [
            provenance.link(
                "license",
                "https://creativecommons.org/licenses/by/4.0/",
                type="text/html",
            ),
            provenance.link(
                "via",
                f"https://app.earthmover.io/{REPO}",
                title="Icechunk repository",
            ),
        ],
    }
    if "gaia:icechunk-snapshot" in a:
        attrs["gaia:icechunk-snapshot"] = a["gaia:icechunk-snapshot"]
    return attrs


def _to_staged(da: xr.DataArray) -> xr.DataArray:
    """Convert *da* to staged units (degC, mm), keeping its name."""
    units = da.attrs.get("units")
    keep = {k: v for k, v in da.attrs.items() if not k.startswith("GRIB_")}
    if units == "K":
        da = da - 273.15
        keep.update(units="degC", comment="Converted from K")
    elif da.name in TO_MM:
        da = da * 1000
        keep.update(units="mm", comment=f"Converted from {units} × 1000")
    da.attrs = {**keep, "cell_methods": CELL_METHODS[da.attrs["GRIB_stepType"]]}
    return da


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: tuple[Variable, ...] = ("tp", "sf", "t2m", "swvl1", "swvl2", "sd"),
) -> xr.Dataset:
    """Load ERA5 hourly single-level fields clipped to *aoi* on the 0.25° grid.

    Every value is stamped at the end of its hour: accumulations (``tp``) cover
    ``(time - 1h, time]``, given as CF ``time_bounds``; instantaneous fields
    (``t2m``) are valid at ``time`` (``cell_methods = "time: point"``).

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    variables : tuple of str
        ERA5 short names, e.g. ``tp`` (total precipitation, mm), ``sf``
        (snowfall, mm of water equivalent), ``t2m``
        (2 m temperature, degC), ``swvl1``/``swvl2`` (volumetric soil water,
        0–7 and 7–28 cm), ``sd`` (snow depth, mm of water equivalent).

    Returns
    -------
    xarray.Dataset
        One variable per ERA5 variable plus CF ``time_bounds``, with provenance
        attributes (see :func:`source_attrs`).
    """
    with _aws_endpoint():
        ds = open_era5()
        last = pd.Timestamp(ds.time.values[-1])
        if pd.Timestamp(start_date) > last:
            raise ValueError(
                f"ERA5 in {REPO} ends at {last:%Y-%m-%dT%H:%MZ}; "
                f"the window starts at {start_date}"
            )
        if pd.Timestamp(end_date) > last:
            log.warning("ERA5 ends at %s, before the window does", last)
        ds = ds[list(variables)].sel(time=slice(start_date, end_date))
        ds = ds.drop_vars([c for c in STATIC if c in ds.coords])
        parts = []
        for var in variables:
            da = io.clip_to_aoi(ds[var], aoi)
            log.info(
                "ERA5 %s: reading %d hour(s) x %d x %d cells",
                var,
                da.sizes["time"],
                da.sizes["y"],
                da.sizes["x"],
            )
            parts.append(_to_staged(da.load()))

    # One time axis: hour-ending bounds for all, cell_methods per variable
    first, *rest = parts
    out = io.with_time_bounds(first, -1, 0, first.attrs["cell_methods"])
    for da in rest:
        out[da.name] = da
    out.time.attrs["bounds"] = "time_bounds"  # dropped when adding variables
    out.attrs = {
        **source_attrs(ds),
        "processing:lineage": "clipped to the AOI on the native 0.25° grid; "
        "K → degC, m of water → mm",
        **provenance.time_extent(out),
    }
    return out
