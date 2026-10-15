import logging
import os
import sys

import boto3
import geopandas as gpd
import numpy as np
import obstore
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import xvec  # noqa: F401  registers the .xvec accessor
from obstore.auth.boto3 import Boto3CredentialProvider
from obstore.store import S3Store, from_url
from zarr.storage import ObjectStore

log = logging.getLogger(__name__)


def load_aoi(vector_path: str) -> gpd.GeoDataFrame:
    """Read a vector file into a GeoDataFrame.

    Parameters
    ----------
    vector_path : str
        Path or URL to a vector file (GeoJSON, Shapefile, etc.).

    Returns
    -------
    geopandas.GeoDataFrame
    """
    return gpd.read_file(vector_path)


def clip_to_aoi(
    da: xr.DataArray,
    aoi: gpd.GeoDataFrame,
    *,
    all_touched: bool = True,
) -> xr.DataArray:
    """Clip an xarray DataArray to an area of interest in its native CRS.

    The AOI is reprojected to the raster's CRS; the raster itself is never
    resampled. Reprojection to a common grid is left to ``gaia align``.

    Parameters
    ----------
    da : xarray.DataArray
        Raster data with rioxarray spatial metadata.
    aoi : geopandas.GeoDataFrame
        Area of interest polygon(s).
    all_touched : bool, optional
        Keep every cell the AOI touches (default True), so later area-weighted
        regridding has full coverage along the AOI boundary.

    Returns
    -------
    xarray.DataArray
    """
    aoi_native = aoi.to_crs(da.rio.crs)
    minx, miny, maxx, maxy = aoi_native.total_bounds
    # Pad by one cell so clip_box never drops partially-covered edge cells
    rx, ry = (abs(r) for r in da.rio.resolution())
    da = da.rio.clip_box(minx - rx, miny - ry, maxx + rx, maxy + ry)
    return da.rio.clip(aoi_native.geometry, all_touched=all_touched)


def with_time_bounds(
    da: xr.DataArray,
    start_hours: float,
    end_hours: float,
    cell_methods: str,
) -> xr.Dataset:
    """Return *da* as a Dataset with a CF ``time_bounds`` variable.

    Each value at ``time`` covers ``(time + start_hours, time + end_hours]``.
    ``cell_methods`` (e.g. ``"time: sum"``) says how values combine over time;
    ``gaia align`` uses both to aggregate products onto a common time step.
    """
    ds = da.to_dataset(name=da.name or "data")
    ds[ds.time.name].attrs["bounds"] = "time_bounds"
    start, end = (
        np.timedelta64(round(h * 3600), "s") for h in (start_hours, end_hours)
    )
    bounds = np.stack([ds.time.values + start, ds.time.values + end], axis=-1)
    ds["time_bounds"] = (("time", "nv"), bounds)
    for var in ds.data_vars.values():
        if "time" in var.dims and var.name != "time_bounds":
            var.attrs["cell_methods"] = cell_methods
    return ds


def _s3_store_for_path(output_path: str):
    # obstore doesn't read endpoint_url from ~/.aws/config, so resolve it
    # via boto3 to support S3-compatible stores (e.g. Ceph) set in AWS_PROFILE.
    endpoint = boto3.Session().client("s3").meta.endpoint_url
    return from_url(
        output_path,
        credential_provider=Boto3CredentialProvider(),
        region="us-west-2",
        endpoint=endpoint,
        virtual_hosted_style_request=False,
    )


def _zarr_store_for_path(output_path: str):
    """Return a Zarr-compatible store for *output_path*.

    Handles both local paths and ``s3://`` URLs.
    """
    if output_path.startswith("s3://"):
        return ObjectStore(_s3_store_for_path(output_path))
    return output_path


def output_exists(output_path: str) -> bool:
    """Return True if anything already exists at *output_path* (local or ``s3://``)."""
    if output_path.startswith("s3://"):
        store = _s3_store_for_path(output_path)
        first_chunk = next(iter(obstore.list(store, chunk_size=1)), [])
        return len(first_chunk) > 0
    return os.path.exists(output_path)


def check_output(output_path: str, overwrite: bool = False):
    """Exit early if S3 credentials are missing or *output_path* already exists.

    Pass ``overwrite=True`` to allow replacing an existing store.
    """
    check_s3_env(output_path)
    if not overwrite and output_exists(output_path):
        sys.exit(
            f"Error: output path '{output_path}' already exists. "
            "Use --overwrite to replace it."
        )


def check_s3_env(output_path: str):
    """Raise early if S3 credentials are missing."""
    if output_path.startswith("s3://") and "AWS_PROFILE" not in os.environ:
        raise ValueError(
            "AWS_PROFILE environment variable not set. "
            "Please set it to the AWS profile name with write permissions "
            "to the target S3 bucket."
        )


_CHUNK_TARGET_BYTES = 10 * 1024 * 1024  # 10 MB
# Every save passes write_empty_chunks=True: zarr otherwise skips chunks that
# hold only the fill value (e.g. spatial_ref = 0, all-NaN chunks), and Kopah
# answers anonymous GETs for missing keys with 403 rather than 404, which
# readers treat as an error instead of as fill value.


def _prepare_for_zarr(ds: xr.Dataset) -> tuple[xr.Dataset, dict]:
    """Load small datasets into memory; re-chunk large ones to ~10 MB chunks.

    Returns the (possibly loaded) dataset and an encoding dict.  When the
    dataset is small enough to fit in a single chunk we set the Zarr chunk
    shape equal to the full variable shape so that ``to_zarr`` doesn't
    re-introduce its own default chunking.
    """
    # Fixed-width numpy strings (<U*) have no Zarr v3 spec; store as variable-length UTF-8
    ds = ds.assign_coords(
        {k: v.astype(object) for k, v in ds.coords.items() if v.dtype.kind == "U"}
    )
    ds = ds.assign(
        {k: v.astype(object) for k, v in ds.data_vars.items() if v.dtype.kind == "U"}
    )
    # Store the CF grid_mapping link in attrs, not encoding: xarray never writes
    # an encoding grid_mapping, yet uses it to drop spatial_ref from the
    # "coordinates" attr, so spatial_ref would reopen as a data variable.
    if "spatial_ref" in ds.coords:
        for var in ds.data_vars.values():
            if {"x", "y"} <= set(var.dims):
                var.encoding.pop("grid_mapping", None)
                var.attrs["grid_mapping"] = "spatial_ref"
    if ds.nbytes <= _CHUNK_TARGET_BYTES:
        ds = ds.load()
        # Tell to_zarr to write each variable as a single chunk
        encoding = {
            name: {"chunks": var.shape}
            for name, var in ds.data_vars.items()
            if var.dims
        }
        return ds, encoding

    return ds, {}


def save_zarr(ds: xr.Dataset | xr.DataArray, output_path: str, overwrite: bool = False):
    """Write an xarray Dataset or DataArray to a Zarr store.

    Datasets smaller than 10 MB are loaded into memory (no chunking on
    disk).  Larger datasets are written with dask chunks as-is.

    Parameters
    ----------
    ds : xarray.Dataset or xarray.DataArray
        The gridded (raster) data to save.
    output_path : str
        Filesystem path or ``s3://`` URL.
    overwrite : bool, optional
        If True, replace any existing store at *output_path*. Default False.
    """
    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset(name=ds.name or "data")

    ds, encoding = _prepare_for_zarr(ds)

    store = _zarr_store_for_path(output_path)
    log.info("Writing Zarr to %s", output_path)
    ds.to_zarr(
        store,
        mode="w" if overwrite else "w-",
        zarr_format=3,
        consolidated=False,
        write_empty_chunks=True,
        encoding=encoding,
    )
    log.info("Saved %s", output_path)


def save_xvec_zarr(DA: xr.DataArray, output_path: str, overwrite: bool = False):
    """Save a station DataArray as an xvec-encoded Zarr store.

    Parameters
    ----------
    DA : xarray.DataArray
        Point-geometry data with a ``station`` / ``geometry`` dimension.
    output_path : str
        Filesystem path or ``s3://`` URL.
    overwrite : bool, optional
        If True, replace any existing store at *output_path*. Default False.
    """
    store = _zarr_store_for_path(output_path)
    log.info("Writing Zarr to %s", output_path)

    # TODO: consider different naming conventions for geometry coordinate and station dimension
    DA = DA.swap_dims({"station": "geometry"})
    DA = DA.xvec.set_geom_indexes("geometry", crs="EPSG:4326")

    encoded = DA.to_dataset().xvec.encode_cf()
    encoded, encoding = _prepare_for_zarr(encoded)
    encoded.to_zarr(
        store,
        mode="w" if overwrite else "w-",
        zarr_format=3,
        consolidated=False,
        write_empty_chunks=True,
        encoding=encoding,
    )
    log.info("Saved %s (xvec-encoded)", output_path)


def save_datatree(dt: xr.DataTree, output_path: str, overwrite: bool = False):
    """Write an :class:`xarray.DataTree` to a Zarr store.

    Each leaf node in the tree that is small enough (≤ 10 MB) is loaded
    into memory so that it is stored as a single chunk.

    Parameters
    ----------
    dt : xarray.DataTree
        Hierarchical dataset.
    output_path : str
        Filesystem path or ``s3://`` URL.
    overwrite : bool, optional
        If True, replace any existing store at *output_path*. Default False.
    """
    # Walk the tree and prepare each leaf node
    nodes: dict[str, xr.Dataset] = {}
    encoding: dict[str, dict] = {}
    for path, node in dt.subtree_with_keys:
        ds = node.dataset
        if ds.sizes:
            ds, enc = _prepare_for_zarr(ds)
            if enc:
                # DataTree.to_zarr expects /-prefixed group paths
                key = f"/{path}" if path != "." else "/"
                encoding[key] = enc
        nodes[path] = ds

    prepared = xr.DataTree.from_dict(nodes)
    store = _zarr_store_for_path(output_path)
    log.info("Writing Zarr to %s", output_path)
    prepared.to_zarr(
        store,
        mode="w" if overwrite else "w-",
        zarr_format=3,
        consolidated=False,
        write_empty_chunks=True,
        encoding=encoding or None,
    )
    log.info("Saved %s", output_path)


def open_datatree(input_path: str) -> xr.DataTree:
    """Open a staged Zarr store (local path or ``s3://`` URL) as a DataTree."""
    store = _zarr_store_for_path(input_path)
    return xr.open_datatree(store, engine="zarr", zarr_format=3, consolidated=False)


def open_s3_zarr(
    bucket: str,
    prefix: str,
    *,
    region: str = "us-west-2",
    anonymous: bool = True,
) -> xr.Dataset:
    """Open a Zarr store on S3 as an xarray Dataset.

    Parameters
    ----------
    bucket : str
        S3 bucket name.
    prefix : str
        Key prefix (path within the bucket) to the Zarr store.
    region : str, optional
        AWS region. Default ``us-west-2``.
    anonymous : bool, optional
        If True, access without credentials. Default True.

    Returns
    -------
    xarray.Dataset
    """
    s3 = S3Store(
        bucket=bucket,
        prefix=prefix,
        region=region,
        skip_signature=anonymous,
    )
    store = ObjectStore(s3, read_only=True)
    # xarray's stub doesn't include zarr v3 stores in filename_or_obj, but the zarr engine accepts them at runtime
    return xr.open_dataset(store, engine="zarr", consolidated=False)  # ty: ignore[invalid-argument-type]
