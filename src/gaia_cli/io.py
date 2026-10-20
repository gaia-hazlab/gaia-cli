import json
import logging
import os
import posixpath
import sys
import warnings
from pathlib import Path

import boto3
import geopandas as gpd
import numpy as np
import obstore
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import xvec  # noqa: F401  registers the .xvec accessor
import zarr
from obstore.auth.boto3 import Boto3CredentialProvider
from obstore.store import LocalStore, from_url
from zarr.storage import ObjectStore

from . import provenance

log = logging.getLogger(__name__)


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
    # B2 also requires the profile's region (e.g. us-west-004) in signatures.
    session = boto3.Session()
    endpoint = session.client("s3").meta.endpoint_url
    return from_url(
        output_path,
        credential_provider=Boto3CredentialProvider(session),
        region=session.region_name or "us-west-2",
        endpoint=endpoint,
        virtual_hosted_style_request=False,
    )


class _UnconditionalObjectStore(ObjectStore):
    """ObjectStore that avoids conditional PUTs.

    zarr creates new keys with ``If-None-Match: *``, which Backblaze B2 rejects
    (501 NotImplemented). Check-then-write is fine here because outputs have a
    single writer and ``check_output`` refuses existing paths up front.
    """

    async def set_if_not_exists(self, key, value) -> None:
        if not await self.exists(key):
            await self.set(key, value)


def _zarr_store_for_path(output_path: str):
    """Return a Zarr-compatible store for *output_path*.

    Handles both local paths and ``s3://`` URLs.
    """
    if output_path.startswith("s3://"):
        return _UnconditionalObjectStore(_s3_store_for_path(output_path))
    return output_path


def object_store(path: str, mkdir: bool = False):
    """An obstore store rooted at *path* (local directory or ``s3://`` URL)."""
    if path.startswith("s3://"):
        return _s3_store_for_path(path)
    return LocalStore(Path(path).resolve(), mkdir=mkdir)


def write_text(path: str, text: str):
    """Write *text* to *path* (local or ``s3://``), creating local directories."""
    parent, name = posixpath.split(path)
    store = object_store(parent or ".", mkdir=True)
    obstore.put(store, name, text.encode())


def read_text(path: str) -> str | None:
    """The text at *path* (local or ``s3://``), or None if there is none."""
    parent, name = posixpath.split(path)
    try:
        return bytes(obstore.get(object_store(parent or "."), name).bytes()).decode()
    except FileNotFoundError:
        return None


def write_json(path: str, doc: dict):
    """Write *doc* as JSON to *path* (local or ``s3://``)."""
    write_text(path, json.dumps(doc, indent=1))


def store_path(output: str, name: str) -> str:
    """Path of the Zarr store *name* in the output folder *output*."""
    return f"{output.rstrip('/')}/{name}.zarr"


def folder_name(output: str) -> str:
    """The basename of the output folder *output* (local or ``s3://``)."""
    if "://" in output:
        return posixpath.basename(output.rstrip("/"))
    return Path(output).resolve().name


def item_path(output_path: str) -> str:
    """Path of the STAC Item written next to the store at *output_path*."""
    return output_path.rstrip("/").removesuffix(".zarr") + ".json"


def _write_item(output_path: str, tree: xr.DataTree, assets: dict | None = None):
    store_name = posixpath.basename(output_path.rstrip("/"))
    write_json(item_path(output_path), provenance.stac_item(tree, store_name, assets))
    log.info("Wrote STAC item %s", item_path(output_path))


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
    """Exit early if *output_path* is on S3 and no AWS credentials are set."""
    if output_path.startswith("s3://") and boto3.Session().get_credentials() is None:
        sys.exit(
            f"Error: no AWS credentials to write to {output_path}. Set AWS_PROFILE "
            "to a profile with write access, or AWS_ACCESS_KEY_ID and "
            "AWS_SECRET_ACCESS_KEY."
        )


_CHUNK_TARGET_BYTES = 10 * 1024 * 1024  # 10 MB
_PACKING = ("dtype", "scale_factor", "add_offset", "_FillValue", "fill_value")


def _packing(var: xr.DataArray) -> dict:
    """*var*'s CF packing, if its encoding stores it as another dtype (e.g. soil's
    integers with a scale_factor), else nothing.

    An explicit encoding replaces the variable's own, so the single-chunk
    encoding must repeat it. Variables written as they are keep NaN fills.
    """
    if np.dtype(var.encoding.get("dtype", var.dtype)) == var.dtype:
        return {}
    return {k: var.encoding[k] for k in _PACKING if k in var.encoding}


# Every save consolidates metadata so anonymous readers never need to LIST the
# store: Backblaze B2 public buckets reject unauthenticated list requests.
warnings.filterwarnings("ignore", message="Consolidated metadata is currently not part")
# Every save passes write_empty_chunks=True: zarr otherwise skips chunks that
# hold only the fill value (e.g. spatial_ref = 0, all-NaN chunks), and Kopah
# answers anonymous GETs for missing keys with 403 rather than 404, which
# readers treat as an error instead of as fill value.


def _prepare_for_zarr(ds: xr.Dataset) -> tuple[xr.Dataset, dict]:
    """Load datasets up to ~10 MB into memory and store each variable as one chunk.

    Returns the (possibly loaded) dataset and an encoding dict. A small
    dataset gets Zarr chunks equal to each variable's full shape, so
    ``to_zarr`` doesn't introduce its own default chunking. A larger one is
    returned unchanged: dask arrays keep their chunks (which Zarr requires to
    be uniform, apart from the last), and in-memory arrays get Zarr's default
    chunking.
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
    # GDAL's Zarr driver reads neither: it takes the CRS from a `_CRS` attribute
    # and nodata from the array's fill_value, which xarray leaves at 0 for
    # integers, writing _FillValue only as an attribute (docs/gdal.md)
    if "spatial_ref" in ds.coords:
        wkt = ds.spatial_ref.attrs.get("crs_wkt")
        for var in ds.data_vars.values():
            if {"x", "y"} <= set(var.dims):
                var.encoding.pop("grid_mapping", None)
                var.attrs["grid_mapping"] = "spatial_ref"
                if wkt:
                    var.attrs["_CRS"] = {"wkt": wkt}
    for var in ds.data_vars.values():
        dtype = np.dtype(var.encoding.get("dtype", var.dtype))
        if dtype.kind in "iu" and "_FillValue" in var.encoding:
            var.encoding["fill_value"] = var.encoding["_FillValue"]
    if ds.nbytes <= _CHUNK_TARGET_BYTES:
        ds = ds.load()
        # Tell to_zarr to write each variable as a single chunk
        encoding = {
            name: {**_packing(var), "chunks": var.shape}
            for name, var in ds.data_vars.items()
            if var.dims
        }
        return ds, encoding

    return ds, {}


def save_datatree(
    dt: xr.DataTree,
    output_path: str,
    overwrite: bool = False,
    append: bool = False,
    assets: dict | None = None,
):
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
    append : bool, optional
        If True, add *dt*'s groups to the existing store at *output_path*,
        keeping its other groups; the root attributes and STAC Item then
        describe the whole store. Default False.
    assets : dict, optional
        Extra assets for the STAC Item, e.g. ``{"config": ...}``.
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
    prepared.attrs = provenance.stamp_tree(prepared)
    store = _zarr_store_for_path(output_path)
    log.info("Writing Zarr to %s", output_path)
    prepared.to_zarr(
        store,
        mode="a" if append else "w" if overwrite else "w-",
        zarr_format=3,
        consolidated=True,
        write_empty_chunks=True,
        encoding=encoding or None,
    )
    if append:
        # The root attributes written above describe only the new groups
        prepared = open_datatree(output_path)
        prepared.attrs = provenance.stamp_tree(prepared)
        root = zarr.open_group(store, mode="r+", zarr_format=3)
        root.attrs.put({str(k): v for k, v in prepared.attrs.items()})
        zarr.consolidate_metadata(store)
    log.info("Saved %s", output_path)
    _write_item(output_path, prepared, assets)


def groups(path: str) -> set[str]:
    """Paths of the groups in the store at *path*, without a leading ``/``."""
    return {g.lstrip("/") for g in open_datatree(path).groups if g != "/"}


def open_datatree(input_path: str) -> xr.DataTree:
    """Open a staged Zarr store (local path or ``s3://`` URL) as a DataTree."""
    store = _zarr_store_for_path(input_path)
    return xr.open_datatree(store, engine="zarr", zarr_format=3, consolidated=None)
