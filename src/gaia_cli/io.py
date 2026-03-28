import os

import geopandas as gpd
import xarray as xr
from obstore.auth.boto3 import Boto3CredentialProvider
from obstore.store import S3Store, from_url
from zarr.storage import ObjectStore


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
    reproject: bool = False,
) -> xr.DataArray:
    """Clip an xarray DataArray to an area of interest.

    For data in a non-EPSG:4326 CRS (e.g. HRRR Lambert Conformal), set
    ``reproject=True`` to first clip to a generous bounding box, reproject
    to EPSG:4326, and then clip to the exact AOI polygon.

    Parameters
    ----------
    da : xarray.DataArray
        Raster data with rioxarray spatial metadata.
    aoi : geopandas.GeoDataFrame
        Area of interest polygon(s).
    reproject : bool, optional
        If True, clip_box → reproject to EPSG:4326 → clip to polygon.
        Default is False (direct clip).

    Returns
    -------
    xarray.DataArray
    """
    if reproject:
        bounds = aoi.total_bounds  # minx, miny, maxx, maxy
        # Generous bounding box to avoid edge artefacts after reprojection
        pad = 5
        da = da.rio.clip_box(
            minx=bounds[0] - pad,
            miny=bounds[1] - pad,
            maxx=bounds[2] + pad,
            maxy=bounds[3] + pad,
            crs="EPSG:4326",
        )
        da = da.rio.reproject("EPSG:4326")

    return da.rio.clip(aoi.geometry)


def _zarr_store_for_path(output_path: str):
    """Return a Zarr-compatible store for *output_path*.

    Handles both local paths and ``s3://`` URLs.
    """
    if output_path.startswith("s3://"):
        s3writable = from_url(
            output_path,
            credential_provider=Boto3CredentialProvider(),
            region="us-west-2",
        )
        return ObjectStore(s3writable)
    return output_path


def check_s3_env(output_path: str):
    """Raise early if S3 credentials are missing."""
    if output_path.startswith("s3://") and "AWS_PROFILE" not in os.environ:
        raise ValueError(
            "AWS_PROFILE environment variable not set. "
            "Please set it to the AWS profile name with write permissions "
            "to the target S3 bucket."
        )


_CHUNK_TARGET_BYTES = 10 * 1024 * 1024  # 10 MB


def _prepare_for_zarr(ds: xr.Dataset) -> tuple[xr.Dataset, dict]:
    """Load small datasets into memory; re-chunk large ones to ~10 MB chunks.

    Returns the (possibly loaded) dataset and an encoding dict.  When the
    dataset is small enough to fit in a single chunk we set the Zarr chunk
    shape equal to the full variable shape so that ``to_zarr`` doesn't
    re-introduce its own default chunking.
    """
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


def save_zarr(ds: xr.Dataset | xr.DataArray, output_path: str):
    """Write an xarray Dataset or DataArray to a Zarr store.

    Datasets smaller than 10 MB are loaded into memory (no chunking on
    disk).  Larger datasets are written with dask chunks as-is.

    Parameters
    ----------
    ds : xarray.Dataset or xarray.DataArray
        The gridded (raster) data to save.
    output_path : str
        Filesystem path or ``s3://`` URL.
    """
    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset(name=ds.name or "data")

    ds, encoding = _prepare_for_zarr(ds)

    store = _zarr_store_for_path(output_path)
    ds.to_zarr(store, mode="w-", zarr_format=3, consolidated=False, encoding=encoding)
    print(f"Data saved to {output_path}")


def save_xvec_zarr(DA: xr.DataArray, output_path: str):
    """Save a station DataArray as an xvec-encoded Zarr store.

    Parameters
    ----------
    DA : xarray.DataArray
        Point-geometry data with a ``station`` / ``geometry`` dimension.
    output_path : str
        Filesystem path or ``s3://`` URL.
    """
    store = _zarr_store_for_path(output_path)

    # TODO: consider different naming conventions for geometry coordinate and station dimension
    DA = DA.swap_dims({"station": "geometry"})
    DA = DA.xvec.set_geom_indexes("geometry", crs="EPSG:4326")

    encoded = DA.to_dataset().xvec.encode_cf()
    encoded, encoding = _prepare_for_zarr(encoded)
    encoded.to_zarr(store, mode="w-", zarr_format=3, consolidated=False, encoding=encoding)
    print(f"Data saved to {output_path} as Xvec-encoded Zarr.")


def save_datatree(dt: xr.DataTree, output_path: str):
    """Write an :class:`xarray.DataTree` to a Zarr store.

    Each leaf node in the tree that is small enough (≤ 10 MB) is loaded
    into memory so that it is stored as a single chunk.

    Parameters
    ----------
    dt : xarray.DataTree
        Hierarchical dataset.
    output_path : str
        Filesystem path or ``s3://`` URL.
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
    prepared.to_zarr(store, mode="w-", zarr_format=3, consolidated=False, encoding=encoding or None)
    print(f"DataTree saved to {output_path}")


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
    return xr.open_dataset(store, engine="zarr", consolidated=False)
