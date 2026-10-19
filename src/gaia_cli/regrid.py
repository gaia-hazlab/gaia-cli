"""Regrid rasters between rectilinear grids in any CRS.

``conservative`` treats every target cell as a polygon and takes the
coverage-weighted mean of the source cells it overlaps, via
``xvec.zonal_stats(method="exactextract")``. NaN source cells are skipped
(xESMF's "conservative_normed"). Other methods delegate to
``rioxarray.reproject_match`` with the rasterio ``Resampling`` of that name.
"""

from __future__ import annotations

from typing import Literal

import geopandas as gpd
import numpy as np
import pyproj
import rioxarray  # noqa: F401  registers the .rio accessor
import shapely
import xarray as xr
import xvec  # noqa: F401  registers the .xvec accessor
from rasterio.enums import Resampling

Method = Literal["conservative", "average", "bilinear", "cubic", "nearest", "mode"]
# A raster variable attribute: the method `gaia align` uses unless --method is given,
# stamped from the dataset module's REGRID_METHOD when it is staged
METHOD_ATTR = "gaia:regrid-method"


def _resolution(da: xr.DataArray) -> tuple[float, float]:
    """(x, y) cell size from coords, falling back to the stored GeoTransform for 1-cell axes."""
    if da.sizes["x"] > 1 and da.sizes["y"] > 1:
        return float(da.x[1] - da.x[0]), float(da.y[1] - da.y[0])
    return da.rio.resolution()


def cell_boxes(da: xr.DataArray) -> gpd.GeoSeries:
    """Cell footprints of *da* (row-major, y then x) in *da*'s CRS."""
    rx, ry = (abs(r) for r in _resolution(da))
    X, Y = np.meshgrid(da.x.values, da.y.values)
    boxes = shapely.box(X - rx / 2, Y - ry / 2, X + rx / 2, Y + ry / 2)
    return gpd.GeoSeries(boxes.ravel(), crs=da.rio.crs)


def cell_areas(da: xr.DataArray) -> np.ndarray:
    """Cell areas of *da* in m², shaped (y, x)."""
    minx, miny, maxx, maxy = da.rio.transform_bounds("EPSG:4326")
    laea = pyproj.CRS.from_proj4(
        f"+proj=laea +lat_0={(miny + maxy) / 2} +lon_0={(minx + maxx) / 2} "
        "+datum=WGS84 +units=m"
    )
    areas = cell_boxes(da).to_crs(laea).area.to_numpy()
    return areas.reshape(da.sizes["y"], da.sizes["x"])


def _spatial_coords(template: xr.DataArray) -> dict:
    """x, y, spatial_ref and any 2D (y, x) coords such as HRRR latitude/longitude."""
    return {
        name: coord
        for name, coord in template.coords.items()
        if (coord.dims and set(coord.dims) <= {"y", "x"}) or name == "spatial_ref"
    }


def regrid_conservative(da: xr.DataArray, template: xr.DataArray) -> xr.DataArray:
    """Area-weighted (first-order conservative) regrid of *da* onto *template*'s grid."""
    other = [d for d in da.dims if d not in ("y", "x")]
    # zonal_stats only needs the raster itself; drop 2D aux coords (e.g. lat/lon)
    source = da.drop_vars(
        [c for c in da.coords if c not in da.dims and set(da[c].dims) & {"y", "x"}]
    )
    cells = cell_boxes(template).to_crs(da.rio.crs)
    stats = source.xvec.zonal_stats(
        cells, x_coords="x", y_coords="y", stats="mean", method="exactextract"
    )
    shape = [da.sizes[d] for d in other] + [template.sizes["y"], template.sizes["x"]]
    data = stats.transpose(*other, "geometry").values.reshape(shape)
    coords = {
        name: coord
        for name, coord in da.coords.items()
        if not ({"y", "x"} & set(coord.dims)) and name != "spatial_ref"
    }
    coords.update(_spatial_coords(template))
    return xr.DataArray(
        data.astype(np.result_type(da.dtype, np.float32)),
        dims=(*other, "y", "x"),
        coords=coords,
        attrs=da.attrs,
        name=da.name,
    )


def regrid(da: xr.DataArray, template: xr.DataArray, method: Method) -> xr.DataArray:
    """Regrid *da* onto *template*'s grid (x, y, CRS) with *method*."""
    # only the template's grid matters; e.g. a scalar time coord would collide with da's
    keep = _spatial_coords(template)
    template = template.drop_vars([c for c in template.coords if c not in keep])
    if method == "conservative":
        return regrid_conservative(da, template)
    da = da.rio.write_nodata(np.nan, encoded=False) if da.dtype.kind == "f" else da
    out = da.rio.reproject_match(template, resampling=Resampling[method])
    # reproject_match's x/y can differ from the template's by float noise
    return out.assign_coords(_spatial_coords(template))
