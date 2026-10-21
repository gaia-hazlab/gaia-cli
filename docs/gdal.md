# Reading stores with GDAL

GDAL's Zarr driver reads staged and aligned stores (GDAL ≥ 3.4; tested with 3.13), and
so does anything built on it: QGIS, rasterio, the `gdal` command line. Each array is a
subdataset:

```bash
gdal raster info pioneer-10m/pioneer-10m.zarr     # lists the subdatasets
gdal raster info 'ZARR:"pioneer-10m/pioneer-10m.zarr":/dem/usgs10m/elevation'
```

Remote stores work through `/vsicurl/` (or `/vsis3/`), using the store's consolidated
metadata:

```bash
gdal raster info 'ZARR:"/vsicurl/https://s3.kopah.uw.edu/gaia/scratch/public/test/test.zarr":/rasters/prism/precipitation'
```

Each raster variable carries its CRS as a `_CRS` attribute and its nodata as the Zarr
`fill_value`, which is what GDAL reads. Time steps become bands, with the time in each
band's `DIM_time_VALUE` metadata.

## Model inputs from an aligned store

After `gaia align --crs utm --res 10 --aoi $AOI` every raster is on the same grid, so
writing a model's files is a format conversion. For example, a float32 ESRI ASCII grid
with -9999 for nodata:

```bash
gdal raster pipeline read 'ZARR:"pioneer-10m/pioneer-10m.zarr":/soil/solus100/cec7' \
  ! unscale --ot Float32 \
  ! write cec7.vrt
gdal raster calc -i A=cec7.vrt --calc A --ot Float32 --nodata -9999 --propagate-nodata cec7.tif
gdal raster convert cec7.tif cec7.asc
```

`unscale` applies the `scale_factor` of packed integers (soil properties), or only
casts otherwise, but keeps the old nodata value (65535 here); `calc` writes every
nodata cell as -9999. [`scripts/landlab-inputs.sh`](https://github.com/gaia-hazlab/gaia-cli/blob/main/scripts/landlab-inputs.sh)
does this for every layer the landlab workflow reads; the
[Pioneer Fire example](../examples/pioneer-landlab.ipynb) runs it.

Some things to know:

- Use `gaia align` rather than `gdal raster reproject` to put datasets on one grid. Both
  use GDAL's warper, and `align` uses its exact transformer (`--et 0`); GDAL's default
  approximation can move 10 m DEM values by up to metres.
- `gdal raster clip` doesn't mask "cells whose center is inside the AOI": a rectangle in
  EPSG:4326 crops to its bounding box. `gaia align --aoi` does.
- Commands that compute statistics write `.aux.xml` files into the store. Set
  `GDAL_PAM_ENABLED=NO` to avoid that.
- Reading a time series band by band refetches the whole chunk for every band, and
  gaia stores small arrays as one chunk. `gdal mdim convert` reads it once.
