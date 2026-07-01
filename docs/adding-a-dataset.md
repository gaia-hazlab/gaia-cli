# Adding a New Dataset to `gaia stage`

Each dataset lives in its own module under `src/gaia_cli/` (e.g. `prism.py`,
`hrrr.py`, `synoptic.py`) and exposes two functions:

- `load(aoi, start_date, end_date, ...) -> xr.DataArray` — fetch + clip the
  data to the AOI, no I/O.
- `stage(vectorPath, start_date, end_date, output_path, ...)` — the CLI
  entrypoint: loads the AOI, calls `load`, and writes a Zarr store via
  `io.py`.

Datasets fall into two shapes: **gridded** (raster) data that goes
`xarray -> zarr`, and **vector timeseries** (point) data that goes
`xvec -> zarr`. Pick the section below that matches your new dataset.

## 1. Gridded (raster) datasets — xarray → Zarr

Examples: `prism.py`, `hrrr.py`.

1. **Open the source** into an `xr.Dataset`/`DataArray` with dimensions that
   include a spatial CRS readable by `rioxarray` (`da.rio.crs`). Use
   whatever access pattern fits the source (STAC + `odc.stac`, Icechunk,
   plain `xr.open_dataset`, etc.) — see `open_prism` / `open_hrrr` for
   examples.
2. **Standardize naming** to match the conventions in the main
   [README](../README.md#design-ideas):
   - Spatial dims: `x`, `y`
   - Time dim: `time` (rename to `<dataset>_time` only if it would collide
     with another time dimension in a combined tree)
   - Variable name + units, e.g. `precipitation` in `mm`
3. **Clip to the AOI** using `io.clip_to_aoi(da, aoi, reproject=...)`. Set
   `reproject=True` if the source CRS isn't already EPSG:4326 (e.g. HRRR's
   Lambert Conformal grid) — this clips to a padded bounding box first,
   reprojects, then clips to the exact polygon.
4. **Implement `load()`** returning the clipped `xr.DataArray`, and
   **`stage()`** wiring it to the CLI (`cyclopts.Parameter` annotated
   args for `--input/-i`, `--start/-s`, `--end/-e`, `--output/-o`, plus any
   dataset-specific flags such as `--catalog` or `--variable`). Call
   `io.check_s3_env(output_path)` before loading, and `io.save_zarr(...)` to
   write the result.
5. Register the command in `cli.py`:
   ```python
   from .newdataset import stage as stage_newdataset
   stage_app.command(stage_newdataset, name="newdataset")
   ```

See `prism.py` and `hrrr.py` for complete reference implementations.

## 2. Vector timeseries datasets — xvec → Zarr

Example: `synoptic.py`.

1. **Fetch station/point data** from the source API and build a tidy
   per-station `xr.DataArray` indexed by `time`, with `station`,
   `elevation`, and a `geometry` (shapely `Point`) coordinate per station
   (see `construct_uniform_datarray`).
2. **Concatenate stations** along a `station` dimension
   (`xr.concat(..., dim="station", join="outer")`), so gaps between
   stations with different coverage become `NaN` rather than raising.
3. **Implement `load()`** returning that `xr.DataArray`, and **`stage()`**
   with the same CLI parameter conventions as gridded datasets. Call
   `io.save_xvec_zarr(DA, output_path)` instead of `save_zarr` — this
   swaps the `station` dim for `geometry`, sets an Xvec geometry index
   (`crs="EPSG:4326"`), and CF-encodes it (`xvec.encode_cf()`) so geometries
   round-trip through Zarr.
4. Register the command in `cli.py`, same as above.

## Wiring into `gaia stage all`

`all.py` combines every dataset into a single `xr.DataTree` and writes one
Zarr store. To make a new dataset compatible with `stage all`:

1. Import your module in `all.py`: `from . import newdataset as _newdataset`.
2. Call `_newdataset.load(aoi, start_date, end_date, ...)` alongside the
   existing `_synoptic.load` / `_prism.load` / `_hrrr.load` calls. Add any
   dataset-specific keyword args (like `hrrr_variable`) as extra parameters
   on `stage_all` if you need to expose them.
3. Convert the result to a named `xr.Dataset`:
   - Gridded: `ds = da.to_dataset(name=da.name or "<varname>")`
   - Vector: apply the same xvec encoding steps as `save_xvec_zarr`
     (`swap_dims({"station": "geometry"})` →
     `xvec.set_geom_indexes("geometry", crs="EPSG:4326")` →
     `to_dataset().xvec.encode_cf()`)
4. Add it to the `xr.DataTree.from_dict({...})` call under an appropriate
   group path — `rasters/<name>` for gridded data, `stations/<name>` for
   vector timeseries.
5. `io.save_datatree(dt, output_path)` handles the rest (per-node chunking
   and Zarr write) automatically — no changes needed there.

Run `pixi run test` after wiring things up; `tests/test_all.py` exercises
the combined `stage all` path.
