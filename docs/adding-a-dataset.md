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
3. **Clip to the AOI** using `io.clip_to_aoi(da, aoi)`. This reprojects the
   AOI into the raster's CRS and keeps every touched cell; the raster stays
   on its **native grid**. Don't reproject in `stage`: putting products on a
   common grid is `gaia align`'s job (see [reprojecting.md](reprojecting.md)).
4. **Describe the time step.** If each value covers a period rather than an
   instant, return `io.with_time_bounds(da, start_hours, end_hours,
   cell_methods)`, e.g. `(-1, 0, "time: sum")` for hour-ending
   accumulations. `gaia align` relies on this to aggregate in time. Document
   the source's conventions in `docs/<dataset>.md` (see
   [prism.md](prism.md), [hrrr.md](hrrr.md)).
5. **Record provenance.** Set the returned Dataset's `attrs` to the
   source's STAC fields: at least `title`, `license`, `providers` and
   `links` (a `license` link plus `via` or `derived_from` to the source),
   with `version` if the source has one and
   `provenance.time_extent(ds)`. If the source publishes a STAC collection,
   start from `provenance.from_stac_collection(url)` (see `hrrr.py`);
   otherwise write them out (see `nlcd.py`). `io.save_*` adds the run
   attributes (`created`, `gaia:cli-version`, `gaia:command`); see
   [provenance-metadata.md](provenance-metadata.md).
6. **Implement `load()`** returning the clipped `xr.Dataset`, and
   **`stage()`** wiring it to the CLI (`cyclopts.Parameter` annotated
   args for `--input/-i`, `--start/-s`, `--end/-e`, `--output/-o`, plus any
   dataset-specific flags such as `--catalog` or `--variable`). Call
   `io.check_s3_env(output_path)` before loading, and `io.save_zarr(...)` to
   write the result.
7. Register the command in `cli.py`:
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
   with the same CLI parameter conventions as gridded datasets. A
   DataArray's `attrs` describe the variable, so provide the source's
   provenance through a separate `source_attrs(da)` function (see
   `synoptic.py`). Also write a `to_dataset(da)` that swaps the `station`
   dim for `geometry`, sets an Xvec geometry index (`crs="EPSG:4326"`),
   CF-encodes it (`xvec.encode_cf()`) so geometries round-trip through
   Zarr, and sets `source_attrs(da)` (see `synoptic.to_dataset`). `stage()`
   then writes `io.save_zarr(to_dataset(DA), output_path)`.
4. Register the command in `cli.py`, same as above.

## Wiring into stage files

`gaia stage --config` ([stage-file.md](stage-file.md)) runs every dataset
listed in `stagefile.DATASETS`, writing one group per dataset. To add yours:

1. Add `"<name>": (newdataset.load, "rasters")` to `DATASETS` in
   `stagefile.py`, with `"stations"` instead of `"rasters"` for vector
   timeseries. The group is `<prefix>/<name>`.
2. Stage-file options come from `load()`'s signature and type hints, so
   every option needs a hint (`int`, `str`, an `Enum`, or a list of `str`)
   and a default. A `load()` without `start_date`/`end_date` is static and
   takes no window.
3. If `load()` returns something other than a finished `xr.Dataset` with
   provenance attrs, convert it in `stagefile._as_dataset` (vector data:
   your module's `to_dataset`).

Run `pixi run test` after wiring things up, and add a case to
`tests/test_online.py`.
