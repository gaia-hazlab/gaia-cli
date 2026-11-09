# Adding a New Dataset to `gaia stage`

Each dataset lives in its own module under `src/gaia_cli/datasets/` (e.g.
`prism.py`, `hrrr.py`, `synoptic.py`), registered in `datasets.DATASETS`. The
module has:

- a **module docstring**: the `gaia stage <name>` help. Its first line is the
  summary in `gaia stage --help`.
- **`load(aoi, [start_date, end_date,] ...) -> xr.Dataset`**: fetch the data and
  clip it to the AOI, with no I/O. Timed datasets take `start_date` and `end_date`
  (inclusive); static layers take neither. The runner passes the window's end less one
  microsecond as `end_date`. `load()` may instead return an `xr.DataTree`, whose
  children are written under the dataset's group (`seismic`: `sites`, `channels`
  and one per sample rate, `waveforms/seismic/50hz`), with the tree's own attrs on
  the group. Any xvec geometry (stations, waveform sites) is CF-encoded for Zarr
  when the store is written (`io.encode_cf`); readers decode it with `io.decode_cf`.
- **`EXTENT`**: a STAC Collection `extent` (`{"spatial": {"bbox": [[w, s, e, n]]},
  "temporal": {"interval": [[start, end]]}}`, `None` for an open end). Copy it from
  the provider's STAC Collection when there is one, otherwise write it by hand,
  with a comment saying where the values come from. Runs skip datasets whose
  extent the AOI or window misses ([stage-file.md](stage-file.md#coverage)).
- **`KIND`**: `"raster"` (group `rasters/<name>`), `"stations"` (group
  `stations/<name>`) or
  `"waveforms"` (group `waveforms/<name>`, with `sites` and `channels` tables and one
  child per sample rate; see [seismic.md](seismic.md)).

and optionally:

- **`VARIANT`**: the name of a `load()` option (an `Enum` with a default) whose
  values are separate products, each written to its own top-level group
  `<name>/<value>` instead of `<GROUPS[KIND]>/<name>`. Stage files key them as
  `<name>/<value>` (see `dem.py`, [stage-file.md](stage-file.md)).
- **`SCENES = True`**: the dataset chooses among scenes, so a stage file may give it
  `pre`/`post` windows with `select` and `count` ([stage-file.md](stage-file.md#choosing-scenes)).
  `load()` then also takes `windows: list[Window] | None = None`, which the runner sets
  (it is not an option), and `start_date`/`end_date` span the windows exactly. Group
  the items into scenes with `stac.group`, pick them with `stac.choose` before loading
  them, and label each with `stac.label_windows` (see `sentinel1_rtc.py`). The
  runner stamps the group `gaia:scenes`, so `gaia align` keeps its acquisition times
  rather than combining them into the other datasets' time steps.
- **`TOKENS`**: environment variables `load()` needs, each with what it is, e.g.
  `{"SYNOPTIC_TOKEN": "a Synoptic API token, https://docs.synopticdata.com"}`.
  `gaia stage` (and `--dry-run`) warns about any that aren't set before loading
  anything; the dataset still raises its own error when it runs.
- **`REGRID_METHOD`**: `gaia align`'s default spatial method for the dataset's
  raster variables, one of `regrid.Method` (`conservative`, `average`,
  `bilinear`, `cubic`, `nearest`). The runner stores it as each variable's
  `gaia:regrid-method`, and `--method` overrides it. Without it, continuous
  variables use `nearest`, the GDAL and rasterio default. Accumulations and means
  want `conservative` (PRISM, HRRR); smooth fields such as elevation want
  `bilinear`. Categorical variables (`flag_values`) always use `mode`.
- **`PAD`**: a `pd.Timedelta`, how far beyond each end of its window the dataset is
  read. Default `window.PAD`, one day, so that periods straddling the window's edges
  are complete for `gaia align`; `seismic` sets `pd.Timedelta(0)`.
- an **`origin`** argument of `load()` (`dict | None`, default None): the stage
  file's `origin` (`{latitude, longitude}`), None with flags only. Like the window
  arguments it is reserved: the runner fills it, and it is never an option or flag.
  Record in the attrs whether the file's origin or a fallback was used.

## Options come from `load()`

There is no per-dataset CLI code. `runner.dataset_command` builds
`gaia stage <name>` from `load()`: the common flags (`-i`, `-s`/`-e` for timed
datasets, `-o`, `--dry-run`, `--overwrite`, `--skip-existing`) plus one flag per
other `load()` argument, `--<option>` and `-<first letter>`. Stage files take the
same options as keys. So every option needs:

- a type hint: `int`, `str`, an `Enum`, or a tuple of `str` or of a `Literal`,
  which limits its values;
- a default;
- an entry in the Parameters section of `load()`'s numpy docstring, which is the
  flag's help.

## 1. Gridded (raster) datasets: xarray → Zarr

Examples: `prism.py`, `hrrr.py`.

1. **Open the source** into an `xr.Dataset`/`DataArray` with a spatial CRS
   readable by `rioxarray` (`da.rio.crs`). Use whatever access pattern fits the
   source (STAC + `odc.stac`, Icechunk, `rioxarray.open_rasterio`, …); see
   `open_prism` / `open_hrrr`. Planetary Computer collections can use `stac.search`
   and `stac.load` (signed, lazy, native grid; see `sentinel1_rtc.py`).
2. **Standardize dims and units, not names**: spatial dims `x`, `y`; time dim
   `time`; SI units (mm, degC). Variables keep the source's names (`ppt`, `tp`).
   If another dataset has the same quantity, add it to the module's `COMMON_NAMES`
   (staged name → common name, e.g. `{"tp": "precipitation"}`): the runner stamps
   it as `gaia:common-name`, which `datasets.cf_criteria()` lets cf-xarray select
   (`ds.cf["precipitation"]`). Set a CF `standard_name` where one fits.
3. **Clip to the AOI** with `io.clip_to_aoi(da, aoi)`. This reprojects the AOI
   into the raster's CRS and keeps every touched cell; the raster stays on its
   **native grid**. Putting products on a common grid is `gaia align`'s job (see
   [reprojecting.md](reprojecting.md)).
4. **Describe the time step.** If each value covers a period rather than an
   instant, return `io.with_time_bounds(da, start_hours, end_hours,
   cell_methods)`, e.g. `(-1, 0, "time: sum")` for hour-ending accumulations.
   `gaia align` relies on this to aggregate in time. Document the source's
   conventions in `dev/<dataset>.md` (see [prism.md](prism.md),
   [hrrr.md](hrrr.md)), and add a short entry to the user docs'
   [datasets page](../docs/datasets.md).
5. **Record provenance.** Set the returned Dataset's `attrs` to the source's STAC
   fields: at least `title`, `license`, `providers` and `links` (a `license` link
   plus `via` or `derived_from` to the source), with `version` if the source has
   one and `provenance.time_extent(ds)`. If the source publishes a STAC
   collection, start from `provenance.from_stac_collection(url)` (see `hrrr.py`);
   otherwise write them out (see `nlcd.py`). `io.save_datatree` adds the run
   attributes (`created`, `gaia:cli-version`, `gaia:command`); see
   [provenance-metadata.md](provenance-metadata.md).

## 2. Vector timeseries datasets: xvec → Zarr

Example: `synoptic.py`.

1. **Fetch station/point data** from the source API and keep every report: one
   `pd.DataFrame` per station, indexed by `time`, with a column per variable,
   plus an `xr.Dataset` of station variables along `station` (ids) with a
   `geometry` (shapely `Point`) coordinate.
2. **Build a CF ragged array** with `stations.ragged(stations, reports)` (see
   `synoptic.construct_ragged` and [stations.md](stations.md)). Don't resample
   onto a grid; readers use `stations.to_grid`.
3. **Return it from `load()`** with provenance attributes and its xvec geometry
   index; the store CF-encodes the geometries when it writes them.

## Registering

Add `"<name>": newdataset` (the module) to `DATASETS` in `datasets/__init__.py`.
`<name>` is the `gaia stage <name>` command, the stage-file key and the group's
last part, hyphenated if it has several words (`burn-severity`); the module name
uses underscores. A test checks every
registered module's `KIND`, `EXTENT`, `VARIANT` and `REGRID_METHOD`.

Run `pixi run test`, and add a case to `tests/test_online.py`.
