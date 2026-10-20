# 01 · Existing work and shared structure

## gaia-cli (`main` @ `d2dae7b`)

- Commands: `gaia stage {synoptic,prism,hrrr,nlcd,all}`, `gaia align`, `gaia prism {fill,rebuild-catalog}`.
- Each dataset module has `load(aoi, start, end, ...) -> xarray` and a cyclopts `stage`
  entrypoint with `-i/-s/-e/-o/--overwrite`. The patterns are documented in
  [../adding-a-dataset.md](../adding-a-dataset.md).
- Rasters are clipped on their native grid (`io.clip_to_aoi`). Period data carries CF
  `time_bounds` + `cell_methods`. `gaia align` regrids to a common grid and time step
  afterwards.
- Stations are a dense `(station, time)` array on a uniform hourly grid. Synoptic data
  is put on the hourly axis as of the latest report (with `report_time`; see
  [../precip-timing.md](../precip-timing.md)). The array is xvec CF-encoded into Zarr.
- `stage all` writes a single DataTree Zarr with `stations/synoptic`, `rasters/prism`
  and `rasters/hrrr`. Metadata is always consolidated, because anonymous Backblaze B2
  reads cannot LIST.
- Open: [PR #7](https://github.com/gaia-hazlab/gaia-cli/pull/7) (OpenTopography DEM, by @Almehedi06) and
  [issue #5](https://github.com/gaia-hazlab/gaia-cli/issues/5). Issue #5 proposes GDAL-style
  subcommands, `gaia stage landlab` to make a consolidated Zarr, and
  `gaia convert landlab` to make Landlab ASCII inputs, possibly via `gdal raster reproject`
  on the staged Zarr.
- CONUS-only today. PRISM, HRRR (dynamical.org) and Annual NLCD are all CONUS.

## gaia-peters-dome (`main` @ `a9401c6`, open PRs #2, #3)

Covered in detail in the first plan. In summary:

- `peters-dome stage <name>` uses the gaia-cli contract and imports `gaia_cli.io`
  directly. PR #3 moves everything into **one Zarr store, one group per dataset**, with
  `--group` and `--skip-existing`.
- **There are three AOIs and six time windows**, kept in `config/event.yml` and
  `scripts/stage_all.sh`. The windows range from 40 minutes for seismic to 18 years for
  DEM strips. Event metadata (origin time, location, DOIs, UTM CRS) also lives in
  `event.yml`.
- New modalities:
  - EarthScope FDSN seismic and infrasound, response-removed and interpolated onto a
    common 20 Hz grid.
  - Planetary Computer Sentinel-1 RTC, Landsat and COP-DEM; Earth Search Sentinel-2;
    PGC ArcticDEM.
  - Open-Meteo ERA5 and ERA5-Land.
  - Glacier, geology and soil context layers.
  - A USGS data release.
- Its datacube is built by a custom script, not `gaia align`, because `align` aggregates
  to the coarsest time step (12-day SAR vs hourly ERA5).
- Synoptic is used via `gaia_cli.synoptic.load`.

## landlab_data_prep (`master` @ `903a864`, no PRs)

Builds Landlab inputs for postfire debris-flow and landslide models.
`landlab-prep-<source> --config base.yaml`; 12 console scripts; about 5.9k lines.

**Config.** There is one YAML file per event, kept out of git. It has `aoi`, `paths`
(`output_dir`, `cache_dir`), `fire` (name, id, post-image date), `raster` (target_res,
resampling), and one block per source.

- Validation is **strict**: unknown keys are errors at every level, and every problem is
  reported at once before any work starts. Each command requires only its own sections,
  but checks every section that is present.
- Time windows are per source: `prism.start/end`, `synoptic.start/end`, `nisar.start/end`
  and `smap.start/end`. HLS uses **named windows**:
  `remote_sensing.windows.pre/post.{start,end}`, validated so that pre ends before post
  starts. Dates are inclusive UTC days.

**Sources**, each in its own subpackage with `config.py` (a strict parser that returns a
frozen dataclass), `core.py` and `cli_run.py`:

| source | access | notes |
|---|---|---|
| DEM | OpenTopography via bmi-topography (API key) | Same approach as gaia-cli PR #7 |
| Burn severity (SBS), dNBR | USGS BAER/MTBS zip downloads, local fallback | Name templating from `fire.name`/`fire.id` |
| Soils | SOLUS100 COGs on GCS (CONUS), STATSGO K-factor (ScienceBase) | `missing_values` codes |
| Land cover | Annual NLCD zip | Same product as `gaia stage nlcd` |
| HLS vegetation indices | NASA CMR STAC (LPCLOUD), `HLSL30_VI`/`HLSS30_VI` v2.0 | Median pre/post composites, dNBR, rdNBR, clear-day counts |
| PRISM | NACSE web service, 800 m or 4 km, ppt/tmin/tmax | Cached by PRISM release version |
| Synoptic | `/stations/precip` with `pmode=intervals` | Interval totals, nearest-station gridding, raw JSON cached with a sha256 in the manifest |
| NISAR GCOV, SME2 | CMR granule search + HTTP range reads of HDF5 | Prefers the Provisional release, falls back to Beta |
| SMAP L4 (SPL4SMGP v008) | CMR + range reads | 3-hourly, 9 km |

**Shared infrastructure** worth reusing:
- `earthdata.py`: Earthdata login (`EARTHDATA_TOKEN` or `~/.netrc`), paged CMR granule
  search, and a `RemoteFile` that reads by HTTP range with a block cache. It already
  serves HLS, NISAR and SMAP.
- `analysis_grid.py`: the single analysis grid. It uses the AOI's UTM zone, with bounds
  snapped outward to `target_res`, and `align_to_grid` is the only resampler. This is
  exactly what `gaia convert landlab` needs.
- `downloads.py`: a cached download with an extraction marker.
- **Provenance.** Every source writes a `*_manifest.json` (request parameters, cache
  files, sha256, method notes) and a `*_index.csv` with one row per scene, interval or
  step: release, track, pass direction, AOI coverage and mean.

**Outputs.** GeoTIFFs and ESRI ASCII on the analysis grid, plus CSVs. There is no Zarr.

## Shared structure

| concern | gaia-cli | gaia-peters-dome | landlab_data_prep | proposed shared form |
|---|---|---|---|---|
| Unit of work | AOI + one window | event: 3 AOIs, 6 windows | event: 1 AOI, per-source windows, pre/post | **`stage.yaml`** ([02](02-stage-file.md)) |
| Per-source code | `load()` + `stage()` | same (gaia contract) | `config.py` + `core.py` + `cli_run.py` | gaia `load()` + `stage()`, with options declared so `stage.yaml` can validate them |
| Config validation | cyclopts types | cyclopts types | strict, report everything at once | strict, report everything at once (port the landlab approach) |
| Time convention | `-s/-e` timestamps | per-modality windows, CF `time_bounds` | inclusive UTC days | one explicit rule ([02](02-stage-file.md#time-conventions)) |
| Grid | native, then `gaia align` | native, then custom cube | one UTM grid up front | native at stage time; `align`/`convert <model>` afterwards |
| Station data | dense hourly xvec Zarr | same, plus seismic at 20 Hz | per-interval nearest-station rasters + raw CSV | CF ragged array of every report + derived dense views ([03](03-storage.md)) |
| Scene metadata | n/a | SAR orbit coordinates on `time` | `*_index.csv` per product | 1-D scene index arrays ([03](03-storage.md)) |
| Provenance | few attrs | `provenance()` attrs + event DOIs | `manifest.json` with request and cache hashes | STAC Item per event (the manifest) + attrs per group |
| Caching | PRISM S3 cache | geology zip only | `cache_dir` for everything | optional `--cache-dir`; raw API responses cached |
| Auth | `SYNOPTIC_TOKEN`, `AWS_PROFILE` | none beyond gaia | `SYNOPTIC_TOKEN`, `EARTHDATA_TOKEN`/netrc, OpenTopography key | the same environment variables, checked up front by each dataset |
| Python / env | pixi, Python 3.15 | pixi, Python 3.12 | pip/conda, Python 3.10–3.13 | pixi; Python 3.15, no obspy ([obspy-dependency.md](obspy-dependency.md)) |

### Duplicated datasets to reconcile

| dataset | gaia-cli | other repo | resolution |
|---|---|---|---|
| PRISM | daily ppt at 800 m from the S3 COG cache, filled from PRISM's file server (`data.prism.oregonstate.edu`) | landlab: the same zips via the NACSE web service (`services.nacse.org`), no shared cache; 800 m by default, 4 km optional; ppt+tmin+tmax | Same upstream data and release numbers in both. Add `--variables` to gaia `prism` and keep the S3 cache. Add 4 km only if a workflow needs it. |
| Synoptic | `/stations/timeseries`: Synoptic's derived cumulative precip (`precip=1`) when the account tier allows it, else gaia differences raw `precip_accum` with a reset heuristic; then the running total is put on an hourly grid as of the latest report | landlab: `/stations/precip?pmode=intervals`, one total per requested interval with `count`, `first_report`, `last_report` and `report_type`; only complete intervals are used | **Keep gaia's timeseries approach.** The precip service's intervals follow report times, not clock hours, and it drops intervals across reporting gaps. gaia's hourly labeling is now fixed (as-of, plus `report_time`); see [../precip-timing.md](../precip-timing.md). N-hour totals for Landlab derive from the staged running total, and gridding goes in `gaia convert landlab`. |
| NLCD | Annual NLCD COG via OSN | landlab: the same product via zip | landlab uses gaia `nlcd`. |
| DEM | PR #7 OpenTopography | landlab: OpenTopography; peters: COP-DEM, ArcticDEM | Start with `dem` = Cop30 (global, no key). Try 3DEP 10 m from Planetary Computer for CONUS landlab runs. Add `--product usgs10m` (OpenTopography, `bmi-topography`, API key) and `arcticdem` only if a workflow needs them ([design goal 2](README.md#design-goals), [04](04-datasets.md#dem-resolution)). |
| Soils | — | landlab: SOLUS100 (CONUS); peters: SoilGrids (global) | `soil --product solus100` first; landlab reads all 7 layers. Peters SoilGrids is only in its unread datacube, so defer it. |
| Optical | — | landlab: HLS VI; peters: S2 + Landsat | Nothing reads landlab's HLS output or peters' Landsat. Peters reads S2, so start with `sentinel2` from Planetary Computer. See [04](04-datasets.md#what-the-workflows-actually-use). |
| SAR | — | landlab: NISAR GCOV; peters: S1 RTC | Nothing reads NISAR. Peters reads S1 RTC, so start with `sentinel1-rtc` from Planetary Computer. |
| Soil moisture | — | landlab: SMAP L4, NISAR SME2; peters: ERA5-Land layer | Nothing reads SMAP, NISAR SME2 or ERA5-Land. Defer. |

## Resolved questions

1. **PRISM variables for landlab:** include tmin and tmax (`prism --variables`, PR 1.2).
2. **Unused `fire-debrisflow-ml` inputs** (slope/aspect, dNBR): stay deferred. The
   STATSGO K-factor became a required landlab field (`kffact`) and is step 1.5,
   `soil --product statsgo`.
3. **Role of the other repos:** gaia-cli is developed independently. gaia-peters-dome
   and landlab_data_prep inform which workflows to support; this plan makes no changes
   to them. After an MVP, their authors can decide whether to adopt gaia-cli.

## Open questions

None.
