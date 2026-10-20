# gaia-cli roadmap: more datasets, more regions, more workflows

Status: open questions answered, ready for implementation, 2026-10-15.

Three codebases now stage data for GAIA hazard workflows:

| repo | author | workflow | region |
|---|---|---|---|
| [gaia-cli](https://github.com/gaia-hazlab/gaia-cli) | @scottyhq | Precipitation comparison (Skagit atmospheric river) | CONUS |
| [gaia-peters-dome](https://github.com/gaia-hazlab/gaia-peters-dome) | @mdenolle | Landslide seismology and remote sensing (2023 Peters Dome avalanche) | Alaska |
| [landlab_data_prep](https://github.com/gaia-hazlab/landlab_data_prep) | @Almehedi06 | Inputs for Landlab postfire debris-flow and landslide models | CONUS (WA) |

All three follow the same pattern. Each takes **one event**: an AOI plus time windows.
It fetches a handful of sources, clips each one, and writes analysis-ready data. They
differ in config format, output format, and gridding. The goal is to make gaia-cli able
to stage the data for all three kinds of workflow. gaia-cli is developed independently:
the other two repos are references, and this plan makes no changes to them. Once there
is an MVP, their authors can decide whether to depend on gaia-cli.

## Documents

| doc | question it answers |
|---|---|
| [01-existing-work.md](01-existing-work.md) | What does each repo do, and what structure do they share? |
| [02-stage-file.md](02-stage-file.md) | What does `stage.yaml` look like: AOIs, time windows (including irregular satellite windows) and dataset options? |
| [03-storage.md](03-storage.md) | Should everything go into one Zarr store, or should timeseries go to Parquet or Iceberg? How does the STAC catalog fit? |
| [04-datasets.md](04-datasets.md) | Which minimal set of datasets and providers covers all three workflows, and where do they work? |
| [05-pr-sequence.md](05-pr-sequence.md) | In what order, one reviewable PR at a time? |
| [obspy-dependency.md](obspy-dependency.md) | Can we stage seismic and infrasound data without obspy? |

## Design goals

These goals inform the decisions below. When two options are otherwise close, choose
the one that serves these goals better.

1. **Minimize added lines of code.** Keep the codebase small enough to read in one
   sitting. Before writing a helper, check whether an existing dependency already does
   the job. Port only the parts of the other repos that a workflow actually uses, and
   prefer one command with a `--product` option over several near-duplicate modules.
2. **Add dependencies and API tokens only when they are critical.** Every new package
   and every new credential makes installs, CI and onboarding harder. Prefer anonymous,
   global sources first, and add a keyed or narrower source when a workflow needs it.
   For example, elevation is essential, so start with Copernicus GLO-30 (Cop30), which
   is global and needs no login. Defer the `bmi-topography` dependency and the
   OpenTopography API key (USGS 10 m and other products) until a workflow needs that
   extra resolution.
3. **Don't preserve backward compatibility yet.** There has been no major release, and
   the few current users can adapt. Rename commands, change flags and change output
   layouts whenever that makes the design simpler. Don't add shims, aliases or
   deprecation paths. Call out breaking changes in the PR description so users can
   update their scripts.

## Decisions

1. **Add a stage file, `stage.yaml`** (one event per file, `gaia stage --config`). It
   supplies AOIs (paths, URLs or inline), default and per-dataset time windows, and
   dataset options. It and the data-selection flags are mutually exclusive, so there is
   no merging. It covers what `gaia-peters-dome`'s `event.yml` + `stage_all.sh` and
   `landlab_data_prep`'s `base.yaml` do. See [02](02-stage-file.md).
2. **Windows can be absolute or relative to the event time, and can be named.** Satellite
   datasets accept `pre`/`post` windows and a "nearest N acquisitions" selection instead
   of one start/end. See [02](02-stage-file.md).
3. **Storage: one Zarr store per event, with CF ragged arrays for observations.**
   Rasters and waveforms stay arrays at native resolution. Irregular station
   observations use the CF DSG contiguous ragged array layout instead of a dense hourly
   grid, so nothing is lost and no dependency is added. The STAC Item gaia already
   writes becomes the event's manifest, with the resolved `stage.yaml` as an asset.
   `-o` is always a folder.
   Move tables to GeoParquet (adds pyarrow) or Iceberg only when a trigger in
   [03](03-storage.md#recommendation) occurs.
4. **Add only datasets a workflow reads, from as few providers as possible.** Many
   datasets the other repos download are never read by their analysis; these are
   deferred. Planetary Computer serves the DEM (Cop30), Sentinel-1 RTC and Sentinel-2
   with no login. EarthScope serves seismic and infrasound. Open-Meteo serves ERA5.
   NASA Earthdata (HLS, NISAR, SMAP) waits until a workflow reads those. See
   [04](04-datasets.md).
5. **Model-specific gridding is a conversion step, not a staging step.**
   `gaia convert landlab` ([issue #5](https://github.com/gaia-hazlab/gaia-cli/issues/5))
   ports `landlab_data_prep`'s single UTM analysis grid and ESRI ASCII export. Staging
   keeps native grids.

6. **No obspy.** Seismic and infrasound use FDSN web services, `pymseed` and
   counts ÷ sensitivity, keeping Python 3.15. See [obspy-dependency.md](obspy-dependency.md).
7. **Order:** Phase 0 foundations, then Phase 1 (`dem` and landlab datasets), then
   Phase 2 (Peters Dome datasets). See [05](05-pr-sequence.md).

## Resolved before implementation

- [x] Stage file name, scope and CLI behavior ([02](02-stage-file.md#resolved-questions))
- [x] Storage: A′ to start, B and C on triggers ([03](03-storage.md#resolved-questions))
- [x] Remote sensing provider: Planetary Computer ([04](04-datasets.md#remote-sensing-provider-choice))
- [x] Python 3.15 vs obspy: no obspy ([obspy-dependency.md](obspy-dependency.md#checks))
- [x] Ragged string columns in xvec + Zarr v3: tested ([03](03-storage.md#resolved-questions))
