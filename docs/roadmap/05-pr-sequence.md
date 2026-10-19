# 05 · PR sequence

Each PR adds one capability or one dataset. A PR may change existing commands, flags or
output layouts when that simplifies the design ([design goal 3](README.md#design-goals)),
and must list those breaking changes in its description. Every dataset PR includes:

- a module with `load`/`stage`, `EXTENT` (a STAC Collection extent) and option declarations
- registration in `cli.py`
- `docs/<dataset>.md` covering source, units, time convention, coverage and known gaps
- an offline test with a small fixture
- an `online`-marked smoke test using the Peters Dome or Bolt Creek reference AOI ([02](02-stage-file.md#proposal))
- a README example

## Phase 0: foundations (no new data)

| # | PR | content | source to port |
|---|---|---|---|
| 0.1 | **Stage file** | `stage.yaml` parser (`pyyaml`), window grammar, inline AOIs, strict validation, `gaia stage --config` for the existing datasets (prism, hrrr, synoptic, nlcd), file and data flags mutually exclusive, `--dry-run` | landlab `config.py`/`dates.py` strictness; peters `event.yml` |
| 0.2 | **Registry + extent** | `EXTENT`/`KIND`/`TIMED` per module (STAC Collection extent, see [04](04-datasets.md#region-coverage)); multi-dataset runner with `--only`, continue-on-error and `--skip-existing`; flag-only runs build the same configuration and use the same runner; `stage all` → `conus-precip` preset | peters `stage_all.sh` behavior |
| 0.3 | **Folder output** | `-o` is a folder: `<name>.zarr` (group writes, consolidated after each write), the resolved `stage.yaml`, and a STAC Item with `data` and `config` assets as the manifest ([03](03-storage.md#stac-item-as-the-manifest)) | peters `common.write_group` / `gdf_to_xvec` fixes |
| 0.4 | **Synoptic → ragged array** | Write every report (time, running total, report type) as a CF DSG contiguous ragged array; add `stations.to_dataframe` and `stations.to_grid` helpers; add `--variables`; update examples. The hourly labeling fix is already done ([../precip-timing.md](../precip-timing.md)). | — |

0.1 and 0.2 go together naturally; 0.3 and 0.4 are independent of them.

## Phase 1: shared and landlab datasets

Ordered by the priorities in [04](04-datasets.md#prioritized-new-datasets).

| # | PR | depends on | notes |
|---|---|---|---|
| 1.1 | `dem --product {cop30,usgs10m}` | 0.2 | Priority 1: both workflows need it. Both products are public VRT mosaics read with rioxarray (Cop30 from OpenTopography, 3DEP 1/3 arc-second from USGS): no key, no new dependency, no STAC helper. Heights keep each product's geoid; vertical datum conversion in `align` is a TODO ([04](04-datasets.md#dem-resolution)). Supersedes PR #7. |
| 1.2 | `prism --variables` | — | Adds tmin and tmax alongside ppt, at 800 m as today |
| 1.3 | `soil --product solus100` | 0.2 | The 7 layers landlab requires. CONUS only. |
| 1.4 | `burn_severity` (SBS) | 0.1 | USGS BAER; filename templating from fire name and id. Skip dNBR, which nothing reads. |

## Phase 2: Alaska landslide (Peters Dome reference)

| # | PR | depends on | notes |
|---|---|---|---|
| 2.1 | `seismic` (EarthScope FDSN) | 0.3 | **No obspy.** FDSN text station metadata + dataselect via `requests`, miniSEED decoding with `pymseed`, counts ÷ sensitivity ([obspy-dependency.md](obspy-dependency.md)). Store at native rate grouped by sample rate, with an optional common grid. Channel metadata as 1-D arrays on a `channel` dimension. Distance and back-azimuth from `origin`. Size guard. Skip `seismic_lp`, which nothing reads. |
| 2.2 | `infrasound` | 2.1 | Same module, counts ÷ sensitivity to Pa, no component dimension |
| 2.3 | `era5` + atmosphere profiles (Open-Meteo) | 0.2 | Document the usage policy and rate limits. Skip ERA5-Land, which nothing reads. |
| 2.4 | `sentinel1-rtc`, `sentinel2` | 0.1 (pre/post windows) | Adds the Planetary Computer STAC helper (`stac.py`). Scene index as 1-D arrays on a `scene` dimension. Skip Landsat, which nothing reads. |

## Phase 3: conversion and reference workflows

This plan changes only gaia-cli. The other repos are references for the workflows to
support, not migration targets. Once there is an MVP, their authors can decide whether
to depend on gaia-cli.

| # | PR | notes |
|---|---|---|
| 3.1 | `gaia convert landlab` ([issue #5](https://github.com/gaia-hazlab/gaia-cli/issues/5)) | Port landlab `analysis_grid.py` (UTM grid snapped to `--res`) and the field map and ASCII writer. Reads a staged store. Consider `gdal raster reproject` on the Zarr. Landlab-specific derivations (Ksat, transmissivity, cohesion from NLCD) stay out unless a later PR adds them. |
| 3.2 | Reference stage files | `examples/peters-dome.yaml` and `examples/bolt-creek.yaml` (the two examples in [02](02-stage-file.md#proposal)), checked by `--dry-run` in CI and run as `online` tests. They show that gaia-cli covers the staging part of both workflows. |

## Later, on request

NASA Earthdata (`hls`, `opera-rtc`, `nisar`, `smap`), higher-resolution elevation
(3DEP 1 m, issue #4), vertical datum conversion in `align`, `arcticdem`, `soilgrids`,
STATSGO K-factor, dNBR, `worldcover`, `daymet`, glacier/ITS_LIVE, `hrrr --domain alaska`,
ComCat events as a vector table, `comcat_id` lookup, clearest-scene selection, a JSON
Schema for `stage.yaml`, GeoParquet tables and Iceberg/Icechunk
([03](03-storage.md#recommendation)). None of these is read by a current workflow; see
[04](04-datasets.md#what-the-workflows-actually-use).

## Ordering

Phase 0 comes first, then **Phase 1, then Phase 2**. Phase 1 is smaller and delivers
`dem`, which both workflows need. Phase 3 needs 0.3
and the datasets each reference workflow uses.

## Resolved questions

1. **Phase order:** Phase 1 before Phase 2.
2. **Waiting on open questions:** none remain in [02](02-stage-file.md#resolved-questions),
   so 0.1 can start.
3. **obspy:** gaia-cli does not depend on obspy. The pymseed and FDSN checks are done
   ([obspy-dependency.md](obspy-dependency.md#checks)).

## Open questions

None.
