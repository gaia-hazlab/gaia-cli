# 04 · Minimal dataset set, providers and coverage

Selection rules:

1. Add a dataset only when a workflow **reads** it. Being downloaded by one of the other
   repos is not enough (see [What the workflows actually use](#what-the-workflows-actually-use)).
2. Prefer anonymous, global providers, and a provider that serves several needed
   datasets ([design goal 2](README.md#design-goals)).
3. Every dataset declares its spatial and temporal coverage as a STAC Collection
   `extent` ([Region coverage](#region-coverage)), so an Alaska AOI fails fast on
   CONUS-only sources.

## What the workflows actually use

Checked on 2026-10-15 by tracing which staged outputs are read by analysis code, not
just downloaded.

**gaia-peters-dome** (`main` @ `a9401c6`). There are no notebooks. The only consumer is
`scripts/make_report.py`, which makes 11 figures. `scripts/build_datacube.py` writes a
`datacube` group, but nothing reads it.

| status | datasets |
|---|---|
| Used in figures | seismic (20 Hz velocity: map, record sections), infrasound (record section, I53 traces), ERA5 via Open-Meteo (daily temperature and precip at the source), Open-Meteo atmosphere profiles, Sentinel-1 RTC (before/after log-ratio), Sentinel-2 L2A (clearest before/after scenes) |
| Context or overlays only | ITS_LIVE speed map, RGI outlines, landslide outlines |
| Event-specific, plotted as published | force-time function, SHALTOP run, SETSM post-event DEM |
| Only in the unread `datacube` | COP-DEM, ArcticDEM mosaic, SoilGrids, geology units |
| Never read | Landsat, ERA5-Land, ArcticDEM strips, Synoptic, long-period seismic (`seismic_lp`), GLIMS, GlaThiDa, geology contacts/faults |
| Not staged at all | WorldCover |

**landlab_data_prep.** Its output fields match Landlab `LandslideProbability` inputs.
The model-facing consumer found locally is `fire-debrisflow-ml`, whose 23 features come
from these layers.

| status | datasets |
|---|---|
| Used (Landlab fields or derived soil properties) | DEM (`USGS10m` via bmi-topography; 10 m analysis grid), burn severity SBS, SOLUS100 (all 7 layers are required), NLCD (cohesion, plant functional type, water mask), STATSGO K-factor (`kffact`, a field in the [contract](https://github.com/gaia-hazlab/landlab_data_prep/blob/master/docs/landlab_workflow_contract.md)) |
| Used outside this repo | PRISM ppt/tmin/tmax, said to be read by `landlab_debrisflow` (not verified; repo not available locally) |
| Never read | HLS, NISAR, SMAP, Synoptic, dNBR |

The landlab DEM path has no hard dependency on 3DEP. `dem_type` is a free string passed
to bmi-topography, with no resolution checks or CONUS logic. The CONUS lock-in comes
from SOLUS100, NLCD, BAER burn severity and PRISM.

**Neither repo reads HLS, NISAR, SMAP or Synoptic.** Synoptic is used only by gaia's own
precipitation comparison.

## Prioritized new datasets

| priority | dataset | provider | used by | notes |
|---|---|---|---|---|
| **1** | `dem` (Cop30, 3DEP 10 m) | OpenTopography, USGS (VRTs) | landlab (core input); peters (cube only) | The only new dataset both workflows need. See [DEM resolution](#dem-resolution). |
| 2 | `prism --variables` (tmin, tmax) | existing | landlab (via `landlab_debrisflow`) | Small change to an existing command |
| 2 | `soil --product solus100` (7 layers) | SOLUS100 COGs on GCS | landlab | CONUS only |
| 2 | `soil --product statsgo` (K-factor) | USGS COG on S3 | landlab | CONUS only; see [../soil.md](../soil.md#statsgo-k-factor) |
| 2 | `burn_severity` (SBS) | USFS BAER image service | landlab | Fire-specific; found by AOI (revised in step 1.4, see [burn_severity.md](../burn_severity.md)) |
| 2 | `seismic`, `infrasound` | EarthScope FDSN | peters | Core of the landslide analysis. See [obspy-dependency.md](obspy-dependency.md). |
| 2 | `era5` + atmosphere profiles | Open-Meteo | peters | One provider, no key |
| 2 | `sentinel1-rtc`, `sentinel2` | Planetary Computer | peters | Anonymous. Peters uses PC for S1 and Earth Search for S2; PC also serves S2 L2A, so one provider covers both. |
| defer | HLS, OPERA RTC-S1, NISAR, SMAP, Landsat, ERA5-Land, ArcticDEM, WorldCover, SoilGrids, dNBR, geology, glacier datasets, `seismic_lp` | — | none | Staged but never read. Add when a workflow reads them. |

Priority 1 comes first because both workflows need it. Within priority 2, the landlab
datasets come first, then Peters Dome ([05](05-pr-sequence.md#ordering)).

### DEM resolution

Cop30 works for landlab's code as-is, but landlab's analysis grid is 10 m, so a 30 m
DEM would be cubic-upsampled. Slope drives landslide probability, and a 30 m source
smooths steep slopes. The ML target (`dem_diff.tif`) is also a DEM difference. So Cop30
is the right global default, but CONUS landlab runs would lose detail.

**Decision (revised in step 1.1):** `dem` has two products, both public VRT mosaics of
COGs read with rioxarray, with no key and no new dependency: `--product cop30` (default,
global, [OpenTopography](https://opentopography.s3.sdsc.edu/raster/COP30/COP30_hh.vrt),
EGM2008 heights) and `--product usgs10m` (USGS 3DEP 1/3 arc-second,
[The National Map](https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/USGS_Seamless_DEM_13.vrt),
NAVD88 heights). This replaces the Planetary Computer plan and PR #7 (`bmi-topography`
and an API key). Heights stay on each product's geoid; converting to a common vertical
datum in `gaia align` is a TODO. See [dem.md](../dem.md).

## By provider

| provider | access / auth | datasets (priority 1–2) | used by | coverage |
|---|---|---|---|---|
| **Planetary Computer** | anonymous SAS signing | `sentinel1-rtc`, `sentinel2` | peters | global |
| OpenTopography, USGS (VRTs) | none | `dem` (`cop30`, `usgs10m`) | landlab, peters | global; `usgs10m` is US |
| **EarthScope FDSN** | none | `seismic`, `infrasound` | peters | global networks; AK network well covered |
| **Open-Meteo** | none; usage policy, rate limits | `era5`, atmosphere profiles | peters | global |
| SOLUS100 (GCS) | none | `soil --product solus100` | landlab | **CONUS** |
| USFS BAER (GeoPlatform image service) | none | `burn_severity` | landlab | USFS fires from 2017, CONUS |
| PRISM (gaia S3 cache / NACSE) | none | `prism` | gaia, landlab | **CONUS** |
| dynamical.org | none | `hrrr` | gaia | **CONUS** |
| USGS (OSN, MRLC) | none | `nlcd` | gaia, landlab | **CONUS** |
| Synoptic | `SYNOPTIC_TOKEN` | `synoptic` | gaia | US incl. Alaska, and beyond |

No priority 1–2 dataset needs a new credential.

### Later, on request

| dataset | provider | note |
|---|---|---|
| `hls`, `opera-rtc`, `nisar`, `smap` | NASA Earthdata (`EARTHDATA_TOKEN` or `~/.netrc`) | landlab has working code in `earthdata.py`, but nothing reads the outputs. Adds a credential for every user and CI. |
| `dem` vertical datum conversion | — | EGM2008 ↔ NAVD88 in `gaia align` (TODO in `align.py`) |
| `arcticdem` (mosaic, strips) | PGC STAC | > 60°N; multitemporal elevation change |
| `soil --product soilgrids` | ISRIC | SoilGrids is global but mostly NaN over ice and rock |
| `worldcover` | Planetary Computer | Global land cover, if a non-CONUS workflow needs NLCD's role |
| `daymet` | ORNL DAAC (Earthdata) / Planetary Computer (to 2020) | Includes Alaska |
| `itslive`, glacier outlines, geology | ITS_LIVE S3, GLIMS, USGS | Context maps only |
| `hrrr --domain alaska` | NOAA AWS GRIB2 | No ARCO store; costly to build |

## Remote sensing provider choice

The usage check settles this. The imagery the workflows actually read is Sentinel-1 RTC
and Sentinel-2 L2A (Peters Dome), both available anonymously. Nothing reads HLS, OPERA,
NISAR or SMAP.

**Decision:** use Planetary Computer for `sentinel1-rtc` and `sentinel2` (`dem` reads public VRTs instead; see [DEM resolution](#dem-resolution)).
That is one provider and one STAC helper with no credential. Revisit NASA Earthdata when
a workflow reads HLS, NISAR or SMAP. The earlier comparison still applies then: HLS is
one harmonized 30 m product, OPERA and NISAR come from the same family, and the OPERA
RTC-S1 archive start (2022-01 or 2023-10) needs checking against any pre-2024 event.

## Region coverage

Each dataset module exports `EXTENT`, a **STAC Collection `extent` object**, rather than a
gaia-specific coverage type:

```python
EXTENT = {
    "spatial": {"bbox": [[-125.02, 24.06, -66.48, 49.94]]},           # [west, south, east, north], EPSG:4326
    "temporal": {"interval": [["1981-01-01T00:00:00Z", None]]},        # None = open-ended
}
```

The schema follows the STAC Collection spec:

- The first `bbox` is the overall extent. Any further bboxes describe sub-regions, for
  example CONUS and Alaska for a US-wide source, and the check uses them when present.
- A bbox that crosses the antimeridian has `west > east`, as in STAC. This matters for
  the Aleutians.
- The temporal interval covers the archive (SMAP from 2015-03, NISAR from 2025-10), so
  a window outside it fails validation too. `None` leaves an end open, and static
  datasets use `[[None, None]]`.

**Where the values come from.** Copy them from the provider's STAC Collection when one
exists, and write them by hand in the same schema when it doesn't. Store them as a
literal in the module with a link to the source, so validation and `--dry-run` need no
network call ([02](02-stage-file.md#validation)).

| dataset | `EXTENT` source |
|---|---|
| `sentinel1-rtc`, `sentinel2` | Planetary Computer collection (`sentinel-1-rtc`, `sentinel-2-l2a`) |
| `hrrr` | dynamical.org collection, which `hrrr.py` already reads for provenance |
| `prism` | the gaia PRISM cache's own collection (`prism_cache/stac.py`) |
| `dem`, `nlcd`, `soil` (SOLUS100, STATSGO), `burn_severity`, `synoptic`, `era5`, `seismic`, `infrasound` | written by hand in the same schema (the provider has no STAC Collection, or none was found yet) |

How the extent is used:

- `gaia stage --config` skips datasets whose extent the AOI or window misses, with a
  warning, and `--dry-run` lists them.
- A single `gaia stage prism` with an Alaska AOI exits with a clear error.
- The check is coarse on purpose: a bbox intersection, not a polygon test. A CONUS bbox
  also covers parts of Canada and Mexico, and those requests fail later at the provider
  with an empty result.
- The staged group's own STAC Item ([03](03-storage.md#stac-item-as-the-manifest))
  records the extent of the data actually written, so input and output use the same
  schema.

Summary by region:

| dataset | CONUS | Alaska | global |
|---|---|---|---|
| prism, hrrr, nlcd, solus100 | ✓ | ✗ | ✗ |
| burn_severity | USFS fires | ✗ | ✗ |
| synoptic | ✓ | ✓ | partial |
| dem (Cop30), sentinel1-rtc, sentinel2, era5, seismic, infrasound | ✓ | ✓ | ✓ |

## Presets

Presets are built-in stage files without AOIs or times:

| preset | datasets | replaces |
|---|---|---|
| `conus-precip` | synoptic, prism, hrrr | `gaia stage all` |
| `landslide` | seismic, infrasound, era5, sentinel1-rtc, sentinel2, dem | peters `stage_all.sh` (the general part) |
| `landlab-postfire` | dem, nlcd, burn_severity, soil, prism | landlab `base.yaml` |

## Resolved questions

1. **3DEP 10 m from Planetary Computer** would work, but `dem` starts with Cop30 only.
   `--product 3dep10m` comes later; PR #7 stays deferred.
2. **Cop30 for landlab:** Cop30 is the global default other workflows need. The
   smoothing of slopes on a 10 m grid is a known limitation, and higher-resolution
   elevation is a goal for later.
3. **Planetary Computer S1 RTC and S2** are expected to cover what Peters Dome uses.
   The `online` test for PR 2.4 confirms it.
4. **Other repos:** unchanged. This plan works only in gaia-cli.

## Open questions

None.
