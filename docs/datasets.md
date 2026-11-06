---
title: Datasets
---

| dataset | what | coverage | grid | time | group |
|---|---|---|---|---|---|
| [`synoptic`](#synoptic) | weather station reports | mostly US, Alaska included | stations | every report | `stations/synoptic` |
| [`prism`](#prism) | PRISM daily precipitation and temperature | CONUS, 1981 on | 800 m | daily | `rasters/prism` |
| [`hrrr`](#hrrr) | NOAA HRRR analysis | CONUS, 2014-10 on | 3 km | hourly | `rasters/hrrr` |
| [`era5`](#era5) | ECMWF ERA5 reanalysis, single levels | global, 1940 on (a few months behind) | 0.25° | hourly | `rasters/era5` |
| [`stage4`](#stage4) | NCEP Stage IV precipitation | CONUS, 2002 on; Alaska, 2020-07 on | 4.8 km | daily, 6-hourly, hourly | `stage4/<accumulation>` |
| [`nlcd`](#nlcd) | Annual NLCD land cover | CONUS, 1985–2024 | 30 m | static | `rasters/nlcd` |
| [`dem`](#dem) | elevation | global (`cop30`), US (`usgs10m`) | 30 m, 10 m | static | `dem/<product>` |
| [`soil`](#soil) | soil properties | CONUS | 100 m, 30 m | static | `soil/<product>` |
| [`burn-severity`](#burn-severity) | BAER soil burn severity | CONUS, fires 2017 on | 20–30 m | static | `rasters/burn-severity` |
| [`sentinel1-rtc`](#sentinel1-rtc) | Sentinel-1 radar backscatter | global, 2014-10 on | 10 m | each pass | `rasters/sentinel1-rtc` |
| [`sentinel2`](#sentinel2) | Sentinel-2 surface reflectance | global, 2015-06 on | 10 m | each pass | `rasters/sentinel2` |

A dataset that the AOI or window misses (for example `prism` in Alaska) is skipped
with a warning. Licenses differ: check each group's `license` attribute before sharing
a store.

Variables keep their source's names (PRISM `ppt`, ERA5 `tp`), in SI units, in staged
and aligned stores. Quantities that several datasets have also carry a common name in
`gaia:common-name`:

| common name | variables |
|---|---|
| `precipitation` (mm) | `prism` `ppt`, `hrrr` `precipitation_surface`, `stage4` `apcp`, `era5` `tp`, `synoptic` `precip_accum` (a running total) |
| `temperature` (degC, 2 m air) | `prism` `tmean`, `hrrr` `temperature_2m`, `era5` `t2m` |

To select by common name, use [cf-xarray](https://cf-xarray.readthedocs.io) with
gaia's criteria:

```python
import cf_xarray as cfxr
from gaia_cli.datasets import cf_criteria

cfxr.set_options(custom_criteria=cf_criteria())
ds.cf["precipitation"]   # ppt in a PRISM group, tp in an ERA5 group
```

`ds.cf[...]` also takes CF standard names, e.g. `lwe_thickness_of_precipitation_amount`.

(synoptic)=
## Synoptic stations

Station reports from the [Synoptic Data](https://synopticdata.com) API. Requires a
`SYNOPTIC_TOKEN`.

- `precip_accum` (default) is staged as a running total in mm. Other Synoptic
  variables, such as `air_temp`, keep their Synoptic units.
- Every report is kept, as a CF ragged array. `gaia_cli.stations.to_grid(ds)` puts the
  reports on a regular hourly (or daily) grid, and `to_dataframe(ds)` gives one row per
  report.
- License: set by each station network.

(prism)=
## PRISM

[PRISM](https://prism.oregonstate.edu) AN81d daily grids, 800 m, EPSG:4269.

- Variables: `ppt` (default, mm), `tmin`, `tmax`, `tmean` (degC).
- A **PRISM day** D covers 12:00 UTC on D−1 to 12:00 UTC on D. The value is stamped
  D 00:00, and `time_bounds` gives the true period.
- Read from a public cache of PRISM's grids, browsable in
  [STAC Browser](https://browser.moregeo.it/external/s3.kopah.uw.edu/gaia/public/prism/collection.json).
  Recent grids are re-released by PRISM for about 6 months; each value records its
  `prism_release` (8 is final).
- License: [PRISM terms of use](https://prism.oregonstate.edu/terms/).

(hrrr)=
## HRRR

NOAA High-Resolution Rapid Refresh analysis from
[dynamical.org](https://dynamical.org/catalog/noaa-hrrr-analysis/), 3 km, Lambert
Conformal Conic, hourly.

- Variable: `precipitation-surface` (default, staged as `precipitation_surface` in mm
  per hour), or any other in `gaia stage hrrr --help`.
- Precipitation is **hour-ending**: the value at 13:00 fell between 12:00 and 13:00
  UTC.
- Wind components are grid-relative, not east/north.
- License: CC-BY-4.0.

(era5)=
## ERA5

ECMWF [ERA5](https://doi.org/10.1002/qj.3803) hourly reanalysis on single levels,
read from [Earthmover](https://earthmover.io)'s Icechunk copy, 0.25°, EPSG:4326.
Requires a free Earthmover account: set `ARRAYLAKE_TOKEN` to an API key (Organization
Settings → API Clients), or run `arraylake auth login` once.

- Variables keep their ECMWF short names. By default: `tp` (total precipitation,
  mm), `sf` (snowfall, mm of water equivalent), `t2m` (2 m temperature, degC),
  `swvl1` and `swvl2` (soil water, 0–7 and 7–28 cm, m³ m⁻³) and `sd` (snow depth,
  mm of water equivalent).
- Others worth adding with `-v`: `d2m`, `u10`, `v10`, `sp` (to compare with HRRR;
  ERA5 winds are east/north), `ssrd` (solar radiation, for snowmelt), `stl1` (soil
  temperature, frozen ground), `cp`/`lsp` (convective and large-scale
  precipitation), `tcwv` (column water vapour). See `gaia stage era5 --help` for
  all of them.
- Not in this copy: daily `t2m` minimum and maximum (`gaia align` can derive them
  from hourly `t2m`), deeper soil water (`swvl3`, `swvl4`), runoff, snowmelt,
  evaporation, snow density and moisture transport (IVT).
- Every value is **hour-ending**, like HRRR: `tp` at 13:00 fell between 12:00 and
  13:00 UTC.
- The free copy is updated quarterly, so the latest months are missing.
- On glaciers ERA5 models a deep snowpack, so `sd` stays in the metres all year.
- `orography` (m) is the height of ERA5's smoothed terrain in each cell, the level
  `t2m` is valid at. In the mountains it can be far below the real terrain; use it
  to correct temperatures for height. `gaia align` drops it.
- License: CC-BY-4.0; keep the Copernicus `acknowledgement` attribute.

(stage4)=
## Stage IV

NCEP/EMC [Stage IV](https://gdex.ucar.edu/datasets/d507005/) radar and gauge
precipitation analyses, about 4.8 km, polar stereographic.

- Staged as `apcp` (GRIB's APCP, total precipitation) in mm, one group per
  `--accumulation`: `daily` (default,
  `stage4/daily`), `sixhourly` or `hourly`. The daily value stamped D covers 12:00 UTC on
  D−1 to 12:00 UTC on D, the same as a PRISM day. 6-hourly and hourly values are
  stamped at the end of their period, like HRRR.
- Hourly grids have no data over parts of the western US (including Washington) before
  2018: those steps are NaN, with a warning. Alaska has no hourly grids.
- The grid follows the AOI: CONUS from 2002, or Alaska from 2020-07-01. An Alaska
  window starting earlier is an error.
- Fetched from the NCAR GDEX archive, one day at a time, into a local cache
  (`./.cache/gaia/stage4` under the directory gaia runs in, or `STAGE4_CACHE`). Days
  missing from the archive are NaN.
- License: CC-BY-4.0.

(nlcd)=
## NLCD

USGS [Annual National Land Cover Database](https://www.usgs.gov/centers/eros/science/annual-national-land-cover-database),
30 m, one year (`--year`, default 2021). Classes are integer codes with CF
`flag_values`/`flag_meanings`. License: U.S. public domain.

(dem)=
## DEM

| `--product` | source | heights | coverage |
|---|---|---|---|
| `cop30` (default) | [Copernicus GLO-30](https://doi.org/10.5069/G9028PQB), surface model, 30 m | EGM2008 | global |
| `usgs10m` | [USGS 3DEP 1/3 arc-second](https://doi.org/10.5069/G98K778D), bare earth, 10 m | NAVD88 | US |

Heights are not converted between vertical datums. Each product is its own group, so
one store can hold both.

(soil)=
## Soil

| `--product` | source | variables |
|---|---|---|
| `solus100` (default) | [SOLUS100](https://doi.org/10.15482/USDA.ADC/25033856), 100 m, CC-BY-4.0 | surface clay, sand, silt, bulk density, pH, CEC, and depth to bedrock |
| `statsgo` | [STATSGO K-factor](https://doi.org/10.5066/P13WAPYV), 30 m from 1:250,000 polygons, public domain | `kffact`, soil erodibility |

(burn-severity)=
## Burn severity

USDA Forest Service [BAER](https://burnseverity.cr.usgs.gov/baer/) soil burn severity
maps, chosen by the AOI. Where several fires intersect it, choose one with `--fire`
(part of its name) and `--year`. Classes: 1 unburned to very low, 2 low, 3 moderate,
4 high. Each fire keeps its own grid. License: U.S. public domain.

(sentinel1-rtc)=
## Sentinel-1 RTC

Sentinel-1 radiometrically terrain-corrected backscatter from
[Planetary Computer](https://planetarycomputer.microsoft.com/dataset/sentinel-1-rtc),
10 m, in the scenes' UTM zone. No account needed.

- `vv` and `vh` (default; `--polarizations` also takes `hh`, `hv`): gamma0 as
  **linear power**; `10 * log10(vv)` gives dB.
- One time step per satellite pass, at the pass's start time. Each step records its
  `platform`, `orbit_state`, `relative_orbit` (track) and `absolute_orbit`, and
  `--orbit-state` and `--relative-orbit` keep only some of them.
- Sentinel-1B stopped in December 2021, so 2022–2024 has a 12-day repeat in most
  places.
- License: CC-BY-4.0.

```bash
gaia stage sentinel1-rtc -i aoi.geojson -s 2023-09-01 -e 2023-09-30 -o out --relative-orbit 131
```

In a stage file, keep only the last pass before and the first pass after an event, on
each track that has both:

```yaml
time: 2023-09-14T02:06:49Z
datasets:
  sentinel1-rtc:
    window:
      pre:  {before: 30d, select: last}
      post: {after: 30d,  select: first}
    match: relative_orbit
```

Each pass then has a `window` coordinate, `pre` or `post`.

(sentinel2)=
## Sentinel-2

Sentinel-2 Level-2A surface reflectance from
[Planetary Computer](https://planetarycomputer.microsoft.com/dataset/sentinel-2-l2a),
at 10 m in the scenes' UTM zone. No account needed.

- `--bands` by common name: `blue`, `green`, `red`, `nir`, `swir16` and `scl` by default;
  also `coastal`, `rededge1`–`3`, `nir08`, `nir09`, `swir22`. Reflectance is 0–1;
  20 m and 60 m bands are upsampled to 10 m.
- `scl` is the scene classification (clouds, shadows, snow, water, …) as CF flags; mask
  clouds with `ds.where(~ds.scl.isin([3, 8, 9, 10]))`.
- One time step per satellite pass, its tiles mosaicked. `--max-cloud-cover` skips
  tiles cloudier than a percentage; that is the cloud cover of the whole 110 km tile,
  not the AOI's.
- Like `sentinel1-rtc`, a stage file can keep the last pass before and the first after
  an event with `pre`/`post` windows, as in the [Sentinel-1 example](#sentinel1-rtc).
- License: Copernicus Sentinel data terms (free, with attribution).

```bash
gaia stage sentinel2 -i aoi.geojson -s 2023-09-01 -e 2023-09-30 -o out --max-cloud-cover 50
```
