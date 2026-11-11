---
title: Datasets
---

| dataset | what | coverage | grid | time | group |
|---|---|---|---|---|---|
| [`synoptic`](#synoptic) | weather station reports | mostly US, Alaska included | stations | every report | `stations/synoptic` |
| [`nrcs`](#nrcs) | NRCS SNOTEL and SCAN station reports: soil moisture, snow | western US, Alaska; SCAN US-wide | stations | hourly | `stations/nrcs` |
| [`prism`](#prism) | PRISM daily precipitation and temperature | CONUS, 1981 on | 800 m | daily | `rasters/prism` |
| [`hrrr`](#hrrr) | NOAA HRRR analysis | CONUS, 2014-10 on | 3 km | hourly | `rasters/hrrr` |
| [`era5`](#era5) | ECMWF ERA5 reanalysis, single levels | global, 1940 on (a few months behind) | 0.25° | hourly | `rasters/era5` |
| [`stage4`](#stage4) | NCEP Stage IV precipitation | CONUS, 2002 on; Alaska, 2020-07 on | 4.8 km | daily, 6-hourly, hourly | `stage4/<accumulation>` |
| [`nlcd`](#nlcd) | Annual NLCD land cover | CONUS, 1985–2024 | 30 m | static | `rasters/nlcd` |
| [`dem`](#dem) | elevation | global (`cop30`), US (`usgs10m`) | 30 m, 10 m | static | `dem/<product>` |
| [`soil`](#soil) | soil properties | CONUS | 100 m, 30 m | static | `soil/<product>` |
| [`burn-severity`](#burn-severity) | BAER soil burn severity | CONUS, fires 2017 on | 20–30 m | static | `rasters/burn-severity` |
| [`sentinel1-rtc`](#sentinel1-rtc) | Sentinel-1 radar backscatter | global, 2014-10 on | 10 m | each pass | `rasters/sentinel1-rtc` |
| [`opera-rtc`](#opera-rtc) | OPERA Sentinel-1 radar backscatter (NASA Earthdata) | land north of 60° S, 2016-04 on | 30 m | each pass | `rasters/opera-rtc` |
| [`nisar-gcov`](#nisar-gcov) | NISAR L-band radar backscatter (NASA Earthdata) | global, 2025-11 on | 10–20 m | each pass | `rasters/nisar-gcov` |
| [`sentinel2`](#sentinel2) | Sentinel-2 surface reflectance | global, 2015-06 on | 10 m | each pass | `rasters/sentinel2` |
| [`seismic`](#seismic) | seismic waveforms (EarthScope FDSN) | global | stations | native rate, raw counts | `waveforms/seismic/<rate>hz` |
| [`infrasound`](#infrasound) | infrasound waveforms (EarthScope FDSN) | global | sensors | native rate, raw counts | `waveforms/infrasound/<rate>hz` |

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
| `gamma0_hh`, `gamma0_hv`, `gamma0_vv`, `gamma0_vh` (linear power) | `sentinel1-rtc` and `opera-rtc` `hh`, `hv`, `vv`, `vh`; `nisar-gcov` `HHHH`, `HVHV`, `VVVV`, `VHVH` |

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

(nrcs)=
## NRCS SNOTEL and SCAN

Hourly station reports from the NRCS [Air and Water Database](https://wcc.sc.egov.usda.gov/awdbRestApi/swagger-ui/index.html):
SNOTEL in the western mountains and Alaska, and SCAN across the US. No token needed.

- `--variables` takes AWDB element codes: `SMS` (default, volumetric soil moisture),
  `STO` (soil temperature), `WTEQ` (snow water equivalent), `SNWD` (snow depth),
  `PREC` (precipitation since October 1), `TOBS` (air temperature), ….
- Sensors below ground are named for their depth in cm: `SMS_5cm`, `SMS_10cm`,
  `SMS_20cm`, `SMS_51cm`, `SMS_102cm` (2, 4, 8, 20, 40 in). The `depth` attribute has
  the exact value.
- Units are SI: soil moisture in m³ m⁻³ (AWDB's percent / 100), mm, °C. Times are UTC
  (AWDB reports local standard time). Data NRCS marks suspect are left out.
- Closed stations are included when their record overlaps the window.
- Every report is kept, as a CF ragged array, like `synoptic`.
- License: U.S. public domain.

Other soil moisture sources, and why NRCS rather than Synoptic, are in the developer
notes ([soil-moisture.md](https://github.com/gaia-hazlab/gaia-cli/blob/main/dev/soil-moisture.md)).

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
- `orography` (m) is the height of the model's terrain: smoothed, roughly each
  cell's mean. Every ERA5 surface value (`t2m`, `sf`, `sd`, ...) is valid at that
  height, not on the ground below it, and regridding onto a finer grid adds no
  terrain. In the mountains it is far below the peaks (around Peters Dome
  1,200–2,100 m against ~2,760 m at the source), so there ERA5 is too warm, puts
  too much of its precipitation down as rain, and spreads orographic precipitation
  over the whole cell. The usual first correction is a lapse rate:
  `t2m − 6.5 °C/km × (elevation − orography)`, with `elevation` from a DEM. Real
  lapse rates vary (inversions in valleys, Alaskan winters), so treat it as an
  estimate. On a map, `elevation − orography` says more than `orography` alone.
- ERA5 is a model run constrained by observations, not an interpolation of stations
  like PRISM. It assimilates surface stations (pressure, and 2 m temperature and
  humidity in its land-surface analysis) and, over the US from 2009, radar–gauge
  precipitation, but its `tp` is a short forecast and its values are cell means.
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

(opera-rtc)=
## OPERA RTC

Sentinel-1 backscatter processed by JPL's
[OPERA](https://hyp3-docs.asf.alaska.edu/guides/opera_rtc_product_guide/) project,
one burst at a time, from ASF DAAC: 30 m, in the bursts' UTM zone. It's an
alternative to `sentinel1-rtc`, with the same options, coordinates and stage-file
windows. Needs a free [Earthdata Login](https://urs.earthdata.nasa.gov) account: set
`EARTHDATA_TOKEN` to a
[user token](https://urs.earthdata.nasa.gov/documentation/for_users/user_token),
which expires after 60 days.

- `vv` and `vh` (default; `--polarizations` also takes `hh`, `hv`): gamma0 as
  **linear power**.
- `mask` (default): the layover and shadow mask of the pass's track, from OPERA's
  static layers, as CF flags (`valid`, `shadow`, `layover`, `layover_and_shadow`,
  `no_data`). `--static` also takes `local_incidence_angle` (from the terrain normal)
  and `incidence_angle` (from the ellipsoid normal), in degrees, `number_of_looks`,
  `rtc_anf_gamma0_to_beta0` and `rtc_anf_gamma0_to_sigma0`. Each is a variable along
  `time`: every pass gets its track's layers, read only for the tracks of the chosen
  passes.
- One time step per pass, its bursts mosaicked, with the same `platform`,
  `orbit_state`, `relative_orbit` and `absolute_orbit` coordinates as
  `sentinel1-rtc`, plus `burst_ids`.
- From 2016-04-14 (OPERA's own processing from 2022), over land north of 60° S.
- License: CC0-1.0 (NASA Earthdata data use policy).

```bash
gaia stage opera-rtc -i aoi.geojson -s 2023-09-01 -e 2023-09-30 -o out --relative-orbit 131
gaia stage opera-rtc -i aoi.geojson -s 2023-09-01 -e 2023-09-30 -o out --static mask --static local_incidence_angle
```

(nisar-gcov)=
## NISAR GCOV

L-band (24 cm) backscatter from [NISAR](https://nisar-docs.asf.alaska.edu/gcov/)'s
geocoded polarimetric covariance product (GCOV), radiometrically terrain corrected,
from ASF DAAC: 10 m in 40 MHz mode (20 m in 20 MHz mode), in the frames' UTM zone.
Needs `EARTHDATA_TOKEN`, as `opera-rtc` does.

- `HHHH` and `HVHV` (default; `--terms` also takes `VVVV`, `VHVH` where quad-pol): the
  diagonal covariance terms of frequency A, gamma0 as **linear power**, like
  `sentinel1-rtc`'s `hh` and `hv`. Most land is imaged HH+HV, where Sentinel-1 is
  usually VV+VH: compare like with like through the common names above.
- `mask` (default): 0 where a sample isn't valid (partly focused), 255 outside the
  swath, else the subswath number. `--layers` also takes `numberOfLooks` and
  `rtcGammaToSigmaFactor` (gamma0 × factor = sigma0).
- `--maturity validated` (default) reads `NISAR_L2_GCOV_V1`, acquired from
  2026-10-21; `--maturity provisional` reads `NISAR_L2_GCOV_PROVISIONAL_V1`, mostly
  2026-06-17 to 2026-10-23 (to be superseded by reprocessing in 2027). A store never
  mixes the two.
- One time step per pass, its frames mosaicked, with the `platform`, `orbit_state`,
  `relative_orbit` and `absolute_orbit` coordinates, plus `frames` and
  `range_bandwidth`.
- Each frame's grids are read remotely from its HDF5 file (about 7 GB), in 4
  processes: a few seconds per pass.
- License: CC0-1.0 (NASA Earthdata data use policy).

```bash
gaia stage nisar-gcov -i aoi.geojson -s 2026-10-21 -e 2026-11-08 -o out
gaia stage nisar-gcov -i aoi.geojson -s 2026-07-01 -e 2026-10-20 -o out --maturity provisional --relative-orbit 172
```

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
(seismic)=
## Seismic

Seismic waveforms from [EarthScope](https://www.earthscope.org)'s FDSN web services,
for the stations inside the AOI (use a circle, e.g. `radius_km: 250`, to limit the
distance).

- `--channel` (default `HH?,BH?`) lists channel codes in order of preference: each
  station keeps every channel of the first band it has, at its lowest location code.
  A station with HH and BH channels gives HHZ, HHN and HHE only. `--network`,
  `--station` and `--location` take FDSN codes with `?` and `*` wildcards.
- Stored as **raw counts at each channel's native rate**, one group per rate
  (`waveforms/seismic/50hz`, `.../100hz`), next to two tables: `sites` (one row per
  station position, with `distance_km` and `back_azimuth` from the stage file's
  `origin`, else the AOI's centroid) and `channels` (one row per channel requested,
  with its orientation, sensitivity and the group holding its data). Gaps are NaN on
  read.
- Samples snap to a grid starting at the window's start; `time_offset` gives each
  channel's true sample times, `time + time_offset`.
- With `node = tree["waveforms/seismic"]`, `gaia_cli.waveforms.to_physical(node,
  "50hz")` gives velocity (m/s) as counts / sensitivity, and
  `gaia_cli.waveforms.to_grid(node, rate=20)` puts every channel on one common grid.
- The window is read exactly, without the one-day padding of other datasets.
- License: set by each seismic network; `waveforms/seismic` lists its `networks`.

(infrasound)=
## Infrasound

Infrasound (pressure) waveforms from the same EarthScope services, stored like
[`seismic`](#seismic), in `waveforms/infrasound/<rate>hz` with `sites` and `channels`
tables.

- `--channel` (default `BDF,HDF`) lists channel codes in order of preference, chosen
  **per sensor position** rather than per station, so every element of an array is
  kept, whether its elements are stations (IM.I53H1 to H8) or location codes (AV.KENI
  01 to 06). Where one sensor has several streams (AK `EP` at 40 Hz and `20` at 20 Hz),
  the fastest is kept.
- `gaia_cli.waveforms.to_physical(node, "40hz")` gives pressure (Pa) as counts /
  sensitivity, the same as obspy's `remove_sensitivity`.
- License: set by each network; `waveforms/infrasound` lists its `networks`.
