---
title: Datasets
---

| dataset | what | coverage | grid | time | group |
|---|---|---|---|---|---|
| [`synoptic`](#synoptic) | weather station reports | mostly US, Alaska included | stations | every report | `stations/synoptic` |
| [`prism`](#prism) | PRISM daily precipitation and temperature | CONUS, 1981 on | 800 m | daily | `rasters/prism` |
| [`hrrr`](#hrrr) | NOAA HRRR analysis | CONUS, 2014-10 on | 3 km | hourly | `rasters/hrrr` |
| [`nlcd`](#nlcd) | Annual NLCD land cover | CONUS, 1985–2024 | 30 m | static | `rasters/nlcd` |
| [`dem`](#dem) | elevation | global (`cop30`), US (`usgs10m`) | 30 m, 10 m | static | `dem/<product>` |
| [`soil`](#soil) | soil properties | CONUS | 100 m, 30 m | static | `soil/<product>` |
| [`burn_severity`](#burn-severity) | BAER soil burn severity | CONUS, fires 2017 on | 20–30 m | static | `rasters/burn_severity` |

A dataset that the AOI or window misses (for example `prism` in Alaska) is skipped
with a warning. Licenses differ: check each group's `license` attribute before sharing
a store.

(synoptic)=
## Synoptic stations

Station reports from the [Synoptic Data](https://synopticdata.com) API. Requires a
`SYNOPTIC_TOKEN`.

- `precip_accum` (default) is staged as `precipitation`, a running total in mm. Other
  Synoptic variables, such as `air_temp`, keep their Synoptic names and units.
- Every report is kept, as a CF ragged array. `gaia_cli.stations.to_grid(ds)` puts the
  reports on a regular hourly (or daily) grid, and `to_dataframe(ds)` gives one row per
  report.
- License: set by each station network.

(prism)=
## PRISM

[PRISM](https://prism.oregonstate.edu) AN81d daily grids, 800 m, EPSG:4269.

- Variables: `ppt` (default, staged as `precipitation`, mm), `tmin`, `tmax`, `tmean`
  (degC).
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

- Variable: `precipitation-surface` (default, staged as `precipitation` in mm per
  hour), or any other in `gaia stage hrrr --help`.
- Precipitation is **hour-ending**: the value at 13:00 fell between 12:00 and 13:00
  UTC.
- Wind components are grid-relative, not east/north.
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
