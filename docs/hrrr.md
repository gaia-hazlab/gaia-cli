# HRRR analysis

Module: `src/gaia_cli/hrrr.py` · Command: `gaia stage hrrr`

## Source

- Product: NOAA High-Resolution Rapid Refresh (HRRR) **analysis**, built by
  dynamical.org from NOAA's open data archive (dataset `noaa-hrrr-analysis`,
  version `0.2.0`, CC-BY-4.0).
- Access: anonymous Icechunk repository at
  `s3://dynamical-noaa-hrrr/noaa-hrrr-analysis/v0.2.0.icechunk/` (us-west-2).
- Reference: [dynamical.org catalog entry](https://dynamical.org/catalog/noaa-hrrr-analysis/).

## Grid

| | |
| --- | --- |
| CRS | Lambert Conformal Conic (2SP), lat_0 = lat_1 = lat_2 = 38.5°, lon_0 = −97.5°, spherical earth R = 6 371 229 m |
| Resolution | 3 km |
| Coordinates | 1D `x`/`y` in metres, plus 2D `latitude`/`longitude` auxiliary coords |
| Staged on | native grid, clipped to the AOI (reprojected into the LCC CRS) with `all_touched=True` |

The CRS has no EPSG code. It comes from the GRIB WKT stored in `spatial_ref`.

## Time

- Hourly, from 2014-10-01 00:00 UTC to the present.
- dynamical.org builds the analysis by "concatenating the first step of each
  forecast along the time dimension. Accumulated variables (e.g.
  precipitation) are read from the second step of the previous hour's
  forecast."

What a value at `time = T` means depends on the variable's `step_type`
attribute:

| `step_type` | Variables | Value at T | Staged bounds |
| --- | --- | --- | --- |
| `instant` | everything except the two below | valid at T, from the analysis of the run starting at T | none; `cell_methods = "time: point"` |
| `avg` | `precipitation_surface`, `snowfall_surface` | average over **(T − 1 h, T]**, from the +1 h step of the run starting at T − 1 h | `time_bounds = (T − 1 h, T]` |

In other words, **precipitation is hour-ending**: the value stamped
`2025-12-01T13:00` is the precipitation that fell from 12:00 to 13:00 UTC.

## Values

- `precipitation_surface`: "Average precipitation rate since the previous
  forecast step", in kg m⁻² s⁻¹ (≡ mm s⁻¹). `hrrr.load` multiplies it by
  3600 s to get **mm per hour**, renames it `precipitation`, and sets
  `cell_methods = "time: sum"`.
- `snowfall_surface`: average rate in m s⁻¹, staged unchanged with
  `cell_methods = "time: mean"`.
- `wind_u_*` / `wind_v_*` are **grid-relative**: velocity along the LCC grid's
  x/y axes, not eastward/northward. After regridding they must be rotated
  before use as true east/north components. Scalar resampling alone is wrong
  for them.
- `categorical_*_surface` are 0/1 flags. Averaging them gives a fraction,
  not a category.
- `geopotential_height_cloud_ceiling` and `percent_frozen_precipitation_surface`
  are NaN by design where there is no ceiling or no precipitation.

## Using with PRISM

To match PRISM day D, which covers (D−1 12:00, D 12:00] UTC (see
[prism.md](prism.md)), sum the 24 hour-ending HRRR values stamped
**D−1 13:00 … D 12:00**. `gaia stage` reads every timed dataset one day
beyond each end of its window ([stage-file.md](stage-file.md#windows)), so for
`2025-12-01/2025-12-20` HRRR covers `2025-11-30T00:00` to `2025-12-21T23:00` and
every PRISM day in the window is complete after `gaia align`. Only the padding
day before the window is incomplete, and `gaia align` sets it to NaN rather than
undercounting it.
