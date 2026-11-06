# ERA5 (plan, step 2.3)

Module: `src/gaia_cli/datasets/era5.py` · Command: `gaia stage era5`

Step 2.3 stages single-level fields only. Pressure-level profiles are a follow-up
(see [Open](#open)).

## What Peters Dome reads

Two things, both from Open-Meteo today ([04](roadmap/04-datasets.md#what-the-workflows-actually-use)):

- daily temperature and precipitation at the source (a time series at one point);
- atmosphere profiles (temperature and wind on pressure levels) for infrasound
  propagation.

ERA5-Land is never read and stays out.

## Source options

| | Earthmover ERA5 (Arraylake) | Google ARCO-ERA5 | Open-Meteo |
|---|---|---|---|
| Format | Icechunk, read with xarray | Zarr on GCS | JSON API, one point per request |
| Auth | free Earthmover account, `ARRAYLAKE_TOKEN` | none | none; non-commercial terms, rate limits |
| Grid | native 0.25°, hourly, 1940– | native 0.25°, hourly, 37 levels | resampled to points |
| Chunking | a temporal and a spatial copy | one hour per chunk, whole globe | n/a |
| Provenance | Icechunk snapshot id, as for `hrrr` | Zarr path | API URL |

**Decision:** Earthmover, repo `earthmover-public/era5`. It works like `hrrr`
(Icechunk session, snapshot id in provenance, clip on the native grid), its temporal
copy makes a small AOI over a short window cheap, and it has pressure levels for the
profiles.

The `gaia` org's marketplace subscription `gaia/era5-subscription` is a mirror of it
(same snapshot and commits). The code names `earthmover-public/era5`, which any
Earthmover account's token can read; the subscription works only for `gaia` org
tokens. `earthmover-public/era5-surface-aws` (in Earthmover's docs; subscribed as
`gaia/era5-surface-aws-subscription`) has only the single-level groups and isn't
used.

The free tier moves data in quarterly (commit "quarterly move single/temporal/status
2026-04-01..2026-06-30", 2026-10-12).

## The store

Checked 2026-11-05, snapshot `VKGPX45BPEQJMVZZXZEG`.

| group | contents | chunks (time, [level,] lat, lon) |
|---|---|---|
| `single/temporal` | 39 single-level variables | (8736, 12, 12) |
| `single/spatial` | the same | (1, 721, 1440) |
| `pressure/temporal` | `t`, `u`, `v`, `z`, `q`, `r`, `w`, `pv` on 13 levels, 1000–50 hPa | (8736, 1, 12, 12) |
| `pressure/spatial` | the same | not checked |
| `500hPa` | `spatial` and `temporal` subgroups, not checked | |

- Dims `valid_time`, `latitude` (90 → −90), `longitude` (0 → 359.75); static `lsm`,
  `z_sfc`, `sdor`, `slor` coordinates. CRS `EPSG:4326` (in `proj:code`).
- Time: hourly, 1940-01-01 to `time_coverage_end` (2026-06-30T23:00 at this
  snapshot). The summary says "updated daily", ERA5T (`expver=0005`) replaced by final
  ERA5 about two months later. The free tier is said to lag ~3 months; the store's
  `time_coverage_end` is the source of truth.
- Attributes: CF-1.7, `license` CC-BY-4.0, an `acknowledgement` (Copernicus wording,
  which CC-BY requires us to keep), `references` (Hersbach et al. 2020), source
  `reanalysis-era5-single-levels` from the CDS.
- Units as GRIB: `tp`, `cp`, `lsp`, `sf` in m (water equivalent), `GRIB_stepType`
  `accum`; radiation in J m⁻²; temperatures in K; `z` is geopotential (m² s⁻²), not
  height.
- Compression: `pcodec`, so the `pcodec` package is a dependency.

Timing for Peters Dome (a 2 × 2-cell box, from Seattle): open 1.3 s; `tp` and `t2m`
for 47 days 2.8 s; `t`, `u`, `v`, `z` on all 13 levels for 4 hours 17 s (52 one-year
chunks).

## Auth and the S3 endpoint

- `arraylake` reads `ARRAYLAKE_TOKEN` (an `ema_…` API key from Organization Settings →
  API Clients; the community tier allows one) and vends S3 credentials for the repo.
  Locally the token can live in `.env`; CI needs the repository secret.
- The bucket config has no `endpoint_url`, so icechunk's AWS SDK uses the endpoint
  from the active AWS profile. With a profile for another store (Kopah), the read
  fails with a DNS error. `AWS_ENDPOINT_URL_S3=https://s3.us-east-1.amazonaws.com`
  fixes it; `AWS_IGNORE_CONFIGURED_ENDPOINT_URLS=true` does not. `era5.py` must set
  this only while opening and reading, so a Kopah `-o` still writes to Kopah
  (`hrrr.py` pins its endpoint for the same reason, but `arraylake` builds the
  storage, and `StorageOptions` has no endpoint).
- `arraylake >=1.8` requires `python <3.15`, so the locked Python is 3.14 (3.13 on
  linux-64).

## What `era5.py` does

- Reads `single/temporal`. `EXTENT` is global from 1940 with an open end (written by
  hand); a window that starts after the store's last hour fails with that date, and
  one that ends after it logs a warning.
- Renames `valid_time` → `time`, `latitude`/`longitude` → `y`/`x`, shifts longitudes
  to −180–180 (Alaska is ~209°E), writes EPSG:4326, and clips with `io.clip_to_aoi`
  (a point AOI keeps the touched 0.25° cells). The static `z_sfc` (surface
  geopotential) is staged as `orography`, ERA5's terrain height in m (÷ 9.80665),
  with no time dimension and `bilinear` regridding. ERA5's terrain is smooth: around
  Peters Dome the cells are 1,200–2,100 m against ~2,760 m at the source, so `t2m`
  is valid well below the summit; `orography` gives the height difference for a
  lapse-rate correction. `gaia align` drops it, as it does any variable without
  time in a timed group. `lsm`, `sdor` and `slor` are dropped.
- `variables` defaults to `tp`, `sf`, `t2m`, `swvl1`, `swvl2`, `sd`, and takes any
  of the 39 (see [Variables](#variables)). Names
  stay ECMWF's; giving one quantity one name across datasets is for `gaia align`, in
  its own PR.
- Units: K → degC; `tp`, `cp`, `lsp`, `sf`, `sd` (m of water) → mm. `GRIB_*`
  attributes are dropped.
- Time: every value is hour-ending, so the group has `time_bounds` `(T − 1 h, T]`;
  `cell_methods` follows `GRIB_stepType`: `accum` → `time: sum`, `max` →
  `time: maximum`, `instant` → `time: point` (align averages these).
- `REGRID_METHOD = "conservative"`, like `hrrr`, so precipitation totals survive
  `gaia align`.
- Provenance: CC-BY-4.0, ECMWF and Earthmover as providers, Hersbach et al. (2020) as
  `sci:doi`/`sci:citation`, the store's Copernicus `acknowledgement`, the snapshot id
  and a `via` link to the repo.

Peters Dome (15 km radius, 18 days, `tp`, `t2m`, `swvl1`, `sd`): 11 s, 51 kB. `sd` is 4–10 m
of water equivalent in September: ERA5 models glacier cells as a deep snowpack, so
`sd` there says nothing about seasonal snow.

## Variables

The defaults compare with PRISM and HRRR (`tp`, `t2m`), split rain from snow (`sf`),
and give soil water and snow (`swvl1`, `swvl2`, `sd`). Each variable adds about half a
second for a small AOI. Others in the store worth adding when a workflow needs them:

| purpose | variables | why |
|---|---|---|
| Rain vs snow (default) | `sf` (snowfall, mm) | `tp − sf` is rain; matters for rain-on-snow and debris-flow triggering |
| Comparing with HRRR | `d2m` (dewpoint), `u10`/`v10` (10 m wind), `sp` (surface pressure) | HRRR has the same quantities; ERA5 winds are east/north, HRRR's are grid-relative |
| Snowmelt energy | `ssrd` (downward solar radiation) | Drives melt along with `t2m` |
| Frozen ground | `stl1` (soil temperature, 0–7 cm) | Frozen soil changes infiltration; relevant in Alaska |
| Storm type | `cp`/`lsp` (convective, large-scale precipitation) | Short intense bursts vs long frontal rain |
| Atmospheric rivers | `tcwv` (total column water vapour) | A rough indicator |

Not in this store, though ERA5 has them:

- `mx2t`/`mn2t` (hourly maximum and minimum 2 m temperature). `gaia align` can take
  the daily minimum and maximum of hourly `t2m`, the closest match to PRISM's `tmin`
  and `tmax` (slightly narrower than true extremes).
- `swvl3`/`swvl4` (soil water at 28–100 and 100–289 cm).
- `ro`, `smlt`, `e` (runoff, snowmelt, evaporation) and `rsn` (snow density, needed to
  turn `sd` into a physical depth).
- Vertically integrated moisture flux (IVT), the usual atmospheric-river measure.

## Open

1. **Profiles (follow-up).** `pressure/temporal` tops out at 50 hPa (~20 km).
   Infrasound ducts in the stratosphere (~40–50 km) need more. Check how the Peters
   Dome analysis uses its profiles before choosing between these and another source
   (ARCO-ERA5 has 37 levels to 1 hPa).
2. Variable names: **decided** to keep ECMWF's short names (`tp`, `t2m`, …) and
   convert only units. Whether to give the same quantity one name across datasets
   (`tp` here, `precipitation` in `hrrr` and `prism`) is planned for `gaia align`, in
   its own PR.
3. The hour-ending convention for `accum` variables is the CDS's; check it against a
   CDS download for one hour.
