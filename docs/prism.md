# PRISM daily precipitation

Module: `src/gaia_cli/prism.py`, `src/gaia_cli/prism_cache/` · Commands:
`gaia stage prism`, `gaia prism fill`

## Source

- Product: PRISM **AN81d** daily time series ("AN" = best-estimate analysis
  that uses all station networks), **800 m / 30 arc-second**, conterminous US.
  Precipitation is the latest version, **D2**.
- Access: unzipped COGs in the gaia PRISM cache (see [below](#the-gaia-prism-cache)),
  loaded with `odc.stac` from STAC items built in memory. File names look
  like `prism_ppt_us_30s_YYYYMMDD`.
- Variables: `-v ppt tmin tmax tmean tdmean vpdmin vpdmax` (default `ppt`,
  written as `precipitation`).
- Reference: [PRISM Datasets (PDF)](https://prism.oregonstate.edu/documents/PRISM_datasets.pdf).

## The gaia PRISM cache

PRISM only serves each grid as a `.zip` package, and blocks a file after it
is downloaded twice in 24 hours (repeated excess can get an IP blocked; see
the [web service docs](https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf)).
So grids are fetched from PRISM once and the COG inside is copied
byte-for-byte to an S3 cache that mirrors PRISM's layout:

```
s3://gaia/public/prism/collection.json                 STAC Collection (root)
s3://gaia/public/prism/us/an/800m/{var}/daily/{yyyy}/catalog.json    STAC Catalog per variable-year
                                                     prism_{var}_us_30s_{yyyymmdd}.tif
                                                     prism_{var}_us_30s_{yyyymmdd}.json   STAC Item
```

Browse it with any STAC client from
<https://s3.kopah.uw.edu/gaia/public/prism/collection.json>.

- **Reading** is anonymous over HTTPS
  (`https://s3.kopah.uw.edu/gaia/public/prism/...`); no credentials needed.
- **Filling** happens on demand: `gaia stage prism` downloads any missing
  grids first, which needs `AWS_PROFILE` with write access to `s3://gaia`.
  Without it, missing grids raise `CacheMiss`. Pre-fill a range with
  `gaia prism fill -s 2025-12-01 -e 2025-12-31 -v ppt tmean`.
- **STAC items:** each `.json` is a STAC 1.1 Item for that grid (projection
  and raster extensions). It records `prism:release`, `updated` (when the
  grid was fetched into the cache), the rest of the zip's `info.txt` as
  `prism:*` properties, the PRISM day as `start_datetime`/`end_datetime`, and
  a `via` link to the upstream zip. For loading, one day's per-variable items
  are merged into one item with an asset per variable.
- **Catalog:** every fill also adds the item to its variable-year
  `catalog.json`, and grows `collection.json`'s extent and summaries. These
  updates aren't locked, so concurrent writers could drop a link;
  `gaia prism rebuild-catalog` regenerates everything from the items.
- **Versions:** daily grids are re-released up to 8 times over ~6 months
  (release 8 is final). For cached grids below release 8,
  PRISM's release-date service (which doesn't count against the download
  limit) is queried, and newer releases are re-fetched, but never within 24 h
  of the previous fetch. Staged variables carry a `prism_release` attribute,
  one value per time step.
- Override the location with `PRISM_CACHE_URL`. Set `PRISM_CACHE_PROFILE` to
  use an AWS profile for the cache that differs from `AWS_PROFILE` (used for
  staged outputs); its `endpoint_url` and `region` come from `~/.aws/config`.
  `PRISM_CACHE_ENDPOINT` overrides the endpoint. With none set the cache is
  Kopah, and an `AWS_PROFILE` for another store leaves it read-only.
- Not done yet: no lock stops two writers fetching the same grid at once
  (TODO in `cache.py`), and nothing keeps recent days current
  automatically (a scheduled `gaia prism fill` could do this).

## Grid

| | |
| --- | --- |
| CRS | EPSG:4269 (NAD83, geographic) |
| Resolution | 0.008333…° (30″), roughly 0.6 × 0.9 km at 48°N |
| Dims after staging | `time`, `y` (latitude), `x` (longitude) |
| Staged on | native grid, clipped to the AOI with `all_touched=True` |

## Time: the "PRISM day"

From the PRISM Datasets document:

> A PRISM day is defined as 1200 UTC-1200 UTC (e.g., 7 AM-7AM EST) […] The
> dataset uses a day-ending naming convention, e.g., a day ending at 1200 UTC
> on 1 January is labeled 1 January.

So the value labelled **D** is the precipitation total over
**(D−1 12:00 UTC, D 12:00 UTC]**. In Pacific time that is 04:00 PST on D−1
to 04:00 PST on D.

- The STAC items (and therefore `time` after `odc.stac.load`) are stamped
  **D 00:00**. The timestamp is only a date label: it is neither the start
  nor the end of the accumulation period.
- `prism.load` writes CF bounds so the period is explicit:
  `time_bounds = (time − 12 h, time + 12 h]` and
  `precipitation.attrs["cell_methods"] = "time: sum"`.
- Google Earth Engine's copy (`OREGONSTATE/PRISM/AN81d`) stamps assets at
  noon UTC instead. Don't mix timestamps from the two sources without
  converting.

## Values

- Variable `ppt` → renamed `precipitation`, units **mm** per PRISM day,
  float32.
- No-data (−9999) is converted to NaN on load.
- In the western US, daily D2 grids are forced to sum to the monthly AN
  grids (M3). Recent days can change until the monthly product is final; see
  the [PRISM update schedule](https://prism.oregonstate.edu/documents/PRISM_update_schedule.pdf).

## Using with HRRR

HRRR precipitation is hourly and hour-ending (see [hrrr.md](hrrr.md)), so
PRISM day D corresponds to the 24 HRRR values stamped
**D−1 13:00 … D 12:00 UTC**. `gaia stage all` therefore starts its HRRR read
11 h before `--start`, and `gaia align` sums HRRR into PRISM days using the
two products' `time_bounds`.
