# Stage IV precipitation

Module: `src/gaia_cli/datasets/stage4.py` · Command: `gaia stage stage4`
(stage-file key and group `stage4/daily`)

## Source

- Product: NCEP/EMC **Stage IV** quantitative precipitation estimates, a mosaic of
  the River Forecast Centers' radar and gauge analyses, about 4.8 km.
- Archive: NCAR GDEX dataset [d507005](https://gdex.ucar.edu/datasets/d507005/)
  ("NCEP/CPC Four Kilometer Precipitation Set, Gauge and Radar",
  [10.5065/D69Z93M3](https://doi.org/10.5065/D69Z93M3)), CC-BY-4.0, read anonymously.
- Only the **24 h** totals are staged so far (`--accumulation daily`). The archive
  also holds 1 h and 6 h files; see [Hourly and 6 h](#hourly-and-6-h).

## Archive layout

```
https://data.gdex.ucar.edu/d507005/stage4/stage4.YYYYMM.tar    one month, uncompressed
└─ ST4.YYYYMMDD                                                 one tar per day, not in date order
   ├─ 01h files, hour-ending 00–23 UTC of that date
   ├─ 06h files ending 00, 06, 12, 18 UTC
   └─ 24h file ending 12 UTC                                    ← staged
```

The name of the 24 h member, per domain and era. The code matches these names; the
dates are approximate:

| period | CONUS | Alaska | format |
|---|---|---|---|
| 2002 to ~2014 | `ST4.YYYYMMDD12.24h.Z` | — | GRIB1, Unix compress |
| ~2015 to 2020-07 | `ST4.YYYYMMDD12.24h.gz` | `st4_ak.….24h.gz` (from 2020-07-01) | GRIB1, gzip |
| ~2020-07-20 on | `st4_conus.….24h.grb2` | `st4_ak.….24h.grb2` | GRIB2 |

From 2020 the daily tars also hold Puerto Rico (`st4_pr`) grids and `.gif` maps, which
are ignored. The 2019–20 tars use PAX headers. Earlier daily tars contain Alaska
`.gif` maps but no Alaska GRIB.

## Fetching and the cache

GRIB has to be read from a local file, but whole monthly tars are large (about
300 MB from 2020). GDEX serves byte ranges, so `fetch_day`:

1. reads only the monthly tar's 512-byte headers (about 31 range requests) to find
   each `ST4.YYYYMMDD`'s offset, and caches that index as JSON;
2. range-gets that one daily tar (about 1 MB in 2002, about 10 MB from 2020);
3. extracts the 24 h member for the domain, decompresses it (`.gz` with `gzip`,
   `.Z` with `unlzw3`) and writes it to the cache.

The cache is `~/.cache/gaia/stage4/` (`STAGE4_CACHE` to change it):
`YYYY/<grib file>` and `index/stage4.YYYYMM.json`. Cached grids are never
re-fetched; Stage IV isn't revised after it is archived. Days missing from the archive
(including a month GDEX hasn't posted yet) are NaN, with a warning. A window with no
data at all is an error.

## Grids

GRIB is read with GDAL's GRIB driver (`rioxarray.open_rasterio`), which needs the
`libgdal-grib` conda package. pip installs of rasterio include the driver. GDAL reads the
CRS from the file:

| | CONUS (HRAP) | Alaska |
|---|---|---|
| CRS | polar stereographic, sphere R = 6371200 m, lat_ts 60°N, lon_0 −105° | same |
| size (x × y) | 1121 × 881 | 460 × 530 |
| spacing | 4762.5 m | 4762.5 m |
| first cell center | x −1902531.67, y −3426604.93 | x −3288463.39, y −1007361.14 |
| outline (EPSG:4326) | −134.09 to −59.91°E, 19.78 to 57.87°N | −178.02 to −122.31°E, 45.81 to 75.72°N |

GRIB1 files store the spacing in whole metres: 4762 in 2002, 4763 in 2020. Their
coordinates as read therefore drift by up to about 800 m from the GRIB2 ones, although
it is the same grid. Every file is given the GRIB2 coordinates above (`stage4.GRIDS`),
so days from different eras stack. A file with any other shape is an error.

**Domain.** The CONUS and Alaska outlines overlap as lon/lat boxes (both cover
western Washington), so the domain is chosen in the grids' own CRS: the first grid
whose rectangle contains the AOI's bounds. An AOI in neither, such as Hawaii or
one spanning both, is an error. `EXTENT` lists both outlines, so the runner's
coarse check already skips AOIs outside both.

## Time

The 24 h total in `….YYYYMMDD12.24h` covers **(D−1 12:00, D 12:00] UTC**, the same
period as a [PRISM day](prism.md#time-the-prism-day). It is stamped like PRISM, D 00:00,
with `time_bounds` (time − 12 h, time + 12 h] and `cell_methods: time: sum`. Stage IV
and PRISM therefore share time coordinates, and `gaia align` can sum HRRR hours into
them.

The time comes from the file name. Inside the GRIB the reference time is the *start*
of the period (D−1 12 UTC), so cfgrib's `time` is a day earlier than the label.

**Alaska** grids start 2020-07-01. A request whose window starts earlier is an error;
the extra day the runner reads before the window (2020-06-30 for a window starting
2020-07-01) is NaN.

## Values

| variable | units | `cell_methods` | `standard_name` |
|---|---|---|---|
| `precipitation` | mm (from kg m⁻²) | `time: sum` | `lwe_thickness_of_precipitation_amount` |

- float32. No-data (9999) is NaN; about 40% of the CONUS grid, outside the RFCs'
  coverage, is always NaN.
- Checked against the archive copies under `/home/orca/data/in_situ/NCEP_stageIV_QPE`
  (2016-01-01, 2025-12-09): identical. For the Bolt Creek AOI, 2025-12-08 to 12-10,
  the daily means were 27.7, 114.1 and 60.4 mm, against PRISM's 30.5, 114.4 and 47.0 mm
  after `gaia align`.

## Hourly and 6 h

To add them, add values to `Accumulation`, extend `member_names` with the `01h`/`06h`
stamps, and use bounds (−1 h, 0) or (−6 h, 0) with the end time as the label, like
HRRR. Alaska has no 1 h product, so hourly with an Alaska AOI must be an error. In the
western RFCs (including NWRFC), 1 h Stage IV can be 6 h totals split by hour, so check
before using it for sub-daily timing.
