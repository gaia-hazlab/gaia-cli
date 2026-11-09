# Stage IV precipitation

Module: `src/gaia_cli/datasets/stage4.py` · Command: `gaia stage stage4 -a <accumulation>`
(stage-file keys and groups `stage4/daily`, `stage4/sixhourly`, `stage4/hourly`)

## Source

- Product: NCEP/EMC **Stage IV** quantitative precipitation estimates, a mosaic of
  the River Forecast Centers' radar and gauge analyses, about 4.8 km.
- Archive: NCAR GDEX dataset [d507005](https://gdex.ucar.edu/datasets/d507005/)
  ("NCEP/CPC Four Kilometer Precipitation Set, Gauge and Radar",
  [10.5065/D69Z93M3](https://doi.org/10.5065/D69Z93M3)), CC-BY-4.0, read anonymously.
- Three accumulations, each its own group: `daily` (24 h, the default),
  `sixhourly` (6 h) and `hourly` (1 h), from the archive's `24h`, `06h` and `01h`
  files. Alaska has no hourly grids; see [Hourly and 6 h](#hourly-and-6-h).

## Archive layout

```
https://data.gdex.ucar.edu/d507005/stage4/stage4.YYYYMM.tar    one month, uncompressed
└─ ST4.YYYYMMDD                                                 one tar per day, not in date order
   ├─ 01h files, hour-ending 00–23 UTC of that date            ← hourly
   ├─ 06h files ending 00, 06, 12, 18 UTC                       ← sixhourly
   └─ 24h file ending 12 UTC                                    ← daily
```

A grid ending at 00 UTC on D is in D's tar, not D−1's. The members' names, per domain
and era (shown for 24 h; `01h` and `06h` follow the same pattern). The code matches
these names; the dates are approximate:

| period | CONUS | Alaska | format |
|---|---|---|---|
| 2002 to ~2014 | `ST4.YYYYMMDD12.24h.Z` | — | GRIB1, Unix compress |
| ~2015 to 2020-07 | `ST4.YYYYMMDD12.24h.gz` | `st4_ak.….24h.gz` (from 2020-07-01) | GRIB1, gzip |
| ~2020-07-20 on | `st4_conus.….24h.grb2` | `st4_ak.….24h.grb2` | GRIB2 |

From 2020 the daily tars also hold Puerto Rico (`st4_pr`) grids and `.gif` maps, which
are ignored. The 2019–20 tars use PAX headers. Earlier daily tars contain Alaska
`.gif` maps but no Alaska GRIB. The 2019-01 to 2019-03 monthly tars (136–164 MB, against
about 90 MB for their neighbours) also hold every grid as a loose file after its daily
tar; indexing reads their many small headers 1 MB at a time, so it takes seconds
rather than thousands of requests.

## Fetching and the cache

GRIB has to be read from a local file, but whole monthly tars are large (about
300 MB from 2020). GDEX serves byte ranges, so `fetch_day`:

1. reads only the monthly tar's 512-byte headers (about 31 range requests) to find
   each `ST4.YYYYMMDD`'s offset, stopping once every day is found, and caches that
   index as JSON;
2. range-gets that one daily tar (about 1 MB in 2002, about 10 MB from 2020);
3. extracts the day's members for the domain and accumulation (one 24 h, four 6 h or
   24 1 h grids), decompresses them (`.gz` with `gzip`, `.Z` with `unlzw3`) and writes
   them to the cache.

The cache is `./.cache/gaia/stage4/`, under the directory gaia runs in, not the
home directory, which is small on shared HPC systems. Set `STAGE4_CACHE` to change
it, for example to share one cache across working directories on scratch storage.
It holds `YYYY/<grib file>` and `index/stage4.YYYYMM.json`. Cached grids are never
re-fetched; Stage IV isn't revised after it is archived. If a day's tar lacks some of
its grids, the tar is fetched again on the next run. Grids missing from the archive
(including a month GDEX hasn't posted yet) are NaN, with one warning that counts them. A
window with no data at all is an error.

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

| accumulation | file | covers | stamped | `time_bounds` |
|---|---|---|---|---|
| `daily` | `….YYYYMMDD12.24h` | (D−1 12:00, D 12:00] UTC | D 00:00, like a [PRISM day](prism.md#time-the-prism-day) | (time − 12 h, time + 12 h] |
| `sixhourly` | `….YYYYMMDDHH.06h`, HH 00/06/12/18 | (T − 6 h, T] | T, the period end | (time − 6 h, time] |
| `hourly` | `….YYYYMMDDHH.01h` | (T − 1 h, T] | T, the period end, like HRRR | (time − 1 h, time] |

All have `cell_methods: time: sum`. Daily Stage IV and PRISM therefore share time
coordinates, hourly Stage IV and HRRR share theirs, and `gaia align` can sum hours into
days. For a window `-s S -e E`, the sub-daily steps are those whose end lies in [S, E]
(as for HRRR), so the first hourly step of a day covers 23:00–00:00 of the day before.
Checked live: for Bolt Creek on 2025-12-09, the 24 hourly and the four 6-hourly basin
means sum to the daily one (114.1 mm).

The time comes from the file name. Inside the GRIB the reference time is the *start*
of the period (D−1 12 UTC), so cfgrib's `time` is a day earlier than the label.

**Alaska** grids start 2020-07-01. A request whose window starts earlier is an error;
the extra day the runner reads before the window (2020-06-30 for a window starting
2020-07-01) is NaN. Alaska has 24 h and 6 h grids only, so `hourly` with an Alaska AOI
is an error.

## Values

| variable | units | `cell_methods` | `standard_name` |
|---|---|---|---|
| `apcp` (`precipitation` after `gaia align`) | mm (from kg m⁻²) | `time: sum` | `lwe_thickness_of_precipitation_amount` |

- float32. No-data (9999) is NaN; about 40% of the CONUS grid, outside the RFCs'
  coverage, is always NaN.
- Checked against the archive copies under `/home/orca/data/in_situ/NCEP_stageIV_QPE`
  (2016-01-01, 2025-12-09): identical. For the Bolt Creek AOI, 2025-12-08 to 12-10,
  the daily means were 27.7, 114.1 and 60.4 mm, against PRISM's 30.5, 114.4 and 47.0 mm
  after `gaia align`.

## Hourly and 6 h

**Hourly gaps in the West.** Before 2018 the hourly mosaic has no data over parts of
the western US, so those grids exist but are NaN there. Sampled at 12 UTC on the first
day of a month (fraction of cells with data):

| | WA (NWRFC) | N. California (CNRFC) | CO, UT, AZ, MT | OK (ABRFC) |
|---|---|---|---|---|
| 1 h, 2002-07 to 2017-12 | none | 2011-01: yes; 2017-01: none | yes | yes |
| 1 h, from 2018-01 | yes | yes | yes | yes |
| 6 h | yes, except 2002-07 | yes | yes | yes |

Washington's first hourly data in a binary search over months was 2018-01 (none in
2017-12). Steps whose grid has no data anywhere over the AOI are left NaN, and `load()`
logs one warning that counts them and, for hourly, names this gap. Daily and 6-hourly
cover these regions. Where hourly exists in the western RFCs (including NWRFC), it can
be 6 h totals split by hour, so check before relying on it for sub-daily timing.
