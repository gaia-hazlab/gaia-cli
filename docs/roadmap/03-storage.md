# 03 · Storage: one Zarr store, Zarr + Parquet, or Icechunk + Iceberg?

## What we store

The shape of the data should drive the format, so start by classifying it:

| kind | examples | shape | typical size per event |
|---|---|---|---|
| Gridded, regular time | HRRR, PRISM, ERA5, SMAP | `(time, y, x)`, dense | MB–GB |
| Gridded, irregular scenes | S1 RTC, HLS, NISAR, ArcticDEM strips | `(time, y, x)`, sparse in time, per-scene metadata | 100 MB–GB |
| Gridded, static | DEM, land cover, soils | `(y, x)` or `(depth, y, x)` | MB–100 MB |
| Waveforms | seismic, infrasound | `(channel, time)`, regular per channel; 1–100 Hz; rates differ across channels | 10–500 MB for minutes to hours |
| Station observations | Synoptic | irregular reports, many variables, gaps, QC flags; stations come and go | 10⁴–10⁶ rows |
| Vectors | glacier outlines, geology, landslide polygons, station metadata | features + attribute columns | KB–MB |
| Scene indexes | per-scene date, track, orbit, cloud %, AOI coverage | rows | KB |

A Peters Dome example: 100 seismic stations × 3 components × 40 min × 20 Hz is about
58 MB of float32. Thirty Synoptic stations hourly for 3 months is about 80 k values.
The tabular data is tiny. Its problems are **shape and fidelity**, not size.

## Today: everything in one Zarr DataTree

`stations/*` are xvec CF-encoded dense arrays on a uniform time grid. This forces lossy
decisions at stage time:

- **Synoptic** running totals are snapped onto an hourly grid. Since the fix in
  [../precip-timing.md](../precip-timing.md), each hour holds the latest report at or
  before it, and a `report_time` coordinate keeps the actual time. Sub-hourly reports
  are still thinned to one per hour, and other variables are dropped.
- **Seismic data** is low-passed and interpolated onto one 20 Hz grid across all
  stations, which discards HH 100 Hz content.
- **Vectors** need xvec CF encoding. Peters Dome hit a real bug: GlaThiDa `LON`/`LAT`
  attribute columns collided with CF `lon`/`lat` arrays on case-insensitive file
  systems. String columns need re-encoding. Every reader needs `xvec.decode_cf()`.
- **Scene metadata** (orbit, track, release) becomes 1-D coordinates on `time`. That
  works, but the scenes cannot be queried without opening the arrays.

The upside is real: **one URL, one `open_datatree()` call, everything aligned in xarray**.
`gaia align` and the planned `convert landlab` operate on it directly. Each store also
gets a STAC Item next to it, and `gaia catalog rebuild` collects those Items into a
catalog ([../stac-catalog.md](../stac-catalog.md)).

## Options

### A. Keep a single Zarr store as it is

Everything stays in the DataTree, including the dense hourly `(station, time)` grid.
This is today's behavior, with the lossy decisions listed above.

- **For:** no work, no new formats, no new dependencies.
- **Against:** it keeps every fidelity problem above.

### A′. One Zarr store with CF ragged arrays for observations (**recommended to start**)

Keep one Zarr DataTree per event, but store irregular observations losslessly, using the
CF Discrete Sampling Geometries (DSG) **contiguous ragged array** representation instead
of a dense grid.

```
s3://gaia/scratch/peters-dome/
  peters-dome.json         # STAC Item = the manifest (see "STAC Item as the manifest")
  stage.yaml               # the resolved stage file, as an asset of the Item
  peters-dome.zarr/        # DataTree, consolidated as today
    rasters/*              # (time, y, x) or (y, x), native grids
    waveforms/seismic/100hz, waveforms/seismic/40hz, ...   # (channel, time) at native rate
    stations/synoptic      # CF DSG timeSeries, contiguous ragged array (below)
    scenes/sentinel1       # 1-D (scene) arrays: id, time, orbit, track, coverage
    vectors/*              # xvec-encoded features, as today
```

`stations/synoptic` as a contiguous ragged array:

| dim | variables |
|---|---|
| `station` | `station_id`, `geometry` (xvec), `elevation`, `network`, `row_size` (reports per station; CF `sample_dimension: obs`) |
| `obs` | `time`, `value` for each variable, `report_type`, `qc` |

- `featureType: timeSeries` is set on the group. `ds.to_dataframe()` on the `obs`
  dimension gives the long table, and `np.repeat(station_id, row_size)` adds the station
  column. A `gaia_cli.stations.to_dataframe(ds)` helper does both.
- Dense hourly or daily grids become **derived views** through a helper,
  `gaia_cli.stations.to_grid(ds, freq="1h", how="sum")`. Nothing writes the dense
  grid to the store. Existing code and examples switch to the ragged layout plus the
  helper. `gaia align` already passes station groups through unchanged.
- Seismic channel metadata (NSLC, location, sample rate, sensitivity, distance, back
  azimuth) is a set of 1-D arrays on a `channel` dimension in each waveform group.
- **For:**
  - It is lossless for observations. That is the main reason for B, without B's
    costs.
  - It uses a CF standard rather than a gaia-specific layout.
  - There is still one store, one `open_datatree()` call and consolidated metadata, so
    public B2/Kopah reads work unchanged.
  - There are **no new dependencies** ([design goal 2](README.md#design-goals)).
  - The writer is small: sort by station and time, then compute `row_size`
    ([design goal 1](README.md#design-goals)).
- **Against:**
  - Readers must know the ragged layout. Tabular users call `to_dataframe()`; xarray
    users go through the `to_grid` helper.
  - Strings are stored as numpy fixed-width unicode (`<U…`), which zarr-python writes
    as a fixed-length UTF-32 data type. That round-trips in xarray (tested; see
    [Resolved questions](#resolved-questions)), but readers other than zarr-python may not
    support it yet.
  - There is no predicate pushdown. At per-event sizes (10⁴–10⁶ rows) the whole table
    fits in memory, so this does not matter yet.
  - Vectors keep the xvec CF-encoding issues listed above (column-name collisions,
    string re-encoding). Few vectors are on the priority list
    ([04](04-datasets.md#prioritized-new-datasets)), so this is tolerable for now.
  - Moving to Iceberg later (C) means rewriting the ragged arrays as Parquet. That is a
    straightforward `to_dataframe().to_parquet()`.

### B. Zarr for arrays, GeoParquet for tables (**later, if a trigger occurs**)

```
s3://gaia/scratch/peters-dome/
  peters-dome.json                  # STAC Item = the manifest
  stage.yaml
  peters-dome.zarr/                 # rasters/*, waveforms/*
  tables/
    synoptic_observations.parquet   # long format: station_id, time, variable, value, qc, ...
    synoptic_stations.parquet       # GeoParquet: station_id, geometry, elevation, network, ...
    seismic_channels.parquet        # GeoParquet: NSLC, geometry, sample rate, sensitivity, distance, baz
    sentinel1_scenes.parquet        # scene index: id, time, orbit, track, coverage, release
  vectors/
    landslide_outlines.parquet      # GeoParquet
```

The layout is the same as A′, with the ragged arrays, scene arrays and vectors moved out
of Zarr into (Geo)Parquet files. Each file is another asset of the STAC Item.

- **For:**
  - The best interoperability: GeoPandas, Polars, R, QGIS, and DuckDB when it is useful.
  - Vectors need no CF-geometry encoding.
  - It matches how landlab already thinks (index CSVs, raw rainfall CSV,
    `stations.gpkg`).
  - **It is forward compatible with Iceberg (option C).** Iceberg data files are
    Parquet, and `pyiceberg`'s `Table.add_files` registers existing Parquet files
    without rewriting them.
- **Against:**
  - **pyarrow is a new and large dependency.** It is not in the current lockfile, and
    GeoPandas needs it to read or write Parquet.
  - There are two formats to keep consistent, with no atomic commit across them. That
    is acceptable for single-writer staging.
  - There are more files per event, though the STAC Item lists them all, so anonymous
    B2 reads still need no LIST.
  - At per-event sizes the query engines that make Parquet attractive (DuckDB, Polars)
    are not needed. Tables load straight into (geo)pandas.

### C. Icechunk repo for arrays + Iceberg tables (Arraylake)

[Arraylake added Iceberg tables in public beta on 2026-09-23](https://www.earthmover.io/blog/iceberg-in-arraylake).
Iceberg namespaces and tables sit next to Icechunk repos in one catalog with shared auth
and storage, exposed via the Iceberg REST catalog to DuckDB, Spark, Polars, PyIceberg
and Trino. Their worked example is exactly ours: station observations in Iceberg joined
against GFS in Icechunk. Geometries are WKB plus bbox columns (the GeoParquet pattern).

- **For:**
  - Versioned, transactional arrays (Icechunk snapshots) and tables (Iceberg commits).
  - Concurrent writers.
  - **Cross-event queries**, e.g. one Synoptic table partitioned by event across every
    FlowAlert event we stage. This is where Iceberg pays off over loose Parquet.
  - gaia-cli already reads Icechunk (HRRR from dynamical.org).
- **Against:**
  - It is beta: no multi-table transactions, no garbage collection, no views, and no
    Azure storage. So there is still **no atomic commit across a table and an array**.
  - Iceberg needs a catalog service in the commit path. Arraylake is the hosted one
    (the Community tier is free), which adds an account and auth dependency for every
    user and a vendor dependency. It is unclear whether it can sit on top of
    Kopah/B2 buckets.
  - Public anonymous reads (the current B2/Kopah workflow) become a catalog question,
    not just a URL.
  - Its value is low **per event**. A single event bundle has one writer and a handful
    of small tables.

## STAC Item as the manifest

gaia-cli already writes a STAC Item next to every store, and `gaia catalog rebuild`
collects them into a `catalog.json` ([../stac-catalog.md](../stac-catalog.md)). The
Item already carries what a separate `manifest.json` would: provenance properties, the
time range, bbox, providers, licenses and links. So there is **no separate manifest**
under any option. The Item is the manifest, and each event becomes one Item with
several assets:

| asset key | href | type | notes |
|---|---|---|---|
| `data` | `./peters-dome.zarr` | `application/vnd+zarr` | As today, with `xarray:open_kwargs` |
| `config` | `./stage.yaml` | `application/yaml` | The **resolved** stage file: AOIs, windows and options after defaults, so the run can be repeated exactly. Written for flag-only runs too. Role `metadata`. |
| `<table>` (B only) | `./tables/<name>.parquet` | `application/vnd.apache.parquet` | With the STAC `table` extension (`table:columns`, `table:row_count`, `table:primary_geometry`) |

- The Item's asset list is the file list, so anonymous B2 reads, which cannot LIST,
  still find every file.
- `gaia catalog rebuild` already skips the inside of `*.zarr` directories, so it finds
  an Item at the top of the event folder unchanged. Only the Item builder
  (`provenance.stac_item`) needs to emit more than one asset.
- Per-dataset provenance stays in each group's attrs, as today. The Item merges it, as
  it does for `gaia stage all`.
- Moving from A′ to B adds Parquet assets to the same Item. The catalog and the event
  layout don't change.

## Recommendation

1. **Start with A′.** Keep one Zarr DataTree per event. Store station observations as
   CF DSG contiguous ragged arrays, waveforms at native rate with channel metadata on a
   `channel` dimension, and scene indexes as 1-D arrays. Derive dense station grids on
   demand. This fixes the fidelity problems with no new dependencies.
2. **The STAC Item is the manifest**, with the resolved `stage.yaml` as an asset.
3. **Move to B when either trigger occurs:**
   - (a) users regularly work on the tables outside xarray (DuckDB, R, QGIS, ML
     pipelines) and `to_dataframe()` becomes friction, or
   - (b) we want Iceberg (C), since Parquet files register into Iceberg without
     rewriting.

   The move is mechanical: write each ragged group with `to_dataframe().to_parquet()`
   and add it as an asset.
4. **Revisit C when either trigger occurs:**
   - (a) we stage many events and want one queryable catalog across them (FlowAlert
     batches, an Alaska landslide inventory), or
   - (b) multiple writers update shared data, e.g. a rolling Synoptic or seismic
     archive like the PRISM cache.
5. **Don't adopt miniSEED/ASDF as the stored form.** Waveforms are arrays in Zarr, with
   NSLC codes, sample rate and start time on each channel, which is everything needed to
   build an `obspy.Stream`. gaia-cli has no obspy dependency, so that conversion stays
   in users' own code ([obspy-dependency.md](obspy-dependency.md)).

## Resolved questions

1. **Users of the dense `stations/*` grid:** only the Skagit notebooks and peters-dome
   `make_report.py`. gaia-cli's own code and examples switch to the ragged layout and
   the `to_grid` helper. Revisit if that causes problems.
2. **Output target:** `-o` is always a folder, never a `.zarr` path. Every run, with
   `stage.yaml` or with flags, writes:

   ```
   <output>/
     <name>.json     # STAC Item (the manifest)
     <name>.zarr/    # the store
     stage.yaml      # resolved configuration
   ```

   `<name>` is the file's `name`, or the folder's basename for a flag-only run.
   `gaia catalog rebuild` finds the Item and checks its `data` asset as it does now.
3. **String columns in ragged arrays:** tested with zarr 3.4.1, xarray 2026.9.0 and
   xvec 0.5.2. Station ids, report types, QC strings (including empty strings), times,
   values, `row_size` and point geometries all round-trip with no custom encoding,
   using the same `swap_dims` + `encode_cf` path as `io.save_xvec_zarr`. The only
   caveat is the fixed-width UTF-32 string type noted under A′.
4. **Arraylake on Kopah/B2:** revisit when a trigger for C occurs.

## Open questions

None.
