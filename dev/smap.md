# SMAP

Module: `src/gaia_cli/datasets/smap.py` · Command: `gaia stage smap` ·
Earthdata helper: `src/gaia_cli/earthdata.py` (shared with [opera-rtc.md](opera-rtc.md)
and [nisar-gcov.md](nisar-gcov.md))

## Source

[SMAP Enhanced L3 Radiometer Global and Polar Grid Daily 9 km EASE-Grid Soil
Moisture](https://nsidc.org/data/spl3smp_e/versions/6) (`SPL3SMP_E`), version 6, from
NSIDC DAAC: L-band (1.4 GHz) radiometer brightness temperatures and the soil moisture
retrieved from them, one file a day.

- **Collection.** Pinned by its CMR concept ID, `C2938664763-NSIDC_CPRD` (provider
  `NSIDC_CPRD`, NSIDC's cloud holdings). It's on CMR STAC as
  `cmr.earthdata.nasa.gov/stac/NSIDC_CPRD/collections/SPL3SMP_E_006`; provenance comes
  from there. From 2015-03-31.
- **Enhanced, 9 km.** The radiometer's footprint is about 36 km. `SPL3SMP_E`
  interpolates the same observations onto a 9 km grid (Backus-Gilbert), and its
  brightness temperatures are what NSIDC's SMAP/Sentinel-1 active-passive product
  (`SPL2SMAP_S`) starts from. The native 36 km product (`SPL3SMP`) has the same layout,
  but isn't read.
- **Granules.** One HDF5 file per UTC day, about 700 MB:
  `SMAP_L3_SM_P_E_<YYYYMMDD>_<CRID>_<counter>.h5`, e.g.
  `SMAP_L3_SM_P_E_20250901_R19240_001.h5`. A day can have more than one processing;
  `earthdata.newest` keeps the latest `ProductionDateTime`.
- **Not used.**
  - `SPL1CTB(_E)` (brightness temperature) and `SPL2SMP(_E)` (soil moisture) hold each
    half orbit as a list of cells (`cell_row`, `cell_column`), not a grid.
  - `SPL4SMGP` (Level-4, model-assimilated, with root-zone soil moisture) stores each
    variable as one whole-globe chunk, so any AOI reads all of it, and GDAL finds no
    transform in it.
- **License.** CC0-1.0 per the STAC collection (NASA Earthdata data use policy).

## Grids (`/Soil_Moisture_Retrieval_Data_AM/`, `/…_PM/`)

Each pass is a group of 2-D arrays on EASE-Grid 2.0 Global (EPSG:6933) at 9 km,
1624 × 3856 cells:

| group | pass | names |
|---|---|---|
| `Soil_Moisture_Retrieval_Data_AM` | descending, 6 AM local | `soil_moisture`, … |
| `Soil_Moisture_Retrieval_Data_PM` | ascending, 6 PM local | `soil_moisture_pm`, … (every name ends `_pm`) |

| arrays | type | notes |
|---|---|---|
| `soil_moisture` (`_dca`, `_scah`, `_scav`), `soil_moisture_error` | float32, −9999 | `cm**3/cm**3`. `soil_moisture` is the DCA retrieval |
| `tb_v_corrected`, `tb_h_corrected` (`tb_3_`, `tb_4_corrected`; `tb_v_`, `tb_h_uncorrected`) | float32, −9999 | Kelvin; "corrected" removes water bodies |
| `retrieval_qual_flag` (and per algorithm), `tb_qual_flag_v`, `_h`, `_3`, `_4`, `surface_flag` | uint16, 65534 | bit flags: `flag_masks` (`"1s, 2s, 4s, …"`, HDF5's short-integer literals) and `flag_meanings` |
| `grid_surface_status` | uint16, 65534 | 0 land, 1 water |
| `surface_temperature`, `vegetation_opacity`, `vegetation_water_content`, `albedo`, `roughness_coefficient`, `boresight_incidence`, water and frozen fractions, `clay_fraction`, `bulk_density` | float32, −9999 | the retrievals' inputs |
| `latitude`, `longitude` (and `_centroid`), `EASE_row_index`, `EASE_column_index`, `tb_time_seconds`, `tb_time_utc` | | not read |

- **Layout.** Arrays are chunked by row (1 × 3856), DEFLATE. The files don't use HDF5
  paged aggregation, unlike NISAR's, and their metadata (with a long ISO XML
  attribute) is spread through the file.
- **No georeferencing.** The dimensions are unnamed and GDAL finds only ground
  control points, so `smap.py` places the grid from EASE-Grid 2.0's constants:
  upper-left corner (−17367530.44516138, 7314540.83063540) m and 9008.055210146 m
  cells. Cell centers from these match the files' own `latitude` and `longitude`
  (checked 2026-11-11 on a 2025-09-01 file).
- **Polar groups.** `…_Polar_AM` and `…_Polar_PM` hold the same retrievals on
  EASE-Grid 2.0 North (EPSG:6931), 2000 × 2000 cells. They aren't read: the global
  grid covers 85.044° S to 85.044° N, which is also the dataset's `EXTENT`.

## Time

- **One time step per UTC day**, stamped at 00:00, with `time_bounds` covering the day
  (`io.with_time_bounds(da, 0, 24, "time: point")`): the granule's temporal range in
  CMR. Each cell is a snapshot taken at some time within that day, when a swath
  crossed it. Over a large AOI a pass's cells come from several swaths, hours apart.
  The files' per-cell times aren't read.
- **The passes are variables**, not time steps: `soil_moisture` (AM) and
  `soil_moisture_pm` along one daily axis, as in the files. Two steps a day with the
  same time stamp would make `time` ambiguous for `sel`, `resample` and `gaia align`.
  `--orbit-state` reads one pass's group (`descending` for AM, `ascending` for PM;
  `sentinel1_rtc.OrbitState`, as `opera-rtc` and `nisar-gcov` use). Every file holds
  both, so it halves the arrays read, not the files.
- **`PAD`** is the default day, like PRISM's and Stage IV's.
- **Missing days**, not in the archive, are NaN (flags 65534), with a warning, as in
  Stage IV; none at all is an error. Days whose swaths miss the AOI are NaN without a
  warning: that is SMAP's usual revisit of 2–3 days at mid-latitudes.

## Values

- Float grids are float32 with −9999 as NaN. No quality masking: the flags are staged
  for users to apply.
- Bit flags stay uint16, with 65534 as their Zarr `fill_value` and `flag_masks` as
  numbers (`[1, 2, 4, 8]`). They get `gaia:regrid-method: nearest`, since `gaia align`
  uses `mode` only for `flag_values`; a TODO in `smap.py` marks where `align` should
  learn `flag_masks`. Continuous variables use `REGRID_METHOD = "average"`.
- Units in UDUNITS (`K`, `m3 m-3`), with CF standard names
  `volume_fraction_of_condensed_water_in_soil` (as ERA5's `swvl1`, `swvl2` and NRCS's
  `SMS`) and `brightness_temperature`. No common names: their depths differ (about
  0–5 cm here, 0–7 cm for `swvl1`, NRCS's sensors at 5 cm and below), and a common
  name would mean editing those datasets too.
- `long_name`, `valid_min`, `valid_max` and `flag_meanings` come from the first file
  read.

## Search

`earthdata.cmr_search` queries CMR's granule search
(`cmr.earthdata.nasa.gov/search/granules.umm_json`) through `fetch.get`, with
`concept_id`, `bounding_box` and `temporal`, following the `CMR-Search-After`
header across pages of 2000. `earthdata.cmr_item` turns each UMM-G record into a
pystac Item: its temporal range, first bounding rectangle, `ProductionDateTime` as
`processing:datetime`, and one `data` asset, the HTTPS `GET DATA` URL.

- **Level-3 files are global**, so the bounding box doesn't narrow them.
- **Ruled out.**
  - CMR STAC (pystac-client): its items carry no production time, which
    `earthdata.newest` needs, nor sizes, and name their assets by file path.
  - `python-cmr`: a dependency for one query, and it refuses a spatial filter
    without a short name and version, so the pinned concept ID can't be used alone.
  - `asf_search` (OPERA's and NISAR's): it searches ASF's holdings only.
- **TODO:** search timings.

## Reads

The token gets each file a CloudFront URL signed for an hour (`earthdata.sign`):
NSIDC's `data.nsidc.earthdatacloud.nasa.gov` answers it with a 303, as ASF does. All
URLs are signed before anything is read. GDAL's multidimensional API then reads them
anonymously, from the `gdal` Python bindings (rasterio has no multidimensional API):

1. `smap.window` turns the AOI's bounds into rows and columns of the grid, one more on
   each side so rioxarray finds the grid's spacing, and `io.clip_to_aoi` drops them.
2. `smap.describe` opens the first file in the main process for each array's dtype,
   fill value and attributes.
3. `smap.read` is one dask task per file: one open
   (`gdal.OpenEx("/vsicurl/<signed url>", gdal.OF_MULTIDIM_RASTER)`), then each
   array's window (`ReadAsArray(array_start_idx=…, count=…)`). An open reads the
   file's metadata, the slow part, so every array of a file comes from one task.
4. Reads use `io.py`'s GDAL settings for remote files: no HEAD requests (the signed
   URLs answer GETs only), HTTP/1.1, no directory listing.

- **Read in worker processes** (`PROCESS_ENV`, `io.compute_in_processes`), as
  `nisar-gcov` is: libhdf5 makes one process's calls one at a time.
  - The runner starts workers only for a non-empty `PROCESS_ENV`, so it spells out
    GDAL's default `CPL_VSIL_CURL_CHUNK_SIZE` (16 kB). NISAR's 4 MiB matches its
    files' HDF5 pages; these files aren't paged, so a large request fetches far more
    than an AOI needs.
  - **TODO:** a file's open and an AOI's read (requests, MB, seconds) through the
    multidimensional API against NISAR's reader (an `HDF5:"/vsicurl/…":<path>`
    subdataset through the HDF5 driver), at the default and 4 MiB request sizes, for
    a few AOI sizes up to global.
  - **TODO:** days read with threads against 4 worker processes.
- **Chunks.** A dataset over 10 MB, which `io` would otherwise store with its dask
  chunks of one day, is chunked by whole days, as many as fit in 10 MB
  (`smap.CHUNK_BYTES`). A day over 10 MB (a large share of the globe) is split into
  bands of rows. Smaller datasets are stored as one chunk per variable, as everywhere.
- **Ruled out**, following `nisar-gcov` ([nisar-gcov.md](nisar-gcov.md)):
  - The HDF5 driver's subdatasets (NISAR's reader). An open reads the whole file's
    metadata, every array and attribute, where the multidimensional API opens only the
    groups and arrays asked for. And SMAP gains nothing from it: GDAL finds no
    transform in these files, and there are no frames for `stac.load_scenes` to
    mosaic.
  - h5py or h5netcdf over fsspec: a second HDF5 read path, subject to the same
    libhdf5 lock.
  - Downloading whole files.
  - Reading from S3 (`s3://nsidc-cumulus-prod-protected`, temporary credentials from
    `data.nsidc.earthdatacloud.nasa.gov/s3credentials`): it works only from AWS
    us-west-2.

## Tests

`tests/test_smap.py` builds stand-in files as in-memory GDAL multidimensional
datasets (the MEM driver) with the files' groups and arrays, and patches
`smap.open_file` to return them, so the reads, windows and fill values run through
the same code as the real files. The AOI is near the grid's upper-left corner, so the
stand-ins need only its first rows and columns. `tests/test_online.py` stages a day
over the Columbia Basin from NSIDC (needs `EARTHDATA_TOKEN`).
