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
- **No georeferencing.** The dimensions are unnamed, there's no CF grid mapping, and
  `/Metadata/GridSpatialRepresentation` names the grid only in prose. `latitude` and
  `longitude` are −9999 wherever there's no retrieval, so GDAL finds only scattered
  ground control points. `smap.py` places the grid from EASE-Grid 2.0's constants, as
  NSIDC's own tutorials do (`NSIDC-Data-Tutorials/notebooks/SMAP`):
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
- `long_name`, `valid_min`, `valid_max` and `flag_meanings` come from the first day's
  file (`smap.describe`).
- **Chunks.** A dataset over 10 MB is chunked by whole days, as many as fit in 10 MB
  (`smap.CHUNK_BYTES`; `io` would otherwise keep its dask chunks, one day each). A day
  over 10 MB (a large share of the globe) is split into bands of rows. Smaller datasets
  are one chunk per variable, as everywhere.

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
- **Timings** (2026-11-11), one request each: 0.2 s for a day over Paradise or a
  month over CONUS (31 files), 0.7 s for a year over CONUS (366 files).

## Reads

Each file is read through its **DMR++ sidecar**, `<file URL>.dmrpp`. It's NSIDC's
OPeNDAP index of the file's HDF5 layout: every array's dtype, attributes and filters,
and each chunk's byte offset and size. One day's file has 35.5 MB of it. With the
index, chunks are plain HTTP byte ranges: no libhdf5, so no worker processes and no
`gdal` bindings.

1. `earthdata.sign` gets each file a CloudFront URL signed for an hour (NSIDC answers
   the token with a 303, as ASF does), all before anything is read.
2. `smap.window` turns the AOI's bounds into the grid's rows and columns, plus one on
   each side so rioxarray finds the grid's spacing. `io.clip_to_aoi` drops the extra
   ones.
3. `smap.open_arrays` calls VirtualiZarr's `DMRPPParser` once per group, since it
   makes one group per call. It skips the arrays not asked for (`smap.ARRAYS` lists
   every array in a group) and opens the result with zarr. `smap.Granule` is the
   store it reads through:
   - `head`, `get` and `get_ranges` serve the `.dmrpp`, fetched once with the token;
   - `get_range_async` reads chunks from the signed URL in aligned 1 MiB blocks
     (`smap.BLOCK`), each fetched once, in a pool of 32 threads.
4. The first file is opened in the main process for each array's dtype, fill value
   and attributes (`smap.describe`), and that open is reused for its read. Each other
   file is one dask task (`smap.read_file`): open, then read every array's window at
   once. Zarr requests a window's chunks concurrently: over CONUS, 4,736 one-row
   chunks, which the blocks turn into about 22 requests.

- **Not `earthaccess.virtualize`.** It works, but its `HTTPStore` sends the token with
  every chunk request, and each request then takes NSIDC's 303 to CloudFront: about
  1 s, against 0.15 s for the signed URL. CONUS took 361 s, at zarr's async
  concurrency 10 or 64 alike. It also fetches and parses the `.dmrpp` once per group.
  Installing earthaccess 0.19 needed the `boto3>=1.43.127` pin loosened.
- **Why a store of our own.** VirtualiZarr's `ManifestStore` drops the query from a
  chunk's URL, and obstore escapes `?` in a path, so a signed URL can't be registered
  as is. Sending the signature as CloudFront cookies instead gets a 403, because it
  also covers the URL's `A-userid` and `backend_params`.
- **Public API.** VirtualiZarr PR #902 (open as of 2026-11) moves the parser's core,
  `virtualizarr.parsers.dmrpp.DMRParser`, into pydap. It keeps `DMRPPParser` as the
  wrapper, which is what `smap.py` calls. The wrapper parses the XML once per group,
  so twice per file (0.4 s each, holding the GIL); calling `DMRParser` directly would
  parse it once.
- **Quirks.**
  - VirtualiZarr leaves Zarr's `fill_value` at 0. The real one is the `_FillValue`
    attribute, base64 of a float64 for the floats, which `smap.fill_value` decodes.
  - Units are an attribute.
  - `tb_time_utc` (strings) can't be parsed; the warning is silenced.
- **Every file needs a `.dmrpp`.** Checked for the 8 files from 2025-09-01 to
  2025-09-08; a missing one is an error. TODO: say so clearly, or fall back.
- **Ruled out**, following [nisar-gcov.md](nisar-gcov.md):
  - h5py or h5netcdf over fsspec: an HDF5 read path, subject to libhdf5's lock, that
    reads the file's metadata request by request (below).
  - Downloading whole files.
  - Reading from S3 (`s3://nsidc-cumulus-prod-protected`, with temporary credentials
    from `data.nsidc.earthdatacloud.nasa.gov/s3credentials`): it works only from AWS
    us-west-2.

### Measurements

The 16 default arrays (8 variables, both passes) of
`SMAP_L3_SM_P_E_20250901_R19240_001.h5` (688 MB). One fresh process per run, from
Seattle over HTTPS, 2026-11-11. Values were identical to GDAL's for all 16 arrays.

| | DMR++ (`smap.load(...).compute()`) | GDAL multidimensional API (reads only) |
|---|---|---|
| Paradise, Mt Rainier | 5.8–6.4 s | 10.9–12.2 s |
| CONUS | 6.7–7.0 s | 17.5 s |
| global | not measured | 37.0 s |
| 6 days over CONUS, `runner.stage` (8 files with the padding) | 23.0 s, 23.4 s (threads) | 113–126 s (4 worker processes; 163 s in threads) |

- **One CONUS day** breaks down as: signing 1.5 s, the `.dmrpp` 1.8 s, two XML
  parses 0.8 s, the manifests about 0.6 s, search and provenance 0.6 s, the chunks
  0.8–1.1 s.
- **The `.dmrpp` is each file's floor**, whatever the AOI: about 13 GB of XML for a
  year of files.
- **Without zarr**, fetching the needed byte ranges directly (merged): 3.5–4.1 s for
  Paradise, 3.6–3.9 s for CONUS.
- **TODO:** measure a global read.

### Before: GDAL's multidimensional API

PR #76 (merged 2026-11-12) read the files with GDAL's multidimensional API, from the
`gdal` Python bindings. Each file was one
`gdal.OpenEx("/vsicurl/<signed url>", gdal.OF_MULTIDIM_RASTER)`, then a read of each
array's window (`ReadAsArray(array_start_idx=…, count=…)`). Reads ran in 4 worker
processes (`PROCESS_ENV`, with 16 kB `CPL_VSIL_CURL_CHUNK_SIZE`), since libhdf5 makes
one process's calls one at a time. See `smap.py` and this file at commit `76d8b7e`.

GDAL reads HDF5 two ways, both through libhdf5 (GDAL 3.13.3, `frmts/hdf5`):

- **HDF5 driver subdatasets** (`HDF5:"/vsicurl/<url>"://<group>/<array>`, as
  `nisar-gcov` reads).
  - Every open walks the whole file (`ReadGlobalAttributes`): each group and dataset
    is opened and each attribute read.
  - Finding no CF grid mapping, it then reads 2-D latitude and longitude whole, for a
    lattice of GCPs (`CreateProjections`).
  - Opening one SMAP array took 11.8 s, 146 requests, 14.8 MB. Walking the 204 arrays
    was about 10 s and 108 requests of that; latitude and longitude were 2.8 s,
    35 requests, 12.8 MB.
- **The multidimensional API.**
  - It opens the file and root group (3 requests), then only what's asked for: the
    groups on an array's path, the array's header, its attributes on request, and
    the chunks that intersect the window.
  - That's still about 110 requests per file. They're each array's header and chunk
    index (a B-tree over 1,624 row chunks), scattered through a file without paged
    aggregation, and fetched one after another.

| AOI (window of 9 km cells), 16 kB blocks | multidimensional API | HDF5 driver subdatasets |
|---|---|---|
| Paradise, Mt Rainier (5 × 6) | 12.2 s · 113 requests · 2.3 MB | 22.8 s · 197 requests · 16 MB |
| CONUS (285 × 625) | 17.5 s · 157 requests · 14 MB | 27.9 s · 241 requests · 28 MB |
| global (1624 × 3856) | 37.0 s · 252 requests · 66 MB | 133.8 s · 1,236 requests · 166 MB |

NISAR's files gain little from the multidimensional API: 3.0 s against 3.7–4.2 s for
a frame's 3 arrays over Paradise, with the same 11 requests. With paged aggregation
and 4 MiB requests, their whole metadata is a few requests. NISAR has no DMR++
sidecar.

## Tests

`tests/test_smap.py` builds stand-in files as in-memory zarr arrays, with the
attributes VirtualiZarr gives the real ones (units, base64 `_FillValue`). It patches
`smap.open_arrays` to return them, so the windows, fill values and flags run through
`smap.load`. `smap.Granule` and `smap.fill_value` have tests of their own: the
`.dmrpp` is fetched once with the token, and each block once. `open_arrays` itself,
with the parser and the real files, is covered by `tests/test_online.py`, which
stages a day over the Columbia Basin from NSIDC (needs `EARTHDATA_TOKEN`).
