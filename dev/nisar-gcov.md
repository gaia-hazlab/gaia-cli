# NISAR GCOV

Module: `src/gaia_cli/datasets/nisar_gcov.py` · Command: `gaia stage nisar-gcov` ·
Earthdata helper: `src/gaia_cli/earthdata.py` (shared with [opera-rtc.md](opera-rtc.md))

## Source

[NISAR GCOV](https://nisar-docs.asf.alaska.edu/gcov/): JPL's L2 geocoded polarimetric
covariance of each L-band frame, radiometrically terrain corrected (area-based, on the
Copernicus DEM), from ASF DAAC. The [NISAR Cookbook](https://asfopensarlab.github.io/NISAR_Cookbook/)
has worked examples.

- **Maturity.** Each maturity is its own CMR collection, and `--maturity` picks one:
  - `validated` (default): `NISAR_L2_GCOV_V1` (`C3622241360-ASF`, DOI 10.5067/NIL2GCOV-V1),
    released 2026-11-06, acquisitions from 2026-10-21.
  - `provisional`: `NISAR_L2_GCOV_PROVISIONAL_V1` (`C2854338529-ASF`, DOI
    10.5067/NIL2GCOV-P1), released 2026-07-20. It holds acquisitions from 2026-06-17
    to about 2026-10-23, and also reprocessed frames back to 2025-11 on some tracks.
    ASF will reprocess these in 2027.

  ASF advises against mixing maturities in one analysis, so `load()` never does. It
  searches with `dataMaturity` and also drops any feature whose `collectionName`
  isn't the one asked for. Both collections are on CMR STAC
  (`cmr.earthdata.nasa.gov/stac/ASF/collections/<name>_1`); provenance comes from
  there.
- **Known issues** ([provisional](https://nisar-docs.asf.alaska.edu/provisional-known-issues/)).
  - The polarimetric phase isn't fully calibrated. That affects only the off-diagonal
    terms, which aren't read.
  - Frequency B has larger radiometric residuals, and it isn't read either.
  - There is a permanent instrument gap from 2026-07-27 to 2026-08-10.
- **Granules.** One granule is one frame, one HDF5 file of about 6.7 GB holding
  frequencies A (the main band) and B (the 5 MHz side band). The name is
  `NISAR_L2_PR_GCOV_<cycle>_<track>_<A|D>_<frame>_<modes>_<pols>_<S>_<start>_<end>_<CRID>_…_<counter>`,
  e.g. `…_029_172_D_065_4005_DHDH_A_20260906T031847_…_P05023_N_F_J_001`. Here `4005`
  is 40+5 MHz, `DHDH` is dual-pol HH+HV, and `P05023` is the provisional CRID.
- **License.** CMR links the NASA Earthdata data use policy, and the STAC collection
  says CC0-1.0.

## Grids (`/science/LSAR/GCOV/grids/frequencyA/`)

| dataset | type | notes |
|---|---|---|
| `HHHH`, `HVHV` (`VVVV`, `VHVH` quad-pol) | float32, `_FillValue` NaN | gamma0, linear power, `units: 1` |
| `mask` | uint8, `_FillValue` 255 | 0: a sample in the averaging set partly focused or invalid; 1…n: subswath number; 255: outside the swath |
| `numberOfLooks` | float32 | |
| `rtcGammaToSigmaFactor` | float32 | gamma0 × factor = sigma0 (about 0.7 on flat ground at ~45°) |
| `xCoordinates`, `yCoordinates`, `projection` | | cell centers; `projection` is a CF grid mapping with `epsg_code` |

- **Grid.** UTM, or polar stereographic beyond 60°. Spacing is 10 m in 40 MHz mode and
  20 m in 20 or 77 MHz. Chunks are 512 × 512, gzip level 1 with shuffle.
- **Metadata** is in `/science/LSAR/identification/` (`orbitPassDirection`,
  `trackNumber`, `frameNumber`, `zeroDopplerStartTime`, `productVersion`). The
  asf_search properties carry all of it (`pathNumber`, `frameNumber`,
  `flightDirection`, `orbit`, `rangeBandwidth`, `mainBandPolarization`,
  `collectionName`, `processingDate`), so `to_item` doesn't open the file.

## Reads

The reads follow OPERA's pattern ([opera-rtc.md](opera-rtc.md#reads)). The token gets
each file a CloudFront URL signed for an hour (`earthdata.sign`, one request per file
even though every grid is an asset). GDAL then reads it anonymously. Each asset's
href becomes the subdataset `HDF5:"/vsicurl/<signed url>"://science/LSAR/GCOV/grids/frequencyA/<name>`,
so `stac.load_scenes` mosaics and clips frames as it does COGs. Tested 2026-11-09 with
GDAL 3.13.3.

- **`libgdal-hdf5`** is a dependency. conda-forge's `libgdal-core` has no HDF5 or
  netCDF driver.
- **The netCDF driver is no good here.** It would read the CF grid mapping, but on
  macOS it can't open `/vsi` files (it needs Linux userfaultfd). The HDF5 driver does
  pick up the CRS, geotransform and nodata from the CF attributes.
- **Read in worker processes** (`PROCESS_ENV`, `io.compute_in_processes`). An open
  reads the file's metadata, and libhdf5 makes one process's calls one at a time, so
  with threads the opens of every pass queue up. The runner instead reads the dataset
  in 4 spawned processes, then writes it from its own. 9 passes over Paradise
  (a 4.5 × 7.8 km AOI): 102 s with threads, 37 s in processes.
- **4 MiB requests.** The files use HDF5 paged aggregation with 4 MiB pages, and the
  metadata is at the front (the first data chunk starts at about 4.8 MB). GDAL's
  default 16 kB `/vsicurl/` requests take 3–5 s to open a file; with
  `CPL_VSIL_CURL_CHUNK_SIZE` at the page size, about 1 s. The workers set it in their
  environment because GDAL fixes the chunk size at a process's first `/vsicurl/`
  read: a `rasterio.Env` around the load (which odc-stac does carry to its reads)
  can't change it, and setting it in the main process would also change every COG
  read. In processes: 33–35 s with the default, 26–28 s with 4 MiB.
- **Four workers.** Eight were barely faster (24–26 s), and one trial with eight
  workers and 1 MiB requests had a file fail to open.
- **Each pass is about 3 opens** (2 terms and the mask). A worker's later grids of a
  file reuse its `/vsicurl/` cache: about 0.4 s.
- **One header per grid.** `earthdata.set_proj` reads one header per (track, frame,
  bandwidth) rather than per granule. Every date of a frame had the same grid:
  checked for tracks 5 and 172 over Paradise, June to October 2026.
- **Mixed grids.** Frames of one pass can be in different UTM zones (track 5's frame
  25 is in zone 11, track 172's frame 65 in 10). The minority zone is reprojected
  (`stac.load`), at the finest spacing among the chosen frames.
- **Ruled out.**
  - h5py or h5netcdf over fsspec or obstore works, and is what the Cookbook and
    `nisar-proposal-figs` use. It would need a separate mosaic and clip path from
    odc-stac's, and is subject to the same libhdf5 lock.
  - h5py's `ros3` driver with a page buffer (`nisarhdf`'s partial read) reads only
    from S3, in us-west-2.
  - Downloading whole 7 GB files (`nisarhdf` over HTTPS, `nisar-sme-stager`).
  - Virtual references (VirtualiZarr's `HDFParser`): loading from them was blocked
    by zarr-developers/VirtualiZarr#978.
