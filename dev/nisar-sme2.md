# NISAR SME2

Module: `src/gaia_cli/datasets/nisar_sme2.py` · Command: `gaia stage nisar-sme2` ·
Reads like [nisar-gcov.md](nisar-gcov.md#reads), whose helpers it shares.

## Source

[NISAR SME2](https://nisar-docs.asf.alaska.edu/sme2/): JPL's L3 soil moisture,
retrieved from a time series of GCOV backscatter, on the global 200 m EASE-Grid 2.0
(EPSG:6933), from ASF DAAC. 5+5 MHz acquisitions are retrieved at 400 m but still
posted at 200 m. The accuracy goal is 0.06 m³/m³ where the cell is unflagged and its
vegetation water content is below 5 kg/m². There are no retrievals over urban areas,
permanent snow, ice or water, frozen ground, or snow.

- **Maturity.** Only `NISAR_L3_SME2_PROVISIONAL_V1` exists (`C2854344945-ASF`, DOI
  10.5067/NIL3SME2-P1), released 2026-07-20. A search with `dataMaturity=VALIDATED`
  found nothing on 2026-11-23. ASF plans a validated collection in 2027, once validated
  GCOV time series are long enough. `load()` pins the provisional collection and drops
  any other `collectionName`. Add a `maturity` option, as `nisar-gcov` has, when the
  validated collection appears.
- **Reprocessing.** Frames are being reprocessed into the same collection (CRID
  P05036, after P05023), so `load()` keeps the newest `processingDate` of each frame
  acquisition.
- **Known issues** ([provisional](https://nisar-docs.asf.alaska.edu/provisional-known-issues/)).
  - Not fully calibrated.
  - Quad-pol areas have 3 to 6 gaps across the swath, where some algorithms give no
    estimate.
  - There is a permanent instrument gap from 2026-07-27 to 2026-08-10.
- **Granules.** One granule is one frame, one HDF5 file of 95–130 MB, named like GCOV's:
  `NISAR_L3_PR_SME2_<cycle>_<track>_<A|D>_<frame>_…`. Neither the files nor ASF
  provide a DMR++ sidecar, so `dmrpp.py` can't be used.
- **License.** The STAC collection says CC0-1.0, and CMR links the EOSDIS data use
  policy.

## Grids (`/science/LSAR/SME2/grids/`)

All are 2-D, chunked 512 × 512, gzip level 1, `_FillValue` −9999 unless noted
(1485 × 2430 cells in the frame sampled).

| dataset | type | notes |
|---|---|---|
| `soilMoisture` | float32 | m³/m³, 0.02–0.60: the mean of the recommended retrievals of the three algorithms, else of all of them |
| `soilMoistureUncertainty` | float32 | m³/m³ |
| `retrievalQualityFlag` | int16 | bit 0: not recommended; bits 1–3: PMI, TSR, DSG applied (the files' own `description`) |
| `surfaceQualityFlag` | int16 | bits, undocumented in the files (below) |
| `ancillaryData/landCover` | int8, fill −127 | class codes, undocumented |
| `ancillaryData/localIncidenceAngle`, `…Uncertainty` | float32 | degrees |
| `ancillaryData/waterbodyFraction` | float32 | 0–1 |
| `algorithmCandidates/{DSG,PMI,TSR}/…` | | each algorithm's `soilMoisture`, uncertainty, flag and parameters; not read |
| `radarData/frequency{A,B}/…` | | `sigma0HH`, `sigma0HV`, looks, noise; not read |
| `xCoordinates`, `yCoordinates`, `EASEGridColumnIndex`, `EASEGridRowIndex`, `projection` | | cell centers, global grid indices, CF grid mapping |

- **Flags.** The product specification (JPL D-107677 Rev C, 2024-11-08) lists
  `retrievalQualityFlag` bits that disagree with the files' description (recommended,
  attempted, success, range, σ anomalies). The files are product specification 1.0.2,
  so `load()` follows their description, as CF `flag_masks`/`flag_meanings`. The
  specification gives `surfaceQualityFlag` bits 0–7: static water, urban,
  precipitation, snow or ice, permanent snow or ice, frozen ground, terrain, and dense
  vegetation. The files also have bit 8 (values of 256 and above, seen over Odessa), so
  `load()` stores no meanings for it until a newer specification confirms them.
- **Grid.** GDAL's HDF5 driver reads the CRS, geotransform and nodata from the CF
  attributes, as it does for GCOV. A frame's edge is a whole number of 200.179 m cells
  from x = 0 and y = 0 (column 28125 of the global grid is x = −58635 cells), so
  odc-stac's default anchor puts its output on the native cells with no resampling.
  Checked 2026-11-23: the staged cells equal the file's, cell for cell.
- **One header per frame.** A frame's window of the global grid was the same on every
  date: tracks 5 and 77, frames 25 and 26, about 14 dates each from June to October 2026.

## Frames are time steps, not mosaicked

Adjacent frames of a pass overlap. Over Odessa, frames 25 and 26 of track 77 overlap
across about 11,000 cells, and they disagree:

| 2026-10-17, overlap | |
|---|---|
| `soilMoisture` equal | 0.1 % of cells |
| difference | mean 0.004, mean absolute 0.015, RMS 0.027 m³/m³; correlation 0.73 |
| `retrievalQualityFlag` equal | 56–64 % (three dates) |
| `surfaceQualityFlag` equal | 100 % |
| `sigma0HH` | median relative difference 2.7 % |

Each frame is retrieved on its own, so a first-wins mosaic (as `nisar-gcov` does)
would pick one retrieval arbitrarily and leave a seam. Instead each frame is its own
time step at its start time: the frames of a pass start about 34 s apart, and the
`frame` coordinate names each one. `stac.group(items, frame)` makes each frame a scene,
so `pre`/`post` windows, `select` and `count` choose frames.

## Reads

The reads use GCOV's pattern: sign each file's URL with the token
(`earthdata.sign`), point each asset at its `HDF5:"/vsicurl/…"://…` subdataset, and
load with `stac.load_scenes`, in worker processes (`PROCESS_ENV`).

- **2 MiB requests.** The files use HDF5 paged aggregation with 2 MiB pages
  (GCOV's are 4 MiB). Opening and reading the four default grids over Odessa (63 × 63
  cells) from one file, with `CPL_DEBUG` (2026-11-23):

  | `CPL_VSIL_CURL_CHUNK_SIZE` | requests | time |
  |---|---|---|
  | default (16 kB) | 19 | 1.6 s |
  | 2 MiB | 6 | 1.5 s |
  | 4 MiB | 5 | 1.75 s |

  The time is about the same from Seattle. Fewer requests should matter more on a
  slower link, as they did for GCOV.
- **Timing.** `gaia stage nisar-sme2` over Odessa, 2026-10-01 to 10-20: 8 frames, 25 s
  end to end, about half of it the reads in 4 processes.
- **Ruled out for now.**
  - DMR++: no sidecars.
  - Downloading whole files (`nisar-sme-stager`): 100 MB a frame is reasonable, but a
    small AOI needs a few MB of it. Washington's first year is 475 frames, 49 GiB.
  - Icechunk stores and processing in AWS us-west-2 (`nisar-sme-stager`) are possible
    later steps, not part of this baseline.
