# 06 · NASA Earthdata: OPERA RTC-S1 bursts and static layers

Plan for a second Sentinel-1 RTC source, OPERA RTC-S1 from ASF DAAC, and the shared
NASA Earthdata access that later datasets (NISAR GCOV first; HLS and SMAP, see
[05](05-pr-sequence.md#later-on-request)) will reuse. Status: 6.1 (`opera-rtc`) is
built ([../opera-rtc.md](../opera-rtc.md)); 6.2 (`opera-rtc-static`) is next.

## Why

- **One processor, one grid.** OPERA's RTC-S1 is JPL's processing of each Sentinel-1
  IW burst at 30 m, on a fixed grid per burst ID. That makes it a stable, documented
  alternative to Catalyst's 10 m GRD-based RTC on Planetary Computer
  ([../sentinel1-rtc.md](../sentinel1-rtc.md)).
- **The layers each pass needs.** Each burst has its own layover/shadow mask. Static
  layers per burst ID (local incidence angle, ellipsoid incidence angle, number of
  looks, gamma0 to beta0/sigma0 factors) give the viewing geometry that PC's collection
  doesn't publish.
- **Earthdata access built once.** EDL bearer-token auth and authenticated GDAL reads
  are needed by every NASA dataset after this one.

## The source

| | RTC-S1 | RTC-S1-STATIC |
|---|---|---|
| CMR short name / concept | `OPERA_L2_RTC-S1_V1` / `C2777436413-ASF` | `OPERA_L2_RTC-S1-STATIC_V1` / `C2795135174-ASF` |
| asf_search | `dataset=OPERA-S1, processingLevel=RTC` | `processingLevel=RTC-STATIC` |
| Temporal extent | 2016-04-14 on (OPERA production since 2022-01; earlier back-filled through HyP3) | one granule per burst ID, dated 2014-04-03 |
| Spatial extent | land north of 60° S | same |
| Granule | one burst: `OPERA_L2_RTC-S1_T131-279929-IW2_20230921T162934Z_<generated>Z_S1A_30_v1.0` | `OPERA_L2_RTC-S1-STATIC_<burst>_20140403_S1A_30_v1.0` |
| COGs | `VV`, `VH` (or `HH`, `HV`) float32, nodata NaN; `mask` uint8 | `local_incidence_angle`, `incidence_angle` (degrees, float32, NaN), `mask`, `number_of_looks`, `rtc_anf_gamma0_to_beta0`, `rtc_anf_gamma0_to_sigma0` |
| Mask | 0 valid, 1 shadow, 2 layover, 3 both, 255 nodata | same codes over the full burst footprint |
| License | CMR says CC0-1.0, and links the NASA Earthdata data use policy | |

The values are linear gamma0 power, about 7 MB per burst per polarization. Grids are
UTM (polar stereographic at high latitude) at 30 m. A burst ID's grid is identical
across dates and identical to its static layers. One burst checked: T073-155356-IW1 on
2023-09-05 and 09-17, plus its static layers, all on EPSG:32604 with origin
(507810, 6865830) and 1636 × 3320 cells.

Over Peters Dome in September 2023, the regional AOI has 34 bursts: track 102 on 09-07
and 09-19 (4 bursts each) and track 131 on 09-09 and 09-21 (13 each). Those are the same
passes as `sentinel1-rtc` gives, so `pre`/`post` with `match: relative_orbit` carries
over.

## Access, as tested

**Search: `asf_search`** (conda-forge `asf_search-base`, which adds `dateparser`,
`ciso8601`, `pytz` and `tenacity`).

- Each result is a GeoJSON feature carrying what the scene index needs:
  - `platform`, `flightDirection`, `pathNumber`, `orbit` (absolute)
  - `polarization`
  - `processingDate`, to keep the newest generation
  - the URL and S3 URL of each file, with sizes
- It filters server side on `relativeOrbit`, `flightDirection` and `intersectsWith`.
- It already knows the NISAR datasets.
- The alternatives were weaker here:
  - CMR STAC items have only `datetime`, with no orbit or platform.
  - `python-cmr` would leave the UMM-G parsing to us.
  - `earthaccess` brings s3fs/botocore and is mostly about downloads and fsspec.
- `asf_search` only covers ASF's holdings. A later non-ASF dataset (HLS, SMAP) would add
  `python-cmr` for search and reuse the auth below.

**Reads: GDAL, from signed URLs.** One authenticated request per needed COG to
`https://cumulus.asf.earthdatacloud.nasa.gov/OPERA/<s3 key>` returns a CloudFront URL
signed for an hour, which GDAL reads anonymously with ranged requests, over HTTP/1.1.
What was tried and ruled out is in [../opera-rtc.md](../opera-rtc.md#reads).

## Design

### `earthdata.py`: provider helper (new, alongside `stac.py`)

- **`TOKENS = {"EARTHDATA_TOKEN": "an Earthdata Login user token, https://urs.earthdata.nasa.gov (expires after 60 days)"}`.**
  Datasets reuse it. `token()` raises a clear error when the variable is unset.
- **`search(aoi, start, end, **asf_filters)`** wraps `asf_search.search`, and
  **`to_item(feature)`** turns each feature into a `pystac.Item`:
  - properties: `platform`, `sat:orbit_state`, `sat:relative_orbit`,
    `sat:absolute_orbit`, `opera:burst_id` and `processing:datetime`
  - one asset per COG, on the cumulus host
  - `stac.group`, `stac.select` and `stac.load_scenes` then work unchanged.
- **`newest(items)`.** Keeps one item per burst and acquisition: the newest
  `processingDate`.
- **`set_proj(items)`.** Fills `proj:code` and `proj:transform` from one COG header per
  burst ID, because the search has no grid info and odc-stac needs it to stay on the
  native grid.
- **`sign(items, assets)`.** Points the assets at the CloudFront URLs ASF signs for
  them, like `stac.sign` does for Planetary Computer. The GDAL settings it relies on
  (no HEAD, HTTP/1.1) are in the one `odc.stac.configure_rio` call, in `io.py`.

**Later.** NISAR GCOV is HDF5, so it would read signed URLs through an fsspec HTTP
filesystem, with h5netcdf. In-region S3 access
(temporary credentials from `cumulus.asf.alaska.edu/s3credentials`) is optional. Both
fit this module and don't need building now.

### `opera-rtc`: bursts per pass (`datasets/opera_rtc.py`)

This mirrors `sentinel1_rtc.py`, so that the two are interchangeable in a stage file.

- `KIND = "raster"`, `SCENES = True`, `REGRID_METHOD = "average"`, and
  `TOKENS = earthdata.TOKENS`. `EXTENT` is copied from the CMR collection, with a
  comment.
- Options are the same as `sentinel1-rtc`: `polarizations`, `orbit_state`,
  `relative_orbit`, `match`, and `windows` set by the runner.
- **Passes.** Bursts are grouped by `(platform, absolute orbit)`, the same as frames,
  and each pass is one time step mosaicked by `stac.load_scenes`. Adjacent bursts and
  subswaths overlap, and the first valid value wins there, as with frames.
- **Variables.**
  - `vv` and `vh` (or `hh` and `hv`) are linear-power gamma0, float32, `units: "1"`.
  - `mask` is per pass, uint8, stored as CF flags (`flag_values` 0–3 and 255,
    `flag_meanings` `valid shadow layover layover_and_shadow no_data`). `gaia align`
    then uses `mode` for it. It is always included, because
    the product guide says to use each product's own mask rather than the static one.
- **Coordinates along `time`.** These are the same as `sentinel1-rtc` (`platform`,
  `orbit_state`, `relative_orbit`, `absolute_orbit` and `window`), plus `burst_ids` in
  place of `scene_ids`.
- **Size.** At 30 m, it is about one ninth of `sentinel1-rtc` per pass.

### `opera-rtc-static`: viewing geometry per track (`datasets/opera_rtc_static.py`)

This is a separate dataset rather than an option of `opera-rtc`.

- It is static, so it takes no window, like `dem` and `nlcd`. `KIND = "raster"`, and
  `REGRID_METHOD = "bilinear"` suits the smooth angles.
- Options:
  - `layers`: a tuple of a `Literal` of the six layers. The default is
    `local_incidence_angle` and `incidence_angle`.
  - `relative_orbit` and `orbit_state` filters.
- **Output.** Values are mosaicked per track along a `relative_orbit` dimension (with
  an `orbit_state` coordinate), because one place has a different geometry on each
  track. With the same CRS choice as `opera-rtc`, the grids match cell for cell, so
  pairing is `static.sel(relative_orbit=ds.relative_orbit)`.
- **Gap.** ASF has no static layers yet for HH/HV bursts (parts of Greenland and Arctic
  Canada).

## PR sequence

Each step is its own PR, and cleanups go in follow-ups.

| # | PR | content |
|---|---|---|
| 6.1 | **`earthdata.py` + `opera-rtc`** | Add `asf_search-base`. The helper (search → Items, newest generation, signing, header-filled `proj:`); the GDAL settings in `io.py`; the burst dataset. Tests: an offline fixture of asf_search GeoJSON for item building, dedup and grouping, `select` with `match`, and signing; an `online` Peters Dome case skipped without `EARTHDATA_TOKEN`; a `docs/datasets.md` entry and `dev/opera-rtc.md` (the access notes above move there). |
| 6.2 | **`opera-rtc-static`** | Static layers per track, reusing 6.1's helper. Offline grouping test, online case, docs. |

Adding the datasets to the Peters Dome example is out of scope here; it is being done
separately.

## Decided

- Names: `opera-rtc` and `opera-rtc-static`.
- Static layers are a separate dataset, along `relative_orbit`.
- Search uses `asf_search`. Reads use GDAL with ranged requests to URLs ASF signs for
  an hour, over HTTP/1.1.
