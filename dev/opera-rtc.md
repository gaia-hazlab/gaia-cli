# OPERA RTC

Module: `src/gaia_cli/datasets/opera_rtc.py` · Command: `gaia stage opera-rtc` ·
Earthdata helper: `src/gaia_cli/earthdata.py` · Plan: [roadmap/06](roadmap/06-nasa-earthdata.md)

## Source

[OPERA RTC-S1](https://hyp3-docs.asf.alaska.edu/guides/opera_rtc_product_guide/):
JPL's radiometric terrain correction of each Sentinel-1 IW burst. ASF DAAC distributes
it as CMR collection `OPERA_L2_RTC-S1_V1` (`C2777436413-ASF`; STAC
`cmr.earthdata.nasa.gov/stac/ASF/collections/OPERA_L2_RTC-S1_V1_1`).

- **Coverage.** It starts on 2016-04-14. OPERA's own production began in 2022-01; the
  earlier years were back-filled through HyP3. It covers land north of 60° S, and new
  bursts appear within about 12 hours.
- **Granules.** One granule is one burst,
  `OPERA_L2_RTC-S1_<T131-279929-IW2>_<acquired>Z_<processed>Z_S1A_30_v1.0`. Its COGs
  are `_VV`, `_VH` (or `_HH`, `_HV`) and `_mask`, plus an `.h5` holding orbit metadata.
- **License.** CMR says CC0-1.0 and links the NASA Earthdata data use policy.

## Search

`earthdata.search` wraps `asf_search.search(dataset="OPERA-S1", processingLevel="RTC")`.
asf_search is used rather than CMR STAC because CMR STAC items carry only `datetime`.
asf_search's GeoJSON has everything the scene index needs:

- `platform`, `flightDirection`, `pathNumber` (track) and `orbit` (absolute)
- `operaBurstID` and `polarization`
- `processingDate`
- the S3 URL of every file

`orbit_state` and `relative_orbit` become asf_search's `flightDirection` and
`relativeOrbit` filters. `earthdata.to_item` turns each feature into a pystac Item with
the same `sat:` properties as Planetary Computer's items, so `stac.group`,
`stac.select` and `stac.load_scenes` work unchanged. A burst and acquisition may have
more than one processing, and `earthdata.newest` keeps the latest `processingDate`.

## Reads

GDAL reads the COGs remotely, subsetting them lazily while the store is written, after
`load()` returns. Tested on 2026-11-09.

- **Hrefs** are `https://cumulus.asf.earthdatacloud.nasa.gov/OPERA/<s3 key>`, built from
  each file's S3 URL. That host serves every granule, old and new. Older granules list
  `datapool.asf.alaska.edu` URLs instead, which 307 to `cumulus.asf.alaska.edu`.
- **Signed before reading** (`earthdata.sign`). After passes are chosen, one request per
  needed COG goes to cumulus with `Authorization: Bearer <EARTHDATA_TOKEN>`. Cumulus
  answers with a 303 to a CloudFront URL signed for one hour, and that URL replaces
  the asset's href. GDAL then reads anonymously, as with Planetary Computer's
  SAS-signed hrefs (`stac.sign`). The token is used only for those requests; it is in
  no href, attribute or written file, and the online test checks that.
  - **One hour to write.** The store has to be written within the signatures' hour:
    the regional Peters Dome run takes about 20 s.
  - **Errors.** A 302 to the Earthdata login page, or any other non-redirect, raises an
    error naming `EARTHDATA_TOKEN`, since an expired token (60 days) is the usual
    cause.
- **GDAL settings** (the one `odc.stac.configure_rio` call, in `io.py`):
  - `CPL_VSIL_CURL_USE_HEAD=NO`. The signed URLs answer GETs only (HEAD gets 403), and
    a HEAD to cumulus is redirected to a HEAD-only signature (`/head-…`).
  - `GDAL_HTTP_VERSION=1.1`. Over HTTP/2, CloudFront occasionally dropped a stream
    under odc-stac's parallel reads. GDAL logged `response_code=0`, didn't retry, and
    read the empty response as a file "not recognized as being in a supported
    format". That failed about every other regional run, each time on a different
    COG, while 192 parallel requests to cumulus and 384 to one signed URL from
    `requests` all succeeded. With HTTP/1.1, 6 of 6 regional runs succeeded; without
    it, 6 of 6 failed.
- **Ruled out:**
  - A global bearer header (`GDAL_HTTP_HEADERS`, or `GDAL_HTTP_AUTH=BEARER` with
    `configure_rio`) is sent to every host, and Planetary Computer's Azure storage
    answers it with 401.
  - A GDAL config file with a `[credentials]` section scoped to cumulus works, but GDAL
    loads it only when it starts. It can't hold the token if it's committed, since
    GDAL doesn't expand `${EARTHDATA_TOKEN}`, and it can't be set from `load()` once an
    earlier dataset has used GDAL.
  - A path-specific header set at runtime (`VSISetPathSpecificOption` through ctypes)
    also works, without the one-hour limit, but needs ctypes. Its failures turned out
    to be the HTTP/2 problem above.
  - `CPL_VSIL_CURL_USE_S3_REDIRECT=NO`, coincident's setting, doesn't help: the
    failures weren't about cached redirects.
  - Downloading whole bursts would give up remote subsetting.

**Headers for the grid.** The search gives no `proj:` fields, so `earthdata.set_proj`
reads one COG header per burst ID, in parallel, after selection. Only the chosen bursts
are read. A burst ID keeps its grid on every date. For example, T073-155356-IW1 on
2023-09-05 and 09-17 is on EPSG:32604 with origin (507810, 6865830) and 1636 × 3320
cells, and so are its static layers. The headers are read from the signed URLs.

## Passes and values

- **Passes.** Bursts are grouped into passes by `(platform, sat:absolute_orbit)`, as
  `sentinel1-rtc` groups frames, and each pass is one time step at its first burst's
  start time. Where adjacent bursts and subswaths overlap, the first valid value wins.
- **Gamma0.** `vv`, `vh` (or `hh`, `hv`) are linear-power gamma0, float32,
  `units: "1"`. The COGs' nodata is NaN, and zero is set to NaN too. The default
  regrid method is `average`, the mean of linear power, as for `sentinel1-rtc`.
- **Mask.** `mask` is each pass's own layover and shadow mask, which the product guide
  prefers over the static one. It is uint8 CF flags: 0 `valid`, 1 `shadow`,
  2 `layover`, 3 `layover_and_shadow` and 255 `no_data`. 255 is also outside the AOI
  and the bursts, as Sentinel-2's SCL keeps 0 `no_data`. The assets declare
  `raster:bands` (`uint8`, nodata 255), so odc-stac reads the mask as uint8 and
  mosaics past the padding.
- **Coordinates along `time`.** These are `platform` (`Sentinel-1A`), `orbit_state`,
  `relative_orbit`, `absolute_orbit`, `burst_ids` (space-separated) and `window`.
- **Size.** At 30 m a pass is about one ninth of `sentinel1-rtc`'s. The regional Peters
  Dome AOI (4 passes, 34 bursts, 3066 × 3443 cells) is 362 MB and was staged in about
  25 s.

## Peters Dome

These are the same passes as `sentinel1-rtc`:

- The source AOI sees descending track 131 on 2023-09-09 and 09-21.
- The regional AOI adds track 102 on 09-07 and 09-19.
- `pre: {before: 30d, select: last}`, `post: {after: 30d, select: first}` and
  `match: relative_orbit` give all four.
