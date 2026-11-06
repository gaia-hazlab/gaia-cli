# Sentinel-2

Module: `src/gaia_cli/datasets/sentinel2.py` · Command: `gaia stage sentinel2` ·
STAC helper: `src/gaia_cli/stac.py`

## Source

Planetary Computer's [`sentinel-2-l2a`](https://planetarycomputer.microsoft.com/dataset/sentinel-2-l2a)
collection: Sentinel-2 MSI Level-2A surface reflectance and scene classification from
ESA's Sen2Cor, as COGs on the 110 km MGRS tiles, in each tile's UTM zone. Global,
2015-06 on. License `proprietary`: the Copernicus Sentinel data terms (free use with
attribution), linked from the group's `links`. Access is anonymous, signed with a SAS
token like `sentinel1-rtc` ([sentinel1-rtc.md](sentinel1-rtc.md#source)).

## Tiles and passes

One STAC item per tile of a datatake. The tiles of one pass (one `s2:datatake_id`)
that the AOI crosses are mosaicked into **one time step**, at the first tile's time
(odc-stac keeps the first valid value where tiles overlap).

- **Reprocessed products.** Planetary Computer lists ESA's reprocessings next to the
  originals: over Peters Dome in September 2023 most tiles appear twice, as baselines
  05.09 and 05.10. Only the newest processing of each tile of each datatake is kept
  (`newest`), and the baseline suffix of `s2:datatake_id` (`_N05.10`) is ignored when
  grouping, so the two processings don't become two passes.
- Each step carries `platform`, `relative_orbit`, `cloud_cover` (the mean
  `eo:cloud_cover` of its tiles, in %), `mgrs_tiles` and `scene_ids` along `time`, and
  `window` with `pre`/`post` windows ([stage-file.md](stage-file.md#choosing-scenes)).

Options: `bands` (common names, default `[blue, green, red, nir, swir16, scl]`) and
`max_cloud_cover` (default 100, no filter), a STAC `query` on `eo:cloud_cover`. That
is the cloud cover of the whole tile, not the AOI's: a mostly cloudy tile can be clear
over a small AOI, and the reverse. No scenes is an error.

| band | asset | pixel size | | band | asset | pixel size |
|---|---|---|---|---|---|---|
| `coastal` | B01 | 60 m | | `nir` | B08 | 10 m |
| `blue` | B02 | 10 m | | `nir08` | B8A | 20 m |
| `green` | B03 | 10 m | | `nir09` | B09 | 60 m |
| `red` | B04 | 10 m | | `swir16` | B11 | 20 m |
| `rededge1`–`3` | B05–B07 | 20 m | | `swir22` | B12 | 20 m |
| `scl` | SCL | 20 m | | | | |

## Values

- **Reflectance**, float32, `units: "1"`. The digital numbers are reflectance × 10000,
  plus 1000 from processing baseline 04.00 (January 2022) on (ESA's `BOA_ADD_OFFSET`),
  so reflectance = (DN − 1000) / 10000 for those and DN / 10000 before, chosen per
  pass from `s2:processing_baseline`. The items carry no `raster:bands` scale or offset.
  DN 0 is NaN.
- **Stored packed**: as uint16 with `scale_factor` 0.0001, `add_offset` −0.1 and
  `_FillValue` 0 (half the size of float32), decoded to float on read. xarray decodes
  them as float64 (Zarr keeps the scale as a JSON number); `.astype("float32")` halves
  that. Values above 6.45 (saturated pixels of pre-2022 baselines) clip.
- **`scl`**, the scene classification, as uint8 CF `flag_values` and `flag_meanings`
  (so `gaia align` uses `mode`): 0 no data (also outside the AOI and tiles), 1
  saturated or defective, 2 topographic or cast shadow, 3 cloud shadow, 4 vegetation,
  5 not vegetated, 6 water, 7 unclassified, 8–9 cloud (medium, high probability), 10
  thin cirrus, 11 snow or ice. Mask clouds with `ds.where(~ds.scl.isin([3, 8, 9, 10]))`.
- Values are instantaneous: there is no `time_bounds`. `REGRID_METHOD = "average"`.

## Grid

The CRS most tiles are in, at the finest requested band's pixel size (10 m by default);
coarser bands are upsampled by nearest neighbour, so each 20 m cell becomes 2 × 2 cells.
With the default bands, a 5 km circle (about 1,000 × 1,000 cells) is about 20 MB per
pass in memory, and 5 MB on disk once packed and compressed.

## Coverage notes

- Revisit is 5 days with two satellites (2A and 2B, then 2B and 2C from 2025), more
  often at high latitude where orbits overlap: Peters Dome had passes on relative
  orbits 86 and 129 every 2–5 days in September 2023.
- Peters Dome: the clear 2023-09-24 pass shows the source area mostly snow-covered
  (SCL 11), matching the paper's note of snowfall right after the event.

## Not yet

- Choosing the clearest scene over the AOI (from `scl`) rather than by tile cloud cover.
