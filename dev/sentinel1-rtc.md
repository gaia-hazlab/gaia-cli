# Sentinel-1 RTC

Module: `src/gaia_cli/datasets/sentinel1_rtc.py` · Command: `gaia stage sentinel1-rtc` ·
STAC helper: `src/gaia_cli/stac.py`

## Source

Planetary Computer's [`sentinel-1-rtc`](https://planetarycomputer.microsoft.com/dataset/sentinel-1-rtc)
collection: Sentinel-1 IW GRD scenes radiometrically terrain corrected by Catalyst
(PCI), as 10 m COGs of gamma0 in each polarization, in the scene's UTM zone. Global,
2014-10 on, updated within days. License CC-BY-4.0 (attribute Copernicus Sentinel data
and Catalyst; the group's `providers` and `links` come from the STAC collection).

Access is anonymous: `stac.token` fetches a SAS token from
`/api/sas/v1/token/sentinel-1-rtc` (valid about a day) and `stac.sign` appends it to the
blob-storage hrefs, so the `planetary-computer` package isn't needed. A token is
reused until it is within an hour of expiry, since the COGs are read while the store
is written, after `load()` returns.

## Scenes and passes

The search returns one STAC item per GRD frame (25 s, about 170 km along track). The frames
of one pass that the AOI crosses are mosaicked into **one time step**:

- Frames are grouped by `(platform, sat:absolute_orbit)`. Grouping by `solar_day`
  would merge an ascending and a descending pass on the same day, which is common at
  high latitude; grouping by item id would split one pass into several steps.
- `time` is the first frame's `datetime`, so a step's time can be a few seconds before
  the frame that covers the AOI.
- Where frames overlap, odc-stac keeps the first valid value.
- Each step carries `platform`, `orbit_state`, `relative_orbit`, `absolute_orbit` and
  `scene_ids` (the frames' item ids, space-separated) as coordinates along `time`. This
  is the scene index ([roadmap 03](roadmap/03-storage.md)): select a track with
  `ds.where(ds.relative_orbit == 131, drop=True)`.

Options filter the search: `polarizations` (default `[vv, vh]`; scenes without all of
them are skipped, e.g. HH/HV scenes over Greenland), `orbit_state` and
`relative_orbit` (STAC `query` filters). No scenes is an error.

## Choosing passes

A stage file can choose passes instead of keeping every one: `pre`/`post` windows with
`select` and `count`, and `match: relative_orbit` to choose on each track separately
([stage-file.md](stage-file.md#choosing-scenes)). `stac.select` picks among the passes
(frames already grouped) before anything is loaded, so only the chosen passes are read,
and each gets a `window` coordinate. Over the regional Peters Dome AOI,
`pre: {before: 30d, select: last}`, `post: {after: 30d, select: first}` and
`match: relative_orbit` give 2023-09-07 and 09-19 on track 102, and 09-09 and 09-21 on
track 131. Without `match`, they give 09-09 (track 131) and 09-19 (track 102): two
viewing geometries, which a change map shouldn't mix.

## Values

- **Linear power** gamma0, float32, `units: "1"`. Convert to dB with
  `10 * np.log10(ds.vv)`. Keeping linear values means `gaia align`'s default
  `average` (`REGRID_METHOD`) is the mean of power, the correct way to aggregate.
- The assets' nodata (`-32768`) and zero or negative values are NaN.
- Values are instantaneous: there is no `time_bounds`.

## Grid

The native grid of the scenes: the CRS most items are in (`stac.native_crs`) at their
own resolution (10 m), with odc-stac's pixel grid matching the COGs'. Items in another
UTM zone (an AOI near a zone edge) are reprojected onto it, with a warning. Clipping is
`io.clip_to_aoi`. Data are read lazily (dask), so the runner's size limit applies
before anything is read: a 100 × 90 km AOI over one month (4 passes) is about 3 GB.

## Coverage notes

- Sentinel-1B failed in December 2021; until Sentinel-1C (operational 2025) most places
  have a 12-day repeat from Sentinel-1A alone.
- Peters Dome (2023-09-14 landslide): the source area is seen only by descending
  relative orbit 131 (2023-09-09 and 09-21); orbit 102 covers the western part of the
  regional AOI.

## Not yet

- `pre`/`post` windows with `select` and `match: relative_orbit` (step 2.4b).
- `gaia align` drops the scene coordinates along `time`.
