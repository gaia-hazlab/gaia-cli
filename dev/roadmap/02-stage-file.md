# 02 · Stage file: AOIs, time windows and per-dataset options

## Problem

`gaia stage <dataset> -i AOI -s START -e END` assumes one AOI and one window. Real
workflows break that assumption in three ways:

1. **Datasets need different AOIs and windows.** Peters Dome uses 40 min of seismic data
   from a 400 km circle, 3 months of met data from a 100 km box, and 3 months of imagery
   from a 15 km box.
2. **Some datasets don't take a time window at all.** DEMs, land cover and soils are
   static, while ArcticDEM strips are multitemporal over 18 years.
3. **Satellite observations are irregular.** "The data between `start` and `end`" is
   rarely what is wanted. Instead it is "the last clear scene before the event and the
   first after", or "a median composite of a pre-fire season and a post-fire season"
   (landlab's HLS `windows.pre/post`), or "two passes on the same relative orbit
   bracketing the event" (Peters Dome S1). Model output (HRRR, ERA5) is regular, so
   start/end works there.

Today these needs live in three places: `stage_all.sh` + `config/event.yml` (Peters
Dome), `base.yaml` (landlab), and `-s/-e` flags (gaia-cli).

## Proposal

A small YAML file, `stage.yaml`, describing one event. The name is generic and ties it
to `gaia stage`. Values resolve in this order:

```
dataset entry in the file  >  file defaults  >  dataset's built-in defaults
```

The file and the data-selection flags (`-i`, `-s`, `-e` and dataset options) are
**mutually exclusive**. A run is described either entirely by `stage.yaml` or entirely
by flags, so there is no merge logic. Only run controls (`-o`, `--only`, `--dry-run`,
`--overwrite`, `--skip-existing`) combine with the file.

```yaml
# stage.yaml
name: peters-dome-2023
time: 2023-09-14T02:06:49Z         # event anchor; enables relative windows. Optional.
origin: {latitude: 63.1455, longitude: -151.1904}     # for station distance/back-azimuth. Optional.
attrs:                             # ACDD global attributes, copied to the root group and STAC Item
  title: Peters Dome rock-ice avalanche, 2023-09-14
  references: https://doi.org/10.1029/2025JF008577

aoi: aoi/peters_dome_local.geojson # default AOI (path or URL)
aois:                              # named AOIs for datasets that need another one
  regional: {bbox: [-152.2, 62.75, -150.2, 63.55]}               # [west, south, east, north], as in STAC
  waveforms: {latitude: 63.1455, longitude: -151.1904, radius_km: 400}

window: 2023-06-01/2023-09-21                  # default window (STAC/OGC interval string)
output: s3://gaia/scratch/peters-dome          # output folder (see 03); -o overrides it

datasets:
  seismic:
    aoi: waveforms
    window: {before: 17min, after: 23min}      # relative to `time`
    channel: "HH?,BH?"                         # FDSN codes; per station, the first band it has
  infrasound:
    aoi: waveforms
    window: {before: 17min, after: 53min}
    channel: "BDF,HDF"
  era5: {aoi: regional}                        # inherits the default window
  sentinel1-rtc:
    window:
      pre:  {before: 30d, select: last}        # most recent acquisition before `time`
      post: {after: 30d,  select: first}       # first acquisition after `time`
    match: relative_orbit                      # pre/post must share an orbit
  sentinel2:
    window:
      pre:  {before: 30d, select: last}
      post: {after: 30d,  select: first}
    max_cloud_cover: 60                        # filter on STAC eo:cloud_cover
  dem: {}                                      # static: no window; Cop30 by default
```

The same format covers the landlab postfire case:

```yaml
name: bolt-creek-fire-2022
time: 2022-09-10                   # ignition, or the post-image date
aoi: aoi/bolt_creek.geojson
datasets:
  dem/usgs10m: {}                  # see 04, DEM resolution
  nlcd: {year: 2021}
  soil: {product: solus100}
  burn_severity: {fire: bolt creek}
  prism: {window: 2022-10-01/2023-05-31, variables: [ppt, tmin, tmax]}
```

### CLI behavior

```bash
gaia stage --config stage.yaml                      # every dataset in `datasets:`
gaia stage --config stage.yaml --only seismic,era5  # a subset
gaia stage --config stage.yaml --dry-run            # resolve and validate only
gaia stage --config stage.yaml -s 2023-09-14        # error: flags and file are exclusive
gaia stage prism -i aoi.geojson -s 2025-12-01 -e 2025-12-20 -o out/   # flags only, no file
```

- Every run writes the same output folder ([03](03-storage.md#stac-item-as-the-manifest)).
  A flag-only run builds the same in-memory configuration as a file, goes through the
  same runner, and writes its resolved `stage.yaml` too. That makes every output
  repeatable with `gaia stage --config <output>/stage.yaml`, and there is one code path.
- `gaia stage all` is replaced by `gaia stage --config` with a built-in `conus-precip`
  preset that lists synoptic, prism and hrrr. No alias is kept
  ([design goal 3](README.md#design-goals)). A preset has no AOI or window, so it is
  the one exception to exclusivity: `gaia stage --config conus-precip -i … -s … -e …`.
  The flags fill what the preset leaves out, so there is still nothing to merge.
- Runs continue past a failed dataset, report every failure at the end, and skip groups
  that already exist (`--skip-existing`). This matches `stage_all.sh`.
- `--dry-run` prints the resolved AOI, window and options for every dataset, plus the
  coverage check ([04](04-datasets.md)), without fetching anything. This is cheap, and
  it is what reviewers and new users will want first.
- The resolved file is written next to the store as `stage.yaml` and listed as an asset
  of the event's STAC Item ([03](03-storage.md#stac-item-as-the-manifest)), so the
  output records exactly what produced it.

### Window grammar

Keep it small. Every form resolves to one or more `(name, start, end)` intervals, plus
an optional selection rule:

| form | meaning | needs `time` |
|---|---|---|
| `{start, end}` | absolute interval | no |
| `"start/end"` | the same, as a STAC/OGC `datetime` interval string; `..` leaves one end open | no |
| `{before, after}` | `[time - before, time + after]` | yes |
| `{start}` / `{end}` combined with `before`/`after` | mixed absolute and relative | yes |
| `{pre: {...}, post: {...}}` | named sub-windows, each in one of the forms above; group gets a `window` dim/coord | if relative |
| `select: all \| first \| last \| nearest` (+ `count: N`) | for scene-based datasets: which acquisitions inside the window to keep | no |
| omitted | dataset default (static datasets ignore windows) | no |

Durations use pandas shorthand (`30d`, `17min`, `6h`), parsed with `pd.Timedelta`, not
offset aliases. Allowed units are `s`, `min`, `h`, `d` and `w`. `pd.Timedelta` has no
months or years, so `M` is rejected instead of meaning month-end; use absolute dates for
seasonal windows. `pd.Timedelta` also accepts ISO 8601 durations (`PT17M`), so those work
too without being the documented form.

Selection is applied **after** search and **before** loading pixels, so "last scene
before the event" never reads a whole season of scenes. Only scene-based datasets
(STAC/CMR) accept `select`. Model and station datasets reject it with a validation error.

### Time conventions

The three repos disagree today, so pick one rule and document it:

- Every interval is half-open in UTC: `[start, end)`.
- A **date-only `end` includes that whole day**, which matches landlab's "inclusive UTC
  days" and current `gaia stage -e 2025-12-20` usage. `end: 2025-12-20` means
  `< 2025-12-21T00:00Z`.
- Naive timestamps are UTC. Offsets are honored and converted.
- Products with period values keep CF `time_bounds`, as today.

### Validation

Port the landlab approach:

- Strict keys, with a typo reported as an error and a suggested key.
- Every problem is collected and reported together before any network call.
- Each dataset declares its options once, as a dataclass or as the `stage()` signature,
  so the file and the CLI validate against the same definition.

Because the file and flags are exclusive, no precedence merging is needed, so
`cyclopts.config.Yaml` adds little. The loader reads YAML with `pyyaml` (already
installed through other packages; it becomes a direct dependency), resolves named AOIs
and windows, and builds the same option objects the flags build.

### Names and syntax borrowed from standards

There is no standard for a multi-source data request file. The nearest relatives are
query APIs: STAC API Item Search, OGC API – EDR, openEO `load_collection` and FDSN web
services. The file borrows their names where a concept matches, so users can carry
knowledge across, and uses its own syntax only where nothing fits.

| concept | file syntax | follows |
|---|---|---|
| Point location | `origin: {latitude, longitude}`. Never a bare `[lat, lon]` pair, which is the opposite order from GeoJSON and `bbox`. Not called `location`, which FDSN uses for the location code (`00`, `10`). | FDSN / QuakeML |
| Bounding box | `bbox: [west, south, east, north]` in EPSG:4326. Not the CDS/ECMWF `area: [N, W, S, E]` order. | STAC, GeoJSON |
| Inline geometry | `intersects:` with a GeoJSON geometry | STAC API |
| Circle | `{latitude, longitude, radius_km}`. Same keys as FDSN, but radius in km instead of FDSN's `maxradius` in degrees. | FDSN (adapted) |
| Absolute interval | `"2023-06-01/2023-09-21"`, `..` for an open end; or `{start, end}` | RFC 3339 / STAC `datetime` |
| Timestamps | RFC 3339 (`2023-09-14T02:06:49Z`); naive means UTC | RFC 3339 |
| Durations | pandas shorthand (`17min`, `30d`) | gaia's choice; see [Window grammar](#window-grammar) |
| Seismic selection | `network`, `station`, `location`, `channel`, with `?` and `*` wildcards | FDSN dataselect / station |
| Cloud filter | `max_cloud_cover`, applied to `eo:cloud_cover` | STAC EO extension |
| Root attributes | `attrs:` keys are ACDD global attributes (`title`, `summary`, `references`, `keywords`, `id`) | ACDD / CF |
| Dataset keys | short gaia names (`sentinel1-rtc`, `dem`). Each group's provenance records the upstream STAC collection id (`sentinel-1-rtc`, `cop-dem-glo-30`). | gaia's choice |

Event-relative windows (`before`/`after`) and `pre`/`post`/`select`/`match` are gaia's
own syntax. ISO 8601 has `duration/end` intervals, but nothing for an offset from an
anchor time or for choosing the last clear scene. This is the part the file adds beyond
the query APIs.

### What stays out of the file

- Grids and resolution for conversion: these go in `gaia convert landlab --res 10`.
- Credentials: environment variables only (`SYNOPTIC_TOKEN`, `EARTHDATA_TOKEN`,
  `AWS_PROFILE`).
- Event-specific science, such as published force inversions: this stays with the
  study's own code, which can write extra groups into the same store.

## Resolved questions

1. **Name:** `stage.yaml`, passed with `gaia stage --config`. It is generic (not every
   use is an event) and ties the file to the command.
2. **One event per file:** yes. Consequences:
   - Each file produces one output folder with one store and one STAC Item, and
     `gaia catalog rebuild` collects the Items into a catalog. Batches are a shell loop
     over files (`for f in events/*.yaml; do gaia stage --config $f; done`).
   - Shared settings (provider options, a common window) are repeated in each file.
     Built-in presets cover the common sets of datasets.
   - Queries across events mean opening several stores. At tens of events that is fine;
     at hundreds it is trigger (a) for Iceberg in [03](03-storage.md#recommendation).
   - Nothing blocks accepting several files or a directory later.
3. **Inline AOIs:** yes. `bbox`, `intersects` and the circle form, alongside paths and
   URLs.
4. **Station distance and back-azimuth** use `origin` when given, otherwise the AOI
   centroid, and the attrs record which was used.
5. **Clearest-scene selection:** not in v1. `max_cloud_cover` plus `first`/`last` only.
6. **Event lookup by `comcat_id`:** not in v1.
7. **JSON Schema for the file:** not in v1.
8. **File and flags:** mutually exclusive (from [03](03-storage.md#resolved-questions)).

## Open questions

None.
