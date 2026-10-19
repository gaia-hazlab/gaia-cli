# Stage file (`stage.yaml`)

A stage file describes one event: its AOIs, time windows and the datasets to stage.
`gaia stage --config` stages every dataset in it into one Zarr store, one group per
dataset (`stations/synoptic`, `rasters/prism`, `rasters/hrrr`, `rasters/nlcd`).
The design and its reasoning are in [roadmap/02-stage-file.md](roadmap/02-stage-file.md).

```bash
gaia stage --config examples/skagit-2025.yaml --dry-run          # validate and print, fetch nothing
gaia stage --config examples/skagit-2025.yaml -o /tmp/skagit.zarr
```

The file replaces the data-selection flags (`-i`, `-s`, `-e` and dataset options),
which can't be combined with `--config`. Only `-o`, `--dry-run` and `--overwrite`
combine with it.

## Example

```yaml
name: peters-dome-2023               # required
time: 2023-09-14T02:06:49Z           # event anchor for relative windows; optional
origin: {latitude: 63.1455, longitude: -151.1904}   # optional
attrs:                               # copied to the store's root group
  title: Peters Dome rock-ice avalanche, 2023-09-14
aois:                                # named AOIs
  regional: {bbox: [-152.2, 62.75, -150.2, 63.55]}
  local: {latitude: 63.1455, longitude: -151.1904, radius_km: 15}
aoi: regional                        # default AOI: a name, path, URL or inline AOI
window: 2023-06-01/2023-09-21        # default window
output: s3://gaia/scratch/peters-dome.zarr   # -o overrides it
datasets:
  hrrr: {aoi: local, window: {before: 6h, after: 1d}, variable: temperature_2m}
  prism: {variables: [ppt, tmax]}
  nlcd: {year: 2021}                 # static: takes no window
```

Each dataset uses its own `aoi` and `window`, otherwise the file's defaults. Its other
keys are the options of that dataset's `load()` function, with the same defaults as
the CLI flags:

| dataset | options |
|---|---|
| `synoptic` | none |
| `prism` | `variables` (list, default `[ppt]`) |
| `hrrr` | `variable` (default `precipitation_surface`) |
| `nlcd` | `year` (default `2021`); static, no window |

## AOIs

| form | example |
|---|---|
| path or URL | `aoi/basin.geojson` (relative to the stage file), `https://…/basin.json` |
| bounding box, `[west, south, east, north]` in EPSG:4326 | `{bbox: [-152.2, 62.75, -150.2, 63.55]}` |
| GeoJSON geometry | `{intersects: {type: Polygon, coordinates: [...]}}` |
| circle, radius in km | `{latitude: 63.1455, longitude: -151.1904, radius_km: 400}` |

## Windows

Every window is a half-open UTC interval `[start, end)`. A date-only `end` includes that
whole day: `2025-12-01/2025-12-20` ends at `2025-12-21T00:00Z`. Naive timestamps are
UTC, and offsets (`2025-12-01T00:00-08:00`) are converted.

| form | meaning |
|---|---|
| `"start/end"`, `..` for an open end | absolute interval |
| `{start, end}` | the same |
| `{before, after}` | `[time - before, time + after]`; a missing side is `time` |
| `{start, after}`, `{before, end}` | mixed absolute and relative |
| `{pre: {...}, post: {...}}` | named sub-windows (for scene-based datasets) |
| `select: all \| first \| last \| nearest`, `count: N` | which scenes to keep (scene-based datasets) |

Durations are a number and a unit (`s`, `min`, `h`, `d`, `w`), such as `17min` or `30d`,
or ISO 8601 (`PT17M`). Months and years aren't durations; use absolute dates.

The current datasets each take one closed window, so `pre`/`post`, `select` and open
ends are rejected for them. Each dataset's own time convention (for example PRISM's
12 UTC day) applies inside the window; see the dataset's doc.

## Validation

Unknown keys (with a suggested key), bad values, missing AOI files and windows that
need a `time` are all reported together, before any data is fetched. `--dry-run` also
reads each AOI once and prints the resolved configuration, including each AOI's bounds.
