---
title: CLI
---

| command | does |
|---|---|
| [`gaia stage <dataset>`](#stage) | stage one dataset |
| [`gaia stage --config`](#stage-files) | stage several datasets from a stage file or preset |
| [`gaia align`](#align) | put staged rasters on one grid and time step, in a new store |
| [`gaia catalog rebuild`](#catalog) | build a STAC catalog of the stores under a prefix |
| [`gaia prism fill`](#prism-cache) | fill the shared PRISM cache (maintainers) |

Every option is listed in the [reference](reference.md) and in `gaia <command> --help`.

(stage)=
## `gaia stage <dataset>`

```bash
AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json

gaia stage prism -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit   # timed: needs -s and -e
gaia stage dem   -i $AOI -o skagit --skip-existing               # static: no window
```

- `-i` is a vector file or URL holding a single Polygon.
- `-s` and `-e` are UTC. A date-only `-e` includes the whole day. Timed datasets read
  one extra day beyond each end, so that periods crossing the edges (such as the PRISM
  day) are complete.
- `-o` is an output folder, local or `s3://`. The store is named after the folder.
- `--skip-existing` adds a dataset to an existing store. `--overwrite` replaces it.
- `--dry-run` checks the arguments and the AOI without fetching anything.

(stage-files)=
## `gaia stage --config`

A stage file lists the AOI, the window and the datasets for one event:

```yaml
name: skagit-2025-12                  # names the store
aoi: https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
window: 2025-12-01/2025-12-20         # whole days, end included
datasets:
  synoptic: {}
  prism: {variables: [ppt, tmax]}
  hrrr: {}
  nlcd: {year: 2021}
  dem/usgs10m: {}                     # the key names the product
```

```bash
gaia stage --config skagit.yaml --dry-run        # validate and print, fetch nothing
gaia stage --config skagit.yaml -o skagit
gaia stage --config skagit.yaml -o skagit --only prism,hrrr
```

- Each dataset's options are the same as its CLI flags. A dataset can set its own `aoi`
  and `window`.
- An AOI can also be a bounding box, `{bbox: [west, south, east, north]}`, or a circle,
  `{latitude: 48.5, longitude: -121.5, radius_km: 20}`.
- A dataset that fails doesn't stop the others. Rerun with `--skip-existing` to retry
  only the failed ones.

A **preset** is a built-in stage file without an AOI or window. `conus-precip` stages
Synoptic, PRISM and HRRR precipitation:

```bash
gaia stage --config conus-precip -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit
```

(align)=
## `gaia align`

```bash
gaia align -i skagit/skagit.zarr -o skagit-aligned
gaia align -i skagit/skagit.zarr -o skagit-10m --crs EPSG:32610 --res 10
```

Writes a new store with every raster on one grid and time step. By default it uses
the coarsest of each: for PRISM and HRRR, HRRR's 3 km grid and PRISM's days. Values are
combined in time according to their `cell_methods` (precipitation is summed). In space
each dataset uses its own default resampling (conservative for PRISM and HRRR,
bilinear for DEMs, mode for categories), which `--method` overrides. The input store is
not changed.

(catalog)=
## `gaia catalog rebuild`

```bash
gaia catalog rebuild -p s3://gaia/scratch/public
```

Collects the STAC Items of every store under a prefix into a `catalog.json`. The
public stores are browsable in
[STAC Browser](https://browser.moregeo.it/external/s3.kopah.uw.edu/gaia/scratch/public/catalog.json)
and on a map in
[stac-map](https://developmentseed.org/stac-map/?href=https://s3.kopah.uw.edu/gaia/scratch/public/catalog.json).

(prism-cache)=
## `gaia prism fill`

```bash
gaia prism fill -s 2025-12-01 -e 2025-12-31 -v ppt tmean
```

PRISM limits downloads, so `gaia` reads PRISM from a shared cache. `gaia stage prism`
fills missing days automatically when you have write access to the cache. Otherwise
only cached days can be staged. `gaia prism fill` pre-fills a range.
