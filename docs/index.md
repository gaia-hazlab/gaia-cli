---
title: gaia-cli
---

`gaia` stages the data a hazard analysis needs for one **event**, an area of interest
(AOI) and a time window, into a single cloud-friendly [Zarr](https://zarr.dev) store.

```bash
gaia stage --config conus-precip -i basin.geojson -s 2025-12-01 -e 2025-12-20 -o skagit
```

Each data source has its own API, file format, grid, time convention and units.
`gaia` handles those differences once, so every workflow can start from the same
analysis-ready store and open it with xarray:

```python
import xarray as xr

tree = xr.open_datatree("skagit/skagit.zarr", engine="zarr")
tree["rasters/prism"].precipitation
```

gaia-cli is new and changes often. The documentation describes the latest commit on
`main`.

## Design

- **One store per event, one group per dataset.** Datasets are grouped by kind:
  `rasters/prism`, `rasters/hrrr`, `stations/synoptic`, `dem/cop30`, `soil/solus100`, ….
  The store opens as an [`xarray.DataTree`](https://docs.xarray.dev/en/stable/user-guide/hierarchical-data.html).
- **Native grids, native time steps.** Staging only clips each dataset to the AOI.
  Nothing is reprojected or resampled, so nothing is lost. To compare products,
  [`gaia align`](cli.md#align) writes a *new* store with every raster on one grid
  and one time step.
- **Every observation is kept.** Station reports are stored as they arrive, as a
  [CF ragged array](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#_contiguous_ragged_array_representation),
  not forced onto an hourly grid.
- **Explicit time.** Each time step carries CF `time_bounds` and `cell_methods`, so
  "PRISM day" (12 to 12 UTC) and "HRRR hour-ending" mean the same thing to every
  reader. Times are UTC.
- **Reproducible.** A YAML [stage file](cli.md#stage-files) describes a run. Every
  run saves the resolved file next to its store, and `gaia stage --config` on that
  file repeats the run.
- **Provenance.** Every group records its source's title, license, version and
  links, and the store records the gaia-cli commit and command that made it. Each
  store gets a [STAC](https://stacspec.org) Item, so stores can be catalogued and
  browsed.

## Output

`-o` names an output folder, local or `s3://`:

```
skagit/
  skagit.zarr/   # the store, one group per dataset
  skagit.json    # STAC Item for the store
  stage.yaml     # the resolved stage file: rerun with `gaia stage --config skagit/stage.yaml`
```

## Install

gaia-cli is not yet on PyPI. Run it from a clone with [pixi](https://pixi.sh):

```bash
git clone https://github.com/gaia-hazlab/gaia-cli
cd gaia-cli
pixi run gaia --help
```

To write to an S3-compatible object store, set `AWS_PROFILE` (or `AWS_ACCESS_KEY_ID`
and `AWS_SECRET_ACCESS_KEY`), and the store's endpoint in `~/.aws/config`.
