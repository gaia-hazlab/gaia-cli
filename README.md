# gaia-cli

[![DOI](https://zenodo.org/badge/1156542089.svg)](https://doi.org/10.5281/zenodo.23771204)


Stage geospatial datasets for an event (an area of interest and a time window) into
one analysis-ready Zarr store: weather stations, PRISM, HRRR, land cover, elevation,
soil and burn severity.

**Documentation: <https://gaia-hazlab.github.io/gaia-cli/>**

```bash
git clone https://github.com/gaia-hazlab/gaia-cli
cd gaia-cli
pixi run gaia stage --config conus-precip \
  -i https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json \
  -s 2025-12-01 -e 2025-12-20 -o skagit
```

See [examples/skagit-precip.ipynb](examples/skagit-precip.ipynb) for an analysis of the result, and
[examples/pioneer-landlab.ipynb](examples/pioneer-landlab.ipynb) for staging Landlab inputs.

## Staging on GitHub Actions

This is a convenience for creating datasets without requiring any local setup, but you must have *write* permissions to this repository to run GitHub actions.

The [Stage Event](.github/workflows/stage-event.yml) workflow collates and stages data required for landlab landslide modeling in Zarr format. It optionally aligns the datasets and and uploads the results to a public object store. Run it from the
[Actions tab](https://github.com/gaia-hazlab/gaia-cli/actions/workflows/stage-event.yml)
or with the GitHub CLI:

```bash
gh workflow run stage-event.yml -R gaia-hazlab/gaia-cli \
  -f config=examples/pioneer.yaml -f name=pioneer \
  -f align="--crs utm --res 10 --aoi examples/pioneer-aoi.geojson" -f landlab_inputs=true
```

See the workflow file for all inputs, and
[stage files](https://gaia-hazlab.github.io/gaia-cli/cli#stage-files) for writing
your own `config`.

## Public datasets

Staged stores shared publicly can be browsed in
[STAC Browser](https://browser.moregeo.it/external/s3.kopah.uw.edu/gaia/scratch/public/catalog.json)
or on a map in
[stac-map](https://developmentseed.org/stac-map/?href=https://s3.kopah.uw.edu/gaia/scratch/public/catalog.json).

## Development

```bash
pixi run test        # tests
pixi run lint        # ruff
pixi run typecheck   # ty
pixi run docs        # preview the user docs (docs/) at http://localhost:3000
```

User documentation lives in [docs/](docs/) and is published to GitHub Pages from
`main`. Developer notes (design, roadmap, per-dataset details) live in [dev/](dev/).
