# gaia-cli

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

## Development

```bash
pixi run test        # tests
pixi run lint        # ruff
pixi run typecheck   # ty
pixi run docs        # preview the user docs (docs/) at http://localhost:3000
```

User documentation lives in [docs/](docs/) and is published to GitHub Pages from
`main`. Developer notes (design, roadmap, per-dataset details) live in [dev/](dev/).
