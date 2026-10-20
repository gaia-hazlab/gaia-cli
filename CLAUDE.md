# Claude Instructions

## Python Environment
Always use `pixi run python ...` when running Python commands. Do NOT create virtual environments (no `python -m venv`, `conda create`, `uv venv`, etc.). All Python execution should go through pixi.

## Common Commands

```bash
# Run tests
pixi run test
pixi run coverage    # tests + line coverage of src/gaia_cli (pytest-cov), fails below 90% (CI)
pixi run test-online # network smoke tests, marked @pytest.mark.online (skipped by default)

pixi run typecheck   # ty check src/

pixi run docs        # preview the MyST user docs (docs/); docs-build builds them

# Lint code & format
pixi run lint        # ruff check --fix + format, src and tests
pixi run lint-check  # the same, check only (CI)
```

## Code Style

This project uses [ruff](https://docs.astral.sh/ruff/) for linting and formatting.

- Run `pixi run lint` to check for issues; fix violations before committing
- Run `pixi run typecheck` to check type annotations
- Do not add type annotations, docstrings, or comments to code you did not change
- Keep changes minimal and focused — avoid refactoring unrelated code

## Project Structure

- Source lives under `src/gaia_cli/` (src layout)

```
src/gaia_cli/
  __main__.py     # python -m gaia_cli entrypoint
  cli.py          # CLI app & subcommand registration
  runner.py       # "gaia stage": runs a stage file, preset or one dataset's flags; builds the dataset commands
  stagefile.py    # stage.yaml: read, validate, resolve; coverage checks
  aoi.py          # AOI forms: vector file, bbox, intersects, circle
  window.py       # times, durations and window grammar
  io.py           # Shared I/O: clipping, Zarr and STAC Item writing, object stores
  provenance.py   # provenance attributes and STAC Items
  catalog.py      # "gaia catalog rebuild"
  align.py        # "gaia align": time aggregation and regridding onto one grid
  regrid.py       # per-array spatial regridding methods
  presets/        # built-in stage files (conus-precip)
  stations.py     # ragged-array helpers for station data: to_dataframe, to_grid
  prism_cache/    # rolling S3 cache of PRISM COGs (no other gaia_cli imports)
  datasets/       # one module per data source (load(), EXTENT, KIND; see dev/adding-a-dataset.md)
    __init__.py   # DATASETS registry: name -> module
    synoptic.py   # Synoptic API station reports (CF ragged array)
    prism.py      # PRISM daily grids (from prism_cache); also `gaia prism fill`
    hrrr.py       # HRRR analysis (Icechunk)
    nlcd.py, dem.py, soil.py, burn_severity.py   # static layers

tests/          # pytest suite; test_online.py needs network
docs/           # user docs (MyST, GitHub Pages): keep concise; reference.md is generated
dev/            # developer notes: design, roadmap, per-dataset details
```


