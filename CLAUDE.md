# Claude Instructions

## Python Environment
Always use `pixi run python ...` when running Python commands. Do NOT create virtual environments (no `python -m venv`, `conda create`, `uv venv`, etc.). All Python execution should go through pixi.

## Common Commands

```bash
# Run tests
pixi run test

pixi run typecheck   # ty check src/

# Lint code & format
pixi run lint        # ruff check src tests
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
  __init__.py   # Public exports
  __main__.py   # python -m gaia_cli entrypoint
  cli.py        # CLI app & subcommand registration
  io.py         # Shared I/O: AOI loading, clipping, Zarr writing
  synoptic.py   # Synoptic API station precipitation
  prism.py      # PRISM daily precipitation (STAC)
  hrrr.py       # HRRR analysis precipitation (Icechunk)
  all.py        # "stage all" – loads every product into an xr.DataTree

tests/
  test_all.py
```


