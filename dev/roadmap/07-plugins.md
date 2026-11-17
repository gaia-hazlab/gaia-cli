# 07 · Dataset plugins, or a growing list of built-in datasets

Status: undecided, 2026-11-16. Step 7.1 (`gaia check`) is useful either way and comes
first.

Groups keep asking for more datasets, some from obscure sources with complex access
code (Stage IV, [#46](https://github.com/gaia-hazlab/gaia-cli/pull/46): 2,247 lines,
which index monthly tars by byte-range reads of their headers and decompress `.Z`
members). This note weighs keeping every dataset in `src/gaia_cli/datasets/` against
letting other packages provide datasets, and asks what `gaia align` needs from a
dataset that gaia doesn't maintain.

## Where things stand

- **The module contract is already plugin-shaped.** A dataset is a module with `load()`,
  `EXTENT` and `KIND`, plus optional attributes read with `getattr` (`VARIANT`, `SCENES`,
  `TOKENS`, `REGRID_METHOD`, `PAD`, `PROCESS_ENV`, `COMMON_NAMES`;
  [../adding-a-dataset.md](../adding-a-dataset.md)). `gaia stage <name>` and the
  stage-file keys are generated from `load()`'s signature. The runner only looks
  modules up by name in `datasets.DATASETS`.
- **The contract is checked only in a test.** `tests/test_stagefile.py` checks `KIND`
  and `EXTENT` for registered modules, and would not run on an outside module.
- **What `load()` returns is not checked.** `gaia align` relies on conventions that a
  dataset can get wrong without any error. Most of these mistakes give wrong results,
  not a failure:

| `align` assumes | if a dataset gets it wrong |
|---|---|
| rasters have `x`, `y` dims and a `spatial_ref` coordinate | the group is left out of the aligned store, without a message |
| period data has CF time bounds (`time.attrs["bounds"]` → `(time, nv)`) | the data is treated as instantaneous, and time windows come out wrong |
| each timed variable has `cell_methods` | it defaults to `time: mean`, so an accumulation is averaged, not summed |
| categorical variables have `flag_values` | they are interpolated instead of resampled with `mode` |
| SI units, provenance attrs (`title`, `license`, `providers`, `links`), the CF ragged layout for stations | nothing checks |

- **Size.** `src/gaia_cli` is about 10,200 lines, 4,550 of them in `datasets/` (17
  datasets). 36 of the 41 commits touching `datasets/` are by one maintainer.
- **Dependencies used by a single dataset:** `unlzw3` (stage4), `virtualizarr` and
  `obspec-utils` (smap), `arraylake` (era5), `icechunk` (hrrr), `pymseed` (seismic,
  infrasound). Everyone installs all of them.

## Option A: keep growing the built-in list

Pros:

- **Refactoring stays cheap.** A change to `io.with_time_bounds`, the provenance attrs
  or the group layout can update every dataset in one PR, which fits
  [design goal 3](README.md#design-goals) (no backward compatibility yet). This is the
  strongest reason to stay put while the design is still moving.
- **One contract, enforced by one CI:** consistent docs (`docs/datasets.md`),
  coverage checks, offline and online tests.
- **Simple provenance.** `gaia:cli-version` names all the code that wrote a store.
- **Shared helpers get reused.** `opera_rtc`, `nisar_gcov` and `smap` share
  `earthdata.py`; the Planetary Computer datasets share `stac.py`.

Cons:

- **The maintainers own code for sources they don't use.** When GDEX changes its tar
  layout, gaia's online tests break and a gaia maintainer fixes it.
- **Dependencies pile up**, for every user, along with the chance of conflicting pins.
- **Every `gaia` run imports all 17 modules** (`datasets/__init__.py`), so startup slows
  and one broken import breaks every command.
- **Design goal 1** (a codebase you can read in one sitting) no longer holds, and
  datasets are most of the growth.
- **Review is a bottleneck**, and the rule "add only datasets a workflow reads"
  ([04](04-datasets.md)) is about gaia's workflows, not other groups'.
- The 90% coverage gate applies to rarely used code too.

## Option B: datasets from other packages (entry points)

A package declares its modules in its `pyproject.toml`:

```toml
[project.entry-points."gaia.datasets"]
stage4 = "gaia_stage4.dataset"
```

and `datasets/__init__.py` adds `importlib.metadata.entry_points(group="gaia.datasets")`
to `DATASETS`, about 15 lines. Nothing else in the runner changes.

Pros:

- The group that needs a source owns it, tests it and releases it on its own schedule.
- Its dependencies stay in its package; gaia's core stays small.
- Groups can prototype and use a dataset without waiting for a gaia review.

Cons:

- **Internal helpers become a public API.** Plugins will import `io.clip_to_aoi`,
  `io.with_time_bounds`, `provenance.*`, `stations.ragged`, `stac.*`, `window` and
  `fetch`, and changing them breaks outside code, which goes against goal 3. Name a
  small surface on purpose (say, `gaia_cli.api`) and treat everything else as private.
- **Provenance has to name the plugin.** The runner should stamp each group with the
  providing distribution and its version (`gaia:dataset-package`). A stage file that
  names a dataset that isn't installed should fail with "provided by package X", not
  "unknown dataset".
- **A broken plugin must not break `gaia`.** Load each entry point in a `try`, and log
  its failure, so `gaia --help` still works.
- Name clashes with built-in datasets (built-ins win, with a warning), and docs spread
  across repositories.
- **Coupling between datasets** has to go first: `smap` imports `OrbitState` from
  `sentinel1_rtc`. Shared types belong in shared modules.

## What align needs: checks, not a new type

A wrapper type around `load()`'s result adds little: the result is already an
`xr.Dataset` or `xr.DataTree`, and the "schema" is a set of conventions on it. Two
plain functions cover it, without a schema library ([design goal 2](README.md#design-goals)):

- **Check a store or a `load()` result.** Recognize each group's kind from its
  structure (raster, stations, waveforms) and check the conventions in the table
  above. Report errors (`align` would give wrong results) and warnings (missing units,
  missing provenance). Run by `gaia check` on any store, including ones written by
  forks or by hand. The runner can later run it on every `load()` result before
  writing, so built-in and outside datasets meet the same bar.
- **Check a module.** Run `KIND`, `EXTENT`, `VARIANT` and `REGRID_METHOD` checks when a
  module is registered, not only in a test. A `typing.Protocol` for the module lets
  `ty` check it statically too.

Expose both as `gaia_cli.testing.check_dataset(module, aoi, window)` so a plugin's own CI
can run the contract.

## Possible middle path

The core keeps datasets that a GAIA workflow reads and that use providers it already
supports (Planetary Computer, NASA Earthdata, EarthScope FDSN, public COGs and VRTs).
Sources with their own access code (custom archives, unusual formats) become plugins.

`stage4` is the natural trial: it is self-contained (it imports only `fetch`, `io`,
`provenance` and `window`), it brings a dependency only it uses (`unlzw3`), and a
contributor outside the core team wrote it. Moving it out would show exactly which
helpers the public API needs.

## Steps

| # | PR | depends on | notes |
|---|---|---|---|
| 7.1 | `gaia check <store>` | — | Summarize each group (kind, grid, time step, variables) and report what would break or mislead `gaia align`. Exit non-zero on errors. Useful with or without plugins. |
| 7.2 | check every `load()` result in the runner | 7.1 | Fix the built-in datasets it flags first. Remove `align`'s silent `time: mean` default and its silent skip of groups that aren't rasters. |
| 7.3 | `gaia_cli.api`: the helpers plugins may use | decision | Only if plugins go ahead. |
| 7.4 | entry-point discovery, `gaia:dataset-package` provenance, module checks at registration | 7.3 | |
| 7.5 | move `stage4` to its own package | 7.4 | The trial. |
