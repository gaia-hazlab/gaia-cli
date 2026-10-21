"""``gaia stage`` – stage datasets into one store.

Every run goes through :func:`execute`: a stage file, a built-in preset filled
in with ``-i``/``-s``/``-e``, or a single dataset's flags (``gaia stage prism``),
which build the same configuration (see :mod:`gaia_cli.stagefile`).
"""

from __future__ import annotations

import dataclasses
import inspect
import logging
import re
import sys
import typing
from pathlib import Path
from typing import Annotated, Any

import geopandas as gpd
import pandas as pd
import xarray as xr
import yaml
from cyclopts import Parameter
from xarray.backends.file_manager import FILE_CACHE

from . import io, provenance, regrid
from .aoi import aoi_key
from .stagefile import (
    CONFIG_FILE,
    DATASETS,
    PRESETS,
    StageConfig,
    StageFileError,
    group,
    is_url,
    load_aois,
    load_config,
    options,
    outside_extent,
    pad,
    parse_config,
    plain,
    presets,
    read_config,
    resolved,
    takes_origin,
    timed,
    write_stage_file,
)
from .window import unknown

log = logging.getLogger(__name__)


def run(
    cfg: StageConfig,
    aois: dict[str, gpd.GeoDataFrame],
    store: str,
    overwrite: bool = False,
    append: bool = False,
    assets: dict | None = None,
) -> dict[str, str]:
    """Load every requested dataset and write each as a group of one store.

    Each group is written, and the store consolidated and its STAC Item
    rewritten, as soon as it loads. A dataset that fails doesn't stop the
    others; returns ``{name: error}`` for the ones that failed. With
    *append*, the groups are added to the existing store.
    """
    failures = {}
    for req in cfg.requests:
        module = DATASETS[req.dataset]
        kwargs = dict(req.options)
        if req.windows:
            w = req.windows[0].padded(pad(module))
            assert w.start is not None and w.end is not None  # checked in load_config
            # load() takes inclusive bounds; the window's end is exclusive. One
            # microsecond before it keeps a waveform's last second, and gives
            # daily and hourly datasets the same steps as one second would
            kwargs["start_date"] = w.start.isoformat()
            kwargs["end_date"] = (w.end - pd.Timedelta(microseconds=1)).isoformat()
        if takes_origin(module):
            kwargs["origin"] = cfg.origin
        log.info(
            "%s: loading (%s)",
            req.name,
            f"{req.windows[0]}, read as {req.windows[0].padded(pad(module))}"
            if req.windows
            else "static",
        )
        try:
            obj = module.load(aois[aoi_key(req.aoi)], **kwargs)
            nodes = _nodes(module, group(req), obj)
            if method := getattr(module, "REGRID_METHOD", None):
                for node in nodes.values():
                    for var in node.data_vars.values():
                        if {"x", "y"} <= set(var.dims):
                            var.attrs.setdefault(regrid.METHOD_ATTR, method)
            tree = xr.DataTree.from_dict(nodes)
            tree.attrs = dict(cfg.attrs)
            io.save_datatree(
                tree, store, overwrite=overwrite, append=append, assets=assets
            )
            append = True  # later groups go into the same store
        except Exception as e:
            log.error("%s: failed: %s", req.name, e)
            log.debug("%s: traceback", req.name, exc_info=True)
            failures[req.name] = f"{type(e).__name__}: {e}"
    # Close cached remote GDAL handles before interpreter shutdown, where their
    # cleanup can crash sys.excepthook on Python 3.14
    FILE_CACHE.clear()
    return failures


def _nodes(module, path: str, obj: xr.Dataset | xr.DataTree) -> dict[str, xr.Dataset]:
    """The store groups for what *module*'s ``load()`` returned: a Dataset at
    the dataset's group *path*, or each child of a DataTree under it (e.g.
    ``waveforms/seismic/50hz``), with the tree's own attrs on *path*.

    Station and waveform data is CF-encoded for the store with the module's
    ``to_dataset``.
    """
    encode = getattr(module, "to_dataset", lambda ds: ds)
    if not isinstance(obj, xr.DataTree):
        return {path: encode(obj)}
    nodes = {
        f"{path}/{name}": encode(child.to_dataset())
        for name, child in obj.children.items()
    }
    if obj.attrs:
        nodes[path] = xr.Dataset(attrs=obj.attrs)
    return nodes


def execute(
    cfg: StageConfig,
    path: str,
    output: str | None,
    only: str | None = None,
    dry_run: bool = False,
    overwrite: bool = False,
    skip_existing: bool = False,
):
    """Check *cfg* (from *path*) and stage it: the one code path for every run.

    Exits with an error listing every problem, or every failed dataset.
    """
    try:
        problems = []
        if only is not None:
            names = [n.strip() for n in only.split(",") if n.strip()]
            # `--only dem` selects every dem/<variant>
            listed = [r.name for r in cfg.requests] + [r.dataset for r in cfg.requests]
            problems += [
                f"--only: {p}"
                for p in unknown(
                    names, list(dict.fromkeys(listed)), "dataset in the file"
                )
            ]
            cfg.requests = [
                r for r in cfg.requests if r.name in names or r.dataset in names
            ]
        if output is None and not dry_run:
            problems.append("output: required, here or with -o")
        elif output is not None and output.rstrip("/").endswith(".zarr"):
            problems.append(
                f"output: -o is a folder for <name>.zarr, {CONFIG_FILE} and the "
                f"STAC Item; drop the .zarr from {output}"
            )
        if overwrite and skip_existing:
            problems.append("--overwrite and --skip-existing can't be combined")
        if problems:
            raise StageFileError(path, problems)
        aois = load_aois(cfg, path)
    except StageFileError as e:
        sys.exit(f"Error: {e}")

    skipped = {
        req.name: reason
        for req in cfg.requests
        if (reason := outside_extent(req, aois[aoi_key(req.aoi)]))
    }
    if dry_run:
        print(
            yaml.safe_dump(
                resolved(cfg, output, aois, skipped),
                sort_keys=False,
                default_flow_style=None,
            ),
            end="",
        )
        return
    assert output is not None
    for name, reason in skipped.items():
        log.warning("%s: skipped, %s", name, reason)
    requests = [r for r in cfg.requests if r.name not in skipped]
    if not requests:
        reasons = "; ".join(f"{n}: {r}" for n, r in skipped.items())
        sys.exit(f"Error: nothing to stage ({reasons})")

    store = io.store_path(output, cfg.name)
    append = skip_existing and io.output_exists(store)
    if skip_existing:
        io.check_s3_env(output)
    else:
        io.check_output(output, overwrite)
    write_stage_file(
        resolved(dataclasses.replace(cfg, requests=requests), output, aois),
        output,
        merge=append,
    )
    if append:
        existing = io.groups(store)
        for req in requests:
            if group(req) in existing:
                log.info("%s: skipped, %s already exists", req.name, group(req))
        requests = [r for r in requests if group(r) not in existing]
        if not requests:
            log.info("Done: every group already exists")
            return
    failures = run(
        dataclasses.replace(cfg, requests=requests),
        aois,
        store,
        overwrite=overwrite,
        append=append,
        assets={"config": provenance.config_asset(f"./{CONFIG_FILE}")},
    )
    if failures:
        lines = "\n".join(f"  - {n}: {e}" for n, e in failures.items())
        written = len(requests) - len(failures)
        hint = " Rerun with --skip-existing to retry them." if written else ""
        sys.exit(
            f"Error: {len(failures)} of {len(requests)} dataset(s) failed; "
            f"{written} written to {store}.{hint}\n{lines}"
        )
    log.info("Done")


def stage(
    # Collects any other tokens, so data flags get a clear error, not "Unknown option"
    *flags: Annotated[str, Parameter(allow_leading_hyphen=True, show=False)],
    config: Annotated[str, Parameter(name=["--config", "-c"])],
    aoi_path: Annotated[str | None, Parameter(name=["--input", "-i"])] = None,
    start_date: Annotated[str | None, Parameter(name=["--start", "-s"])] = None,
    end_date: Annotated[str | None, Parameter(name=["--end", "-e"])] = None,
    output_path: Annotated[str | None, Parameter(name=["--output", "-o"])] = None,
    only: Annotated[str | None, Parameter(name=["--only"])] = None,
    dry_run: Annotated[bool, Parameter(name=["--dry-run"], negative="")] = False,
    overwrite: Annotated[bool, Parameter(name=["--overwrite"], negative="")] = False,
    skip_existing: Annotated[
        bool, Parameter(name=["--skip-existing"], negative="")
    ] = False,
):
    """Stage every dataset in a stage file, or a built-in preset, into one store.

    The file gives the AOIs, time windows and options for each dataset (see
    dev/stage-file.md). It replaces ``-i``, ``-s``, ``-e`` and dataset
    options, which can't be combined with it. A preset (``conus-precip``) is
    a built-in file without an AOI or window: give those with ``-i``, ``-s``
    and ``-e``.

    Datasets whose extent the AOI or window misses are skipped with a
    warning. A dataset that fails doesn't stop the others; every failure is
    reported at the end.

    Parameters
    ----------
    config : str
        Path or http(s) URL of a ``stage.yaml`` file, or the name of a
        preset. A file read from a URL gives its AOIs as URLs or inline.
    aoi_path : str
        With a preset only: path to a vector file defining the AOI.
    start_date : str
        With a preset only: start (inclusive), a date or UTC time.
    end_date : str
        With a preset only: end; a date includes the whole day.
    output_path : str
        Output folder for ``<name>.zarr``, its STAC Item ``<name>.json`` and
        the resolved ``stage.yaml``. Default: the file's ``output``.
    only : str
        Comma-separated datasets from the file to stage, e.g. ``prism,hrrr``.
    dry_run : bool
        Validate the file, read the AOIs and print the resolved configuration,
        including datasets skipped for coverage, without fetching any data.
    overwrite : bool
        Write into an output folder that isn't empty, replacing its store.
        Default False.
    skip_existing : bool
        Add to an existing store, skipping datasets whose group it already
        has. Use it to retry the datasets that failed. Default False.

    Examples
    --------
    Check, then stage, the Skagit Basin December 2025 atmospheric river:
        $ gaia stage --config examples/skagit-2025.yaml --dry-run
        $ gaia stage --config examples/skagit-2025.yaml -o skagit

    Station, PRISM and HRRR precipitation for any CONUS AOI:
        $ AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
        $ gaia stage --config conus-precip -i $AOI -s 2025-12-01 -e 2025-12-20 -o skagit
    """
    preset = PRESETS / f"{config}.yaml"
    is_preset = config in presets() and not Path(config).exists()
    data_flags = [
        f
        for f, v in (("-i", aoi_path), ("-s", start_date), ("-e", end_date))
        if v is not None
    ]
    if flags or (data_flags and not is_preset):
        sys.exit(
            f"Error: --config can't be combined with {' '.join([*data_flags, *flags])}; "
            f"set the AOI, window and options in {config} instead"
        )
    try:
        if is_preset:
            if aoi_path is None:
                raise StageFileError(
                    config, ["the preset has no AOI; give one with -i (and -s, -e)"]
                )
            doc = read_config(str(preset))
            doc["aoi"] = aoi_path
            # Like a flag-only run, the store is named after the output folder
            if output_path is not None:
                doc["name"] = io.folder_name(output_path)
            window = {"start": start_date, "end": end_date}
            if any(window.values()):
                doc["window"] = {k: v for k, v in window.items() if v is not None}
            cfg = parse_config(doc, config, Path.cwd())
        elif is_url(config):
            cfg = load_config(config)
        elif not Path(config).exists():
            raise StageFileError(
                config,
                [f"file not found, and not a preset ({', '.join(presets())})"],
            )
        else:
            cfg = load_config(config)
    except StageFileError as e:
        sys.exit(f"Error: {e}")
    execute(
        cfg,
        config,
        output_path or cfg.output,
        only=only,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
    )


# The flags of every `gaia stage <dataset>` command, around its load() options
INPUT = Annotated[
    str,
    Parameter(
        name=["--input", "-i"],
        help="Vector file or URL holding the AOI, a single Polygon.",
    ),
]

START = Annotated[
    str,
    Parameter(name=["--start", "-s"], help="Start (inclusive), a date or UTC time."),
]

END = Annotated[
    str,
    Parameter(name=["--end", "-e"], help="End; a date includes the whole day."),
]

OUTPUT = Annotated[
    str,
    Parameter(
        name=["--output", "-o"],
        help="Output folder for the Zarr store (named after the folder), its "
        "STAC Item and the resolved stage.yaml.",
    ),
]

DRY_RUN = Annotated[
    bool,
    Parameter(
        name=["--dry-run"],
        negative="",
        help="Validate the flags, read the AOI and print the resolved "
        "configuration without fetching any data.",
    ),
]

OVERWRITE = Annotated[
    bool,
    Parameter(
        name=["--overwrite"],
        negative="",
        help="Write into an output folder that isn't empty, replacing its store.",
    ),
]

SKIP_EXISTING = Annotated[
    bool,
    Parameter(
        name=["--skip-existing"],
        negative="",
        help="Add to an existing store, unless it already has this dataset's group.",
    ),
]

SHORT_FLAGS = ("-i", "-s", "-e", "-o")


def dataset_command(name: str) -> typing.Callable:
    """The ``gaia stage <name>`` command: the common flags around the options
    of the dataset's ``load()``, each also ``-<first letter>``.

    Its help is the module's docstring and ``load()``'s Parameters section.
    """
    module = DATASETS[name]
    P = inspect.Parameter

    def param(n, hint, default=P.empty):
        return P(n, P.POSITIONAL_OR_KEYWORD, default=default, annotation=hint)

    params = [param("aoi_path", INPUT)]
    if timed(module):
        params += [param("start_date", START), param("end_date", END)]
    params.append(param("output_path", OUTPUT))
    for opt, (hint, default) in options(module.load).items():
        flags = [f"--{opt.replace('_', '-')}"]
        if f"-{opt[0]}" not in SHORT_FLAGS:
            flags.append(f"-{opt[0]}")
        kinds = [typing.get_origin(h) for h in (hint, *typing.get_args(hint))]
        many = list in kinds or tuple in kinds
        params.append(
            param(
                opt, Annotated[hint, Parameter(flags, consume_multiple=many)], default
            )
        )
    params += [
        param("dry_run", DRY_RUN, False),
        param("overwrite", OVERWRITE, False),
        param("skip_existing", SKIP_EXISTING, False),
    ]
    sig = inspect.Signature(params)

    def command(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)
        bound.apply_defaults()
        stage_dataset(name, **bound.arguments)

    command.__signature__ = sig  # ty: ignore[unresolved-attribute]
    doc = inspect.cleandoc(module.load.__doc__ or "")
    section = re.search(
        r"^Parameters\n-+\n.*?(?=^\w+\n-+\n|\Z)", doc, re.MULTILINE | re.DOTALL
    )
    command.__doc__ = "\n\n".join(
        [inspect.cleandoc(module.__doc__ or ""), *([section[0]] if section else [])]
    )
    return command


def stage_dataset(
    name: str,
    aoi_path: str,
    output_path: str,
    start_date: str | None = None,
    end_date: str | None = None,
    dry_run: bool = False,
    overwrite: bool = False,
    skip_existing: bool = False,
    **options,
):
    """Stage one dataset from its CLI flags (``gaia stage prism -i … -s …``).

    Builds the stage file the flags describe and runs it like any other. The
    store is named after the output folder.
    """
    key = name
    if variant_opt := getattr(DATASETS[name], "VARIANT", None):
        # `gaia stage dem -p usgs10m` is the key dem/usgs10m
        key = f"{name}/{plain(options.pop(variant_opt))}"
    doc: dict[str, Any] = {
        "name": io.folder_name(output_path),
        "aoi": aoi_path,
        "datasets": {key: {k: plain(v) for k, v in options.items()}},
    }
    if timed(DATASETS[name]):
        doc["window"] = {"start": start_date, "end": end_date}
    try:
        cfg = parse_config(doc, f"gaia stage {name}", Path.cwd())
    except StageFileError as e:
        sys.exit(f"Error: {e}")
    execute(
        cfg,
        f"gaia stage {name}",
        output_path,
        dry_run=dry_run,
        overwrite=overwrite,
        skip_existing=skip_existing,
    )
