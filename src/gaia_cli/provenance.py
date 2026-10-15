"""Provenance attributes for staged Zarr stores, named after STAC fields.

Each product's ``load()`` sets *source* attributes (``license``, ``providers``,
``links``, ...) on the Dataset it returns. The ``io.save_*`` functions add
*run* attributes (``created``, ``gaia:cli-version``, ``gaia:command``) to the
store's root group. See docs/provenance-metadata.md.
"""

from __future__ import annotations

import datetime as dt
import functools
import importlib.metadata
import json
import logging
import shlex
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd
import requests
import xarray as xr

log = logging.getLogger(__name__)

STAC_EXTENSIONS = [
    "https://stac-extensions.github.io/processing/v1.2.0/schema.json",
    "https://stac-extensions.github.io/version/v1.2.0/schema.json",
    "https://stac-extensions.github.io/scientific/v1.0.0/schema.json",
]
# Fields and link relations copied from an upstream STAC collection
COLLECTION_FIELDS = (
    "title",
    "description",
    "license",
    "providers",
    "version",
    "sci:doi",
    "sci:citation",
)
COLLECTION_LINK_RELS = {"license", "about", "cite-as"}
RUN_KEYS = ("created", "gaia:cli-version", "gaia:command", "stac_extensions")
MULTIPLE_LICENSES_NOTE = (
    "Products in this store have different licenses; check each group's "
    "`license` and `links` before sharing."
)

# Set by the CLI entrypoint; None when gaia_cli is used as a library
command: str | None = None


def set_command(argv: Sequence[str]):
    """Record the CLI invocation (``gaia`` plus *argv*) as ``gaia:command``."""
    global command
    command = shlex.join(["gaia", *argv])


@functools.cache
def cli_version() -> str:
    """Git commit of the running gaia-cli, else its package version."""
    pkg = Path(__file__).resolve().parent
    try:
        out = subprocess.run(
            ["git", "-C", str(pkg), "rev-parse", "--show-toplevel", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        )
        toplevel, sha = out.stdout.split()
        # Only trust git when this is gaia-cli's own checkout, not an
        # unrelated repo that happens to contain the environment
        if Path(toplevel).resolve() / "src" / "gaia_cli" == pkg:
            return sha
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    try:
        dist = importlib.metadata.distribution("gaia-cli")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"
    # Recorded by `pip install git+https://...`
    direct_url = json.loads(dist.read_text("direct_url.json") or "{}")
    return direct_url.get("vcs_info", {}).get("commit_id") or dist.version


def _iso(t) -> str:
    return pd.Timestamp(t).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_attrs() -> dict:
    """Attributes describing this run of gaia-cli."""
    attrs = {
        "created": _iso(dt.datetime.now(dt.UTC)),
        "gaia:cli-version": cli_version(),
        "stac_extensions": STAC_EXTENSIONS,
    }
    if command is not None:
        attrs["gaia:command"] = command
    return attrs


def stamp(attrs: Mapping) -> dict:
    """*attrs* plus fresh run attributes, for a store's root group."""
    return {**attrs, **run_attrs()}


def stamp_tree(tree: xr.DataTree) -> dict:
    """Root attributes for *tree*: run attributes, plus a note when the
    groups carry different licenses (the root then has no ``license``)."""
    licenses = {
        node.attrs["license"]
        for node in tree.subtree
        if node is not tree and "license" in node.attrs
    }
    attrs = dict(tree.attrs)
    if len(licenses) > 1:
        attrs["gaia:license-note"] = MULTIPLE_LICENSES_NOTE
    return stamp(attrs)


def derived_from(attrs: Mapping, href: str) -> dict:
    """Root attributes for a store derived from the store at *href*.

    The input's run attributes move to ``gaia:inputs`` so that a chain such
    as ``stage`` → ``align`` keeps every run's commit and time.
    """
    out = {k: v for k, v in attrs.items() if k not in (*RUN_KEYS, "gaia:inputs")}
    run = {k: attrs[k] for k in RUN_KEYS if k in attrs and k != "stac_extensions"}
    out["gaia:inputs"] = [*attrs.get("gaia:inputs", []), {"href": href, **run}]
    out["links"] = [*attrs.get("links", []), link("derived_from", href)]
    return out


def link(rel: str, href: str, **fields) -> dict:
    """A STAC link object."""
    return {"rel": rel, "href": href, **fields}


def time_extent(obj: xr.Dataset | xr.DataArray) -> dict:
    """STAC ``start_datetime``/``end_datetime`` covering *obj*'s time axis.

    Uses the CF ``time_bounds`` when present, so a period-ending value
    contributes the start of its period.
    """
    if "time" not in obj.dims or obj.sizes["time"] == 0:
        return {}
    bounds = obj.time.attrs.get("bounds")
    if isinstance(obj, xr.Dataset) and bounds in obj:
        start, end = obj[bounds].values.min(), obj[bounds].values.max()
    else:
        start, end = obj.time.values.min(), obj.time.values.max()
    return {"start_datetime": _iso(start), "end_datetime": _iso(end)}


def collection_attrs(collection: Mapping | None, href: str) -> dict:
    """Source attributes copied from a STAC *collection* found at *href*."""
    collection = collection or {}
    attrs = {k: collection[k] for k in COLLECTION_FIELDS if k in collection}
    links = [
        lk
        for lk in collection.get("links", [])
        if lk.get("rel") in COLLECTION_LINK_RELS
    ]
    attrs["links"] = [*links, link("derived_from", href, type="application/json")]
    return attrs


def from_stac_collection(href: str) -> dict:
    """Source attributes from the upstream STAC collection at *href*.

    A failed fetch is logged and leaves just a ``derived_from`` link, so it
    never stops a staging run.
    """
    try:
        resp = requests.get(href, timeout=30)
        resp.raise_for_status()
        collection = resp.json()
    except (requests.RequestException, ValueError) as e:
        log.warning("could not read STAC collection %s (%s)", href, e)
        collection = None
    return collection_attrs(collection, href)
