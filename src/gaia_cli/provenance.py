"""Provenance attributes for staged Zarr stores, named after STAC fields.

Each product's ``load()`` sets *source* attributes (``license``, ``providers``,
``links``, ...) on the Dataset it returns. ``io.save_datatree`` adds
*run* attributes (``created``, ``gaia:cli-version``, ``gaia:command``) to the
store's root group. See dev/provenance-metadata.md.
"""

from __future__ import annotations

import datetime as dt
import functools
import importlib.metadata
import json
import logging
import shlex
import subprocess
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pyproj
import requests
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import xvec  # noqa: F401  registers the .xvec accessor

log = logging.getLogger(__name__)

# STAC extension schema -> field prefix that marks its use. An extension's
# schema rejects an object that declares it without any of its fields.
STAC_EXTENSIONS = {
    "https://stac-extensions.github.io/processing/v1.2.0/schema.json": "processing:",
    "https://stac-extensions.github.io/version/v1.2.0/schema.json": "version",
    "https://stac-extensions.github.io/scientific/v1.0.0/schema.json": "sci:",
    "https://stac-extensions.github.io/xarray-assets/v1.0.0/schema.json": "xarray:",
}
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
ZARR = "application/vnd+zarr"
YAML = "application/yaml"
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


def extensions(*field_sets: Iterable[str]) -> list[str]:
    """Schema URLs of the STAC extensions used by any of the given fields."""
    fields = [f for fs in field_sets for f in fs]
    return [
        url
        for url, prefix in STAC_EXTENSIONS.items()
        if any(f.startswith(prefix) for f in fields)
    ]


def run_attrs() -> dict:
    """Attributes describing this run of gaia-cli."""
    attrs = {
        "created": _iso(dt.datetime.now(dt.UTC)),
        "gaia:cli-version": cli_version(),
    }
    if command is not None:
        attrs["gaia:command"] = command
    return attrs


def stamp(attrs: Mapping) -> dict:
    """*attrs* plus fresh run attributes, for a store's root group."""
    out = {k: v for k, v in attrs.items() if k != "stac_extensions"}
    out.update(run_attrs())
    if exts := extensions(out):
        out["stac_extensions"] = exts
    return out


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
    contributes the start of its period. Station reports (``time`` along
    ``obs``) give the first and last report.
    """
    if "time" not in obj.coords or obj.time.size == 0:
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


def _bbox(ds: xr.Dataset) -> list[float] | None:
    """WGS84 bounds of *ds*'s raster grid or CF-encoded point geometries."""
    if "geometry_container" in ds:  # xvec.encode_cf output
        decoded = ds.xvec.decode_cf()
        # xvec restores the index of a geometry dimension (stations), not of a
        # geometry coordinate along another one (waveform channels)
        if "geometry" in decoded.xindexes:
            crs = decoded.xindexes["geometry"].crs
        else:
            mapping = ds.geometry_container.attrs.get("grid_mapping", "spatial_ref")
            crs = pyproj.CRS.from_cf(ds[mapping].attrs)
        geoms = gpd.GeoSeries(decoded.geometry.values, crs=crs).to_crs(4326)
        return list(geoms.total_bounds)
    if {"x", "y"} <= set(ds.dims) and ds.rio.crs is not None:
        return list(ds.rio.transform_bounds("EPSG:4326"))
    return None


def bbox(tree: xr.DataTree) -> list[float] | None:
    """WGS84 bounding box covering the data in every group of *tree*.

    Computed from the data rather than the AOI, since clipping keeps every
    touched cell and stations need not fill the AOI.
    """
    boxes = [b for node in tree.subtree if (b := _bbox(node.to_dataset()))]
    if not boxes:
        return None
    west, south, east, north = zip(*boxes, strict=True)
    return [float(min(west)), float(min(south)), float(max(east)), float(max(north))]


def _unique(items: list[dict], key) -> list[dict]:
    seen, out = set(), []
    for item in items:
        if (k := key(item)) not in seen:
            seen.add(k)
            out.append(item)
    return out


def config_asset(href: str) -> dict:
    """The ``config`` asset: the resolved stage file at *href*."""
    return {
        "href": href,
        "type": YAML,
        "title": "Resolved stage file (rerun with gaia stage --config)",
        "roles": ["metadata"],
    }


def stac_item(
    tree: xr.DataTree, store_name: str, assets: Mapping[str, dict] | None = None
) -> dict:
    """A STAC Item for the store *store_name*, written next to it.

    Properties are the root attributes; when the groups hold the provenance
    (``gaia stage``), their time ranges, providers and links are merged,
    and a single shared ``license`` is lifted to the Item. *assets* are
    added after the ``data`` asset.
    """
    groups = [node.attrs for node in tree.subtree]
    props = {
        str(k): v
        for k, v in tree.attrs.items()
        if k not in ("links", "stac_extensions")
    }
    starts = [a["start_datetime"] for a in groups if "start_datetime" in a]
    ends = [a["end_datetime"] for a in groups if "end_datetime" in a]
    if starts and ends:  # RFC 3339 UTC strings sort chronologically
        props.update(datetime=None, start_datetime=min(starts), end_datetime=max(ends))
    else:
        props["datetime"] = props.get("created")
    licenses = {a["license"] for a in groups if "license" in a}
    if "license" not in props and len(licenses) == 1:
        props["license"] = licenses.pop()
    providers = [p for a in groups for p in a.get("providers", [])]
    if providers:
        props["providers"] = _unique(providers, lambda p: p["name"])

    # The description's first words tell stores from `gaia align` apart in a catalog
    kind = "gaia aligned store" if "gaia_align" in props else "gaia store"
    groups_desc = ", ".join(n.path.lstrip("/") for n in tree.leaves if n is not tree)
    if "description" in props:
        props["description"] = f"{kind}: {props['description']}"
    elif groups_desc:
        props["description"] = f"{kind} with groups: {groups_desc}"
    else:
        props["description"] = kind

    box = bbox(tree)
    item = {
        "type": "Feature",
        "stac_version": "1.1.0",
        "id": store_name.removesuffix(".zarr"),
        "geometry": None,
        "properties": props,
        "links": _unique(
            [lk for a in groups for lk in a.get("links", [])],
            lambda lk: (lk["rel"], lk["href"]),
        ),
        "assets": {
            "data": {
                "href": f"./{store_name}",
                "type": ZARR,
                "title": "Zarr store (open with xarray.open_datatree)",
                "roles": ["data"],
                "xarray:open_kwargs": {"engine": "zarr", "consolidated": True},
            },
            **(assets or {}),
        },
    }
    asset_fields = [f for a in item["assets"].values() for f in a]
    item["stac_extensions"] = extensions(props, asset_fields)
    if box:
        west, south, east, north = box
        item["bbox"] = box
        item["geometry"] = {
            "type": "Polygon",
            "coordinates": [
                [
                    [west, south],
                    [east, south],
                    [east, north],
                    [west, north],
                    [west, south],
                ]
            ],
        }
    return item
