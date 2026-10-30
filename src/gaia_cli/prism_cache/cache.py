"""A rolling S3 cache of unzipped PRISM daily COGs.

Grids are copied byte-for-byte from the upstream zip packages to keys that
mirror the upstream layout::

    us/an/{res}/{var}/daily/{yyyy}/prism_{var}_us_{code}_{yyyymmdd}.tif
    us/an/{res}/{var}/daily/{yyyy}/prism_{var}_us_{code}_{yyyymmdd}.json

The ``.json`` file is a STAC Item for the grid, holding the PRISM release
number (``prism:release``), when it was fetched (``updated``) and provenance.
It is written after the ``.tif``, so its presence marks a complete entry.
Grids before release 8 are provisional; ``ensure`` re-downloads them when
PRISM publishes a newer release.
"""

from __future__ import annotations

import datetime as dt
import json
import logging
import os
import warnings
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse

import obstore
from obstore.store import S3Store

from . import stac, upstream
from .stac import Entries

log = logging.getLogger(__name__)

DEFAULT_URL = "s3://gaia/public/prism"
DEFAULT_ENDPOINT = "https://s3.kopah.uw.edu"
# Never download the same grid again within this window (PRISM's limit).
MIN_REFETCH_INTERVAL = dt.timedelta(hours=24)


class CacheMiss(Exception):
    """Requested grids are missing and this cache cannot write them."""


def date_range(start: dt.date, end: dt.date) -> list[dt.date]:
    return [start + dt.timedelta(days=i) for i in range((end - start).days + 1)]


def _as_date(value: str | dt.date) -> dt.date:
    if isinstance(value, dt.date):
        return value
    return dt.date.fromisoformat(value)


def _now() -> dt.datetime:
    return dt.datetime.now(dt.UTC)


def _profile_endpoint(profile: str) -> str | None:
    """S3 endpoint URL configured for an AWS *profile* (None if unknown)."""
    import boto3
    import botocore.exceptions

    try:
        client = boto3.Session(profile_name=profile).client("s3")
    except botocore.exceptions.BotoCoreError:
        return None
    return client.meta.endpoint_url.rstrip("/")  # AWS's own endpoint if unset


class PrismCache:
    """Read and fill the PRISM cache.

    Parameters
    ----------
    store : obstore store
        Store rooted at the cache prefix (e.g. ``s3://gaia/public/prism``).
    href_base : str
        Public URL of the cache prefix, used for asset hrefs.
    writable : bool
        Whether this process may download from PRISM and write to *store*.
    """

    def __init__(self, store, href_base: str, writable: bool = False):
        self.store = store
        self.href_base = href_base.rstrip("/")
        self.writable = writable

    @classmethod
    def from_env(cls, writable: bool | None = None) -> PrismCache:
        """Build a cache from the environment.

        ``PRISM_CACHE_URL`` (default ``s3://gaia/public/prism``) locates the
        cache. ``PRISM_CACHE_PROFILE`` names the AWS profile for the cache,
        separate from ``AWS_PROFILE`` used for staged outputs; its
        ``endpoint_url`` and ``region`` are used unless ``PRISM_CACHE_ENDPOINT``
        overrides the endpoint. With neither set the cache is on Kopah
        (``https://s3.kopah.uw.edu``).

        Writing uses boto3 credentials and is enabled, unless *writable* is
        given, when ``PRISM_CACHE_PROFILE`` is set, when ``AWS_PROFILE`` names
        a profile whose ``endpoint_url`` is the cache endpoint, or when
        ``AWS_ACCESS_KEY_ID`` is set and ``AWS_ENDPOINT_URL`` is the cache
        endpoint (so credentials for another store, e.g. Backblaze, leave the
        cache read-only). Otherwise reading is anonymous.
        """
        url = os.environ.get("PRISM_CACHE_URL", DEFAULT_URL)
        profile = os.environ.get("PRISM_CACHE_PROFILE")
        endpoint = os.environ.get("PRISM_CACHE_ENDPOINT")
        if endpoint is None:
            endpoint = (profile and _profile_endpoint(profile)) or DEFAULT_ENDPOINT
        endpoint = endpoint.rstrip("/")
        if profile is None and os.environ.get("AWS_PROFILE"):
            if _profile_endpoint(os.environ["AWS_PROFILE"]) == endpoint:
                profile = os.environ["AWS_PROFILE"]
            else:
                log.info(
                    "AWS_PROFILE=%s does not point at the PRISM cache (%s); "
                    "reading the cache anonymously. Set PRISM_CACHE_PROFILE to "
                    "fill it.",
                    os.environ["AWS_PROFILE"],
                    endpoint,
                )
        if writable is None:
            # In CI the env keys are for the output store, which may not be the cache
            env_keys = os.environ.get("AWS_ACCESS_KEY_ID") and (
                os.environ.get("AWS_ENDPOINT_URL", "").rstrip("/") == endpoint
            )
            writable = bool(profile or env_keys)
        parsed = urlparse(url)
        bucket, prefix = parsed.netloc, parsed.path.strip("/")
        common = {
            "bucket": bucket,
            "prefix": prefix,
            "endpoint": endpoint,
            "virtual_hosted_style_request": False,
        }
        if writable:
            import boto3
            from obstore.auth.boto3 import Boto3CredentialProvider

            # Signed requests need the store's region (e.g. B2's us-west-004)
            session = boto3.Session(profile_name=profile)
            store = S3Store(
                credential_provider=Boto3CredentialProvider(session),
                region=session.region_name or "us-west-2",
                **common,
            )
        else:
            store = S3Store(skip_signature=True, region="us-west-2", **common)
        return cls(store, f"{endpoint}/{bucket}/{prefix}", writable)

    # -- keys -------------------------------------------------------------

    @staticmethod
    def key(var: str, date: dt.date, res: str = "800m") -> str:
        """Cache key without extension."""
        return f"us/an/{res}/{var}/daily/{date:%Y}/{upstream.stem(var, date, res)}"

    def href(self, var: str, date: dt.date, res: str = "800m") -> str:
        """Public URL of a cached COG."""
        return f"{self.href_base}/{self.key(var, date, res)}.tif"

    # -- reading ----------------------------------------------------------

    def _get_json(self, path: str) -> dict | None:
        try:
            return json.loads(bytes(obstore.get(self.store, path).bytes()))
        except FileNotFoundError:
            return None

    def read_collection(self) -> dict | None:
        """Return the cache's root STAC Collection, or None if not written."""
        return self._get_json(stac.COLLECTION_KEY)

    def read_item(self, var: str, date: dt.date, res: str = "800m") -> dict | None:
        """Return the STAC item for one grid, or None if not cached."""
        return self._get_json(f"{self.key(var, date, res)}.json")

    def _existing_items(self, var: str, years: Iterable[int], res: str) -> set[str]:
        keys: set[str] = set()
        for year in years:
            prefix = f"us/an/{res}/{var}/daily/{year}/"
            for batch in obstore.list(self.store, prefix=prefix):
                keys.update(o["path"] for o in batch if o["path"].endswith(".json"))
        return keys

    def read_items(
        self, var: str, dates: list[dt.date], res: str = "800m"
    ) -> dict[dt.date, dict | None]:
        """STAC items for many dates: one listing per year, then parallel
        reads of the items that exist."""
        present = self._existing_items(var, {d.year for d in dates}, res)
        wanted = [d for d in dates if f"{self.key(var, d, res)}.json" in present]
        with ThreadPoolExecutor(max_workers=16) as pool:
            found = dict(
                zip(wanted, pool.map(lambda d: self.read_item(var, d, res), wanted))
            )
        return {d: found.get(d) for d in dates}

    # -- writing ----------------------------------------------------------

    def fill(self, var: str, date: dt.date, res: str = "800m") -> dict | None:
        """Download one grid into the cache and return its STAC item.

        Returns None if the grid is not available upstream.
        """
        # TODO: two writers filling the same grid at once would both download
        # it (counting twice against PRISM's 24h limit). Unlikely with few
        # users; revisit with a lock object written via a conditional put
        # (obstore mode="create"), which kopah's Ceph supports.
        if not self.writable:
            raise CacheMiss("cache is read-only")
        try:
            grid = upstream.download(var, date, res)
        except upstream.NotReleased:
            return None
        key = self.key(var, date, res)
        item = stac.grid_item(grid.tif, grid.meta, self.href_base, key, _now())
        obstore.put(self.store, f"{key}.tif", grid.tif)
        self._put_json(f"{key}.json", item)
        self._register(key, item)
        return item

    def _put_json(self, path: str, doc: dict):
        obstore.put(self.store, path, json.dumps(doc, indent=1).encode())

    def _register(self, key: str, item: dict):
        """Link a newly cached item from its year catalog and the collection."""
        # TODO: read-modify-write without a lock; concurrent writers can drop
        # links. `gaia prism rebuild-catalog` repairs this.
        base = self.href_base
        cat_key = stac.catalog_key(key)
        catalog = self._get_json(cat_key) or stac.new_catalog(base, cat_key, item)
        if stac.add_link(catalog, "item", f"{base}/{key}.json", stac.GEOJSON):
            self._put_json(cat_key, catalog)

        collection = self._get_json(stac.COLLECTION_KEY) or stac.new_collection(
            base, item
        )
        child = stac.add_link(collection, "child", f"{base}/{cat_key}")
        extended = stac.extend_collection(collection, item)
        if child or extended:
            self._put_json(stac.COLLECTION_KEY, collection)

    def rebuild_catalog(self) -> int:
        """Regenerate every catalog and the collection from the cached items.

        Returns the number of items found.
        """
        paths = [
            o["path"]
            for batch in obstore.list(self.store, prefix="us/")
            for o in batch
            if o["path"].endswith(".json") and not o["path"].endswith("catalog.json")
        ]
        with ThreadPoolExecutor(max_workers=16) as pool:
            found = {
                p: it
                for p, it in zip(paths, pool.map(self._get_json, paths))
                if it is not None
            }

        base = self.href_base
        catalogs: dict[str, dict] = {}
        collection = None
        for path, item in sorted(found.items()):
            key = path.removesuffix(".json")
            # Bring items written by older versions up to date
            if stac.relink_item(item, base, key):
                self._put_json(path, item)
            cat_key = stac.catalog_key(key)
            if cat_key not in catalogs:
                catalogs[cat_key] = stac.new_catalog(base, cat_key, item)
            stac.add_link(catalogs[cat_key], "item", f"{base}/{path}", stac.GEOJSON)
            if collection is None:
                collection = stac.new_collection(base, item)
            stac.add_link(collection, "child", f"{base}/{cat_key}")
            stac.extend_collection(collection, item)

        for cat_key, catalog in catalogs.items():
            self._put_json(cat_key, catalog)
        if collection is not None:
            self._put_json(stac.COLLECTION_KEY, collection)
        return len(found)

    # -- the main entry point ---------------------------------------------

    def ensure(
        self,
        variables: Iterable[str],
        start: str | dt.date,
        end: str | dt.date,
        res: str = "800m",
        revalidate: bool = True,
    ) -> Entries:
        """Make sure grids for *variables* over ``[start, end]`` are cached.

        Missing grids are downloaded (requires a writable cache). Cached
        provisional grids (release < 8) are re-downloaded when PRISM has
        published a newer release, unless *revalidate* is False or the grid
        was fetched in the last 24 hours.

        Returns ``{var: {date: STAC item}}`` for every grid that is cached.
        Dates PRISM has not released yet are left out.

        Raises
        ------
        CacheMiss
            If grids are missing and the cache is read-only.
        """
        start, end = _as_date(start), _as_date(end)
        dates = date_range(start, end)
        today = _now().date()
        result: Entries = {}
        for var in variables:
            if var not in upstream.VARIABLES:
                raise ValueError(f"Unknown PRISM variable {var!r}")
            cached = self.read_items(var, dates, res)
            # Dates in the future can't exist upstream; don't ask for them
            missing = [d for d, it in cached.items() if it is None and d < today]
            provisional = [
                d
                for d, it in cached.items()
                if it is not None and stac.release(it) < upstream.FINAL_RELEASE
            ]
            log.info(
                "PRISM %s: %d of %d grid(s) in cache %s (%d provisional, %d missing)",
                var,
                len(dates) - len(missing),
                len(dates),
                self.href_base,
                len(provisional),
                len(missing),
            )

            stale = []
            if revalidate and provisional:
                log.info(
                    "PRISM %s: checking %d provisional grid(s) for newer releases (%s)",
                    var,
                    len(provisional),
                    urlparse(upstream.SERVICE_URL).netloc,
                )
                latest = upstream.releases(var, min(provisional), max(provisional), res)
                for d in provisional:
                    item = cached[d]
                    assert item is not None
                    rel = latest.get(d)
                    if (
                        rel is not None
                        and rel.release > stac.release(item)
                        and _now() - stac.updated(item) >= MIN_REFETCH_INTERVAL
                    ):
                        stale.append(d)

            if missing and not self.writable:
                raise CacheMiss(
                    f"{len(missing)} PRISM {var} grid(s) between {missing[0]} and "
                    f"{missing[-1]} are not cached. Set PRISM_CACHE_PROFILE (or "
                    "AWS_PROFILE) to a profile with write access to the cache to "
                    "fetch them (e.g. `gaia prism fill`)."
                )
            if stale and not self.writable:
                warnings.warn(
                    f"{len(stale)} cached PRISM {var} grid(s) have newer releases "
                    "upstream; using the cached versions (cache is read-only).",
                    stacklevel=2,
                )
                stale = []

            todo = sorted(missing + stale)
            if todo:
                log.info(
                    "PRISM %s: fetching %d grid(s) from %s into the cache "
                    "(%d missing, %d re-released; downloads are serialized, "
                    "%.0f s apart, per PRISM's limits)",
                    var,
                    len(todo),
                    urlparse(upstream.DATA_URL).netloc,
                    len(missing),
                    len(stale),
                    upstream.REQUEST_PAUSE,
                )
            for i, d in enumerate(todo, 1):
                log.info("PRISM %s: [%d/%d] fetching %s", var, i, len(todo), d)
                item = self.fill(var, d, res)
                if item is not None:
                    cached[d] = item

            result[var] = {d: it for d, it in cached.items() if it is not None}
        return result
