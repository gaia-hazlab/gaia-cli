"""Reads of NASA Earthdata HDF5 files through their DMR++ sidecars.

A DAAC's DMR++ file (``<file URL>.dmrpp``, OPeNDAP) indexes an HDF5 file's
layout: each array's dtype, attributes and filters, and each chunk's byte
offset and size. VirtualiZarr's ``DMRPPParser`` turns it into a Zarr store
whose chunks are HTTP byte ranges of the file: no libhdf5, so reads run in
threads. ``smap`` and ``smap-l4`` read NSIDC's files this way (dev/smap.md#reads).
"""

from __future__ import annotations

import asyncio
import base64
import threading
import warnings
from collections.abc import Mapping, Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from typing import Any, cast

import numpy as np
import zarr

from . import earthdata, fetch

# Chunks are read as byte ranges in aligned blocks, each fetched once, so
# neighboring chunks (an array's sit together in the file) share a request
BLOCK = 2**20
FETCHES = ThreadPoolExecutor(32)


class _Bytes:
    """What obspec's ``get`` returns, for bytes already in memory."""

    def __init__(self, data: bytes):
        self.data = data

    def buffer(self) -> bytes:
        return self.data


class Granule:
    """One file as an obspec store for VirtualiZarr's ``DMRPPParser``.

    The parser reads ``<url>.dmrpp`` with ``head`` and ``get`` or
    ``get_ranges``: it's fetched once, with the token, and kept, so each group
    parsed reuses it. The chunks are read with ``get_range_async``, as byte
    ranges of the URL the DAAC signs for the file, got with the first chunk:
    signed URLs last an hour, so each file is signed when it's read.
    VirtualiZarr's ManifestStore drops a URL's query, and obstore escapes one
    in a path, so a signed URL needs a store of its own. Each aligned
    :data:`BLOCK` is fetched once, and the reads share it.
    """

    def __init__(self, url: str):
        self.url = url
        self.signed: str | None = None
        self.dmrpp: bytes | None = None
        self.blocks: dict[int, Future] = {}
        self.lock = threading.Lock()
        self.sign_lock = threading.Lock()

    def _dmrpp(self) -> bytes:
        with self.lock:
            if self.dmrpp is None:
                resp = fetch.get(
                    f"{self.url}.dmrpp",
                    headers={"Authorization": f"Bearer {earthdata.token()}"},
                )
                if resp.status_code == 404:
                    # Every SMAP file sampled had one (dev/smap.md#reads)
                    name = self.url.rsplit("/", 1)[-1]
                    raise RuntimeError(
                        f"no DMR++ index for {name} ({name}.dmrpp, HTTP 404), "
                        "which gaia reads it through; leave its time step out "
                        "of the window"
                    )
                resp.raise_for_status()
                self.dmrpp = resp.content
            return self.dmrpp

    def _signed(self) -> str:
        with self.sign_lock:
            if self.signed is None:
                self.signed = earthdata.signed_url(self.url)
            return self.signed

    def head(self, path: str) -> dict:
        return {"size": len(self._dmrpp())}

    def get(self, path: str) -> _Bytes:
        return _Bytes(self._dmrpp())

    def get_ranges(
        self,
        path: str,
        *,
        starts: Sequence[int],
        ends: Sequence[int] | None = None,
        lengths: Sequence[int] | None = None,
    ) -> list[bytes]:
        if ends is None:
            assert lengths is not None
            ends = [s + n for s, n in zip(starts, lengths, strict=True)]
        data = self._dmrpp()
        return [data[s:e] for s, e in zip(starts, ends, strict=True)]

    def _get(self, block: int) -> bytes:
        headers = {"Range": f"bytes={block * BLOCK}-{(block + 1) * BLOCK - 1}"}
        resp = fetch.get(self._signed(), headers=headers)
        resp.raise_for_status()
        return resp.content

    async def get_range_async(
        self,
        path: str,
        *,
        start: int,
        end: int | None = None,
        length: int | None = None,
    ) -> bytes:
        if end is None:
            assert length is not None
            end = start + length
        first, last = start // BLOCK, (end - 1) // BLOCK
        with self.lock:
            for b in range(first, last + 1):
                if b not in self.blocks:
                    self.blocks[b] = FETCHES.submit(self._get, b)
            futures = [self.blocks[b] for b in range(first, last + 1)]
        data = b"".join([await asyncio.wrap_future(f) for f in futures])
        return data[start - first * BLOCK : end - first * BLOCK]


def open_arrays(
    url: str, paths: Sequence[str], arrays: Mapping[str, Sequence[str]]
) -> dict[str, zarr.Array]:
    """The arrays at *paths* (``/<group>/<name>``) of the file at *url*, read
    through its DMR++ sidecar. *arrays* names every array in each group, so
    the parser skips the chunk lists of those not asked for. VirtualiZarr's
    parser makes one group per call."""
    from obspec_utils.registry import ObjectStoreRegistry
    from virtualizarr.parsers import DMRPPParser

    store = Granule(url)
    # The file and its .dmrpp, by their folder
    registry = ObjectStoreRegistry({url.rsplit("/", 1)[0]: cast(Any, store)})
    groups: dict[str, list[str]] = {}
    for path in paths:
        group, name = path.strip("/").split("/")
        groups.setdefault(group, []).append(name)
    out = {}
    for group, names in groups.items():
        skip = [a for a in arrays[group] if a not in names]
        parser = DMRPPParser(group=group, skip_variables=skip)
        with warnings.catch_warnings():  # e.g. the string arrays it can't parse
            warnings.simplefilter("ignore")
            manifest = parser(f"{url}.dmrpp", registry)
        z = zarr.open_group(manifest, mode="r", zarr_format=3)
        out.update({f"/{group}/{n}": cast(zarr.Array, z[n]) for n in names})
    return out


def fill_value(attrs: dict, dtype: np.dtype) -> float:
    """An array's ``_FillValue``: VirtualiZarr leaves the floats' as base64 of
    their bytes (its Zarr fill_value is 0)."""
    fill = attrs["_FillValue"]
    if isinstance(fill, str):
        raw = base64.b64decode(fill)
        fill = np.frombuffer(raw, "<f8" if len(raw) == 8 else dtype)[0]
    return float(fill)


def read(
    arrays: Mapping[str, zarr.Array], paths: Sequence[str], rows: slice, cols: slice
) -> list[np.ndarray]:
    """*rows* and *cols* of each of *arrays* at *paths*, their chunks all at
    once."""
    with ThreadPoolExecutor(len(paths)) as pool:
        return list(pool.map(lambda p: np.asarray(arrays[p][rows, cols]), paths))


def read_file(
    url: str,
    paths: Sequence[str],
    arrays: Mapping[str, Sequence[str]],
    rows: slice,
    cols: slice,
) -> list[np.ndarray]:
    """*rows* and *cols* of each array at *paths* in the file at *url*: one
    DMR++ fetch, then the arrays' chunks."""
    return read(open_arrays(url, paths, arrays), paths, rows, cols)
