"""SMAP L-band brightness temperature and soil moisture (9 km), clipped to an AOI.

Each day in the window is a time step of NSIDC's enhanced Level-3 radiometer
product (SPL3SMP_E version 006): the day's morning (descending) and evening
(ascending) passes on the 9 km EASE-Grid 2.0, the evening's variables with the
files' ``_pm`` suffix. Days missing from the archive are NaN, with a warning.
Read remotely from each day's HDF5 file. Needs EARTHDATA_TOKEN.

Examples
--------
    $ gaia stage smap -i palouse.geojson -s 2025-09-01 -e 2025-09-10 -o palouse-smap
    $ gaia stage smap -i palouse.geojson -s 2025-09-01 -e 2025-09-10 -o palouse-smap --orbit-state descending --variables soil_moisture
"""

from __future__ import annotations

import asyncio
import base64
import logging
import math
import threading
import warnings
import xml.etree.ElementTree as ET
from concurrent.futures import Future, ThreadPoolExecutor
from typing import Any, Literal, cast

import dask
import dask.array
import geopandas as gpd
import numpy as np
import pandas as pd
import pystac
import xarray as xr
import zarr

from .. import earthdata, fetch, io, provenance, regrid
from .sentinel1_rtc import OrbitState

log = logging.getLogger(__name__)

# SPL3SMP_E version 006 in NSIDC's cloud, pinned
CONCEPT_ID = "C2938664763-NSIDC_CPRD"
STAC = "https://cmr.earthdata.nasa.gov/stac/NSIDC_CPRD/collections/SPL3SMP_E_006"
# CMR's, but only to the global grid's edges (85.044°): the polar groups aren't read
EXTENT = {
    "spatial": {"bbox": [[-180.0, -85.044, 180.0, 85.044]]},
    "temporal": {"interval": [["2015-03-31T00:00:00Z", None]]},
}
KIND = "raster"
# `gaia align` default for the measurements; the integer flags use nearest (load)
REGRID_METHOD = "average"
TOKENS = earthdata.TOKENS
# Each file's chunk index comes from its DMR++ sidecar (<file>.dmrpp, OPeNDAP),
# and its chunks are read as byte ranges in threads: no libhdf5, so no worker
# processes. Ranges are fetched in aligned blocks, each once, so neighboring
# rows (each array's chunks sit together in the file) share a request
BLOCK = 2**20
FETCHES = ThreadPoolExecutor(32)

# EASE-Grid 2.0 Global at 9 km: its upper-left corner and cell size (m), and
# its rows and columns
CRS = "EPSG:6933"
X0, Y0 = -17367530.44516138, 7314540.83063540
CELL = 9008.055210146
SHAPE = (1624, 3856)

# Each pass's group in the files, and the suffix of its variables' names
GROUPS = {
    OrbitState.descending: ("Soil_Moisture_Retrieval_Data_AM", ""),
    OrbitState.ascending: ("Soil_Moisture_Retrieval_Data_PM", "_pm"),
}

# The grids of each pass's group, by their AM names (the cells' latitude,
# longitude, row, column and observation time aren't read)
Variable = Literal[
    "albedo",
    "albedo_dca",
    "albedo_scah",
    "albedo_scav",
    "boresight_incidence",
    "bulk_density",
    "clay_fraction",
    "freeze_thaw_fraction",
    "grid_surface_status",
    "radar_water_body_fraction",
    "retrieval_qual_flag",
    "retrieval_qual_flag_dca",
    "retrieval_qual_flag_scah",
    "retrieval_qual_flag_scav",
    "roughness_coefficient",
    "roughness_coefficient_dca",
    "roughness_coefficient_scah",
    "roughness_coefficient_scav",
    "soil_moisture",
    "soil_moisture_dca",
    "soil_moisture_error",
    "soil_moisture_scah",
    "soil_moisture_scav",
    "static_water_body_fraction",
    "surface_flag",
    "surface_temperature",
    "surface_water_fraction_mb_h",
    "surface_water_fraction_mb_v",
    "tb_3_corrected",
    "tb_4_corrected",
    "tb_h_corrected",
    "tb_h_uncorrected",
    "tb_qual_flag_3",
    "tb_qual_flag_4",
    "tb_qual_flag_h",
    "tb_qual_flag_v",
    "tb_v_corrected",
    "tb_v_uncorrected",
    "vegetation_opacity",
    "vegetation_opacity_dca",
    "vegetation_opacity_scah",
    "vegetation_opacity_scav",
    "vegetation_water_content",
]
VARIABLES: tuple[Variable, ...] = (
    "soil_moisture",
    "soil_moisture_error",
    "tb_v_corrected",
    "tb_h_corrected",
    "retrieval_qual_flag",
    "tb_qual_flag_v",
    "tb_qual_flag_h",
    "surface_flag",
)

# The files' units in UDUNITS
UNITS = {"Kelvin": "K", "cm**3/cm**3": "m3 m-3"}
SOIL_MOISTURE = "volume_fraction_of_condensed_water_in_soil"
STANDARD_NAMES = {
    "soil_moisture": SOIL_MOISTURE,
    "soil_moisture_dca": SOIL_MOISTURE,
    "soil_moisture_scah": SOIL_MOISTURE,
    "soil_moisture_scav": SOIL_MOISTURE,
    "tb_v_corrected": "brightness_temperature",
    "tb_h_corrected": "brightness_temperature",
    "tb_v_uncorrected": "brightness_temperature",
    "tb_h_uncorrected": "brightness_temperature",
}

# Chunks of whole time steps, as many as fit in 10 MB (the size io stores as
# one chunk); a time step larger than that is split into bands of rows
CHUNK_BYTES = 10 * 1024**2

# Element names in a DMR++ file
DAP = "{http://xml.opendap.org/ns/DAP/4.0#}"


def window(aoi: gpd.GeoDataFrame) -> tuple[slice, slice]:
    """The rows and columns of the grid's cells *aoi*'s bounds touch, and one
    more on each side, within the grid: rioxarray needs two cells to find the
    grid's spacing, and ``io.clip_to_aoi`` drops the extra ones."""
    minx, miny, maxx, maxy = aoi.to_crs(CRS).total_bounds
    rows = slice(
        max(math.floor((Y0 - maxy) / CELL) - 1, 0),
        min(math.ceil((Y0 - miny) / CELL) + 1, SHAPE[0]),
    )
    cols = slice(
        max(math.floor((minx - X0) / CELL) - 1, 0),
        min(math.ceil((maxx - X0) / CELL) + 1, SHAPE[1]),
    )
    if rows.stop <= rows.start or cols.stop <= cols.start:
        raise ValueError("AOI is outside SMAP's grid, which reaches 85.044° N and S")
    return rows, cols


def coords(rows: slice, cols: slice) -> dict[str, np.ndarray]:
    """Cell-center ``y`` and ``x`` of *rows* and *cols*."""
    return {
        "y": Y0 - CELL * (np.arange(rows.start, rows.stop) + 0.5),
        "x": X0 + CELL * (np.arange(cols.start, cols.stop) + 0.5),
    }


def chunks(ny: int, nx: int) -> dict[str, int]:
    """Chunks of a float32 dataset of *ny* by *nx* cells (see CHUNK_BYTES)."""
    step = 4 * ny * nx
    if step <= CHUNK_BYTES:
        return {"time": CHUNK_BYTES // step, "y": -1, "x": -1}
    return {"time": 1, "y": CHUNK_BYTES // (4 * nx), "x": -1}


class SignedRanges:
    """Byte ranges of one file at a signed URL, for VirtualiZarr's chunk reads.

    VirtualiZarr reads a chunk with ``get_range_async(path, start=, end=)``; its
    ManifestStore drops the query from a chunk's URL, and obstore escapes one in
    a path, so a signed URL needs a store of its own. Each aligned
    :data:`BLOCK` is fetched once, and the reads share it.
    """

    def __init__(self, url: str):
        self.url = url
        self.blocks: dict[int, Future] = {}
        self.lock = threading.Lock()

    def _get(self, block: int) -> bytes:
        headers = {"Range": f"bytes={block * BLOCK}-{(block + 1) * BLOCK - 1}"}
        resp = fetch.get(self.url, headers=headers)
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


def open_arrays(url: str, signed: str, paths: list[str]) -> dict[str, zarr.Array]:
    """The arrays at *paths* of the file at *url*, read through its DMR++
    sidecar (``<url>.dmrpp``, fetched with the token) from its *signed* URL."""
    from virtualizarr.parsers.dmrpp import DMRParser

    resp = fetch.get(
        f"{url}.dmrpp", headers={"Authorization": f"Bearer {earthdata.token()}"}
    )
    resp.raise_for_status()
    root = ET.fromstring(resp.content)
    store = SignedRanges(signed)
    groups: dict[str, list[str]] = {}
    for path in paths:
        group, name = path.strip("/").split("/")
        groups.setdefault(group, []).append(name)
    out = {}
    for group, names in groups.items():
        # Only the arrays asked for: the parser skips the others' chunk lists
        others = [
            v.attrib["name"]
            for g in root.iter(f"{DAP}Group")
            if g.get("name") == group
            for v in g
            if "name" in v.attrib and v.attrib["name"] not in names
        ]
        parser = DMRParser(root=root, data_filepath=url, skip_variables=others)
        with warnings.catch_warnings():  # e.g. the string arrays it can't parse
            warnings.simplefilter("ignore")
            # SignedRanges has the one method the reads use
            manifest = parser.parse_dataset(object_store=cast(Any, store), group=group)
        z = zarr.open_group(manifest, mode="r", zarr_format=3)
        out.update({f"/{group}/{n}": cast(zarr.Array, z[n]) for n in names})
    return out


def fill_value(attrs: dict, dtype: np.dtype) -> float:
    """The array's ``_FillValue``: VirtualiZarr leaves the floats' as base64 of
    their bytes (its Zarr fill_value is 0)."""
    fill = attrs["_FillValue"]
    if isinstance(fill, str):
        raw = base64.b64decode(fill)
        fill = np.frombuffer(raw, "<f8" if len(raw) == 8 else dtype)[0]
    return float(fill)


def cf_attrs(attrs: dict) -> dict:
    """An array's attributes for the store: units in UDUNITS, and its bit
    flags' masks as numbers."""
    out = {
        k: attrs[k]
        for k in ("long_name", "valid_min", "valid_max", "flag_meanings")
        if k in attrs
    }
    if "flag_masks" in attrs:  # e.g. "1s, 2s, 4s": HDF5's short integers
        out["flag_masks"] = [
            int(m.strip().removesuffix("s")) for m in attrs["flag_masks"].split(",")
        ]
    if unit := attrs.get("units"):
        out["units"] = UNITS.get(unit, unit)
    return out


def describe(
    url: str, signed: str, paths: list[str]
) -> dict[str, tuple[np.dtype, float, dict]]:
    """The dtype, fill value and attributes of each array at *paths*."""
    out = {}
    for path, array in open_arrays(url, signed, paths).items():
        attrs = dict(array.attrs)
        out[path] = (array.dtype, fill_value(attrs, array.dtype), cf_attrs(attrs))
    return out


def read(
    url: str, signed: str, paths: list[str], rows: slice, cols: slice
) -> list[np.ndarray]:
    """*rows* and *cols* of each array at *paths* in the file at *url*: one
    DMR++ fetch, then the arrays' chunks, all at once."""
    arrays = open_arrays(url, signed, paths)
    with ThreadPoolExecutor(len(paths)) as pool:
        return list(pool.map(lambda p: np.asarray(arrays[p][rows, cols]), paths))


def file_day(item: pystac.Item) -> pd.Timestamp:
    """The UTC day of *item*'s file."""
    day = pd.Timestamp(item.properties["start_datetime"][:10])
    assert isinstance(day, pd.Timestamp)  # not NaT
    return day


def source_attrs(items: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every file."""
    attrs = provenance.from_stac_collection(STAC)
    attrs["description"] = (
        "SMAP enhanced L3 radiometer global daily 9 km EASE-Grid soil moisture "
        "(SPL3SMP_E), with its brightness temperatures, from NSIDC DAAC."
    )
    attrs["version"] = "SPL3SMP_E_006"
    attrs["processing:lineage"] = (
        "the newest processing of each day's file; its AM (descending) and PM "
        "(ascending, _pm) groups' cells clipped to the AOI on the native 9 km "
        "EASE-Grid 2.0; -9999 and days missing from the archive → NaN"
    )
    attrs["links"] += [
        provenance.link(
            "derived_from", f"{STAC}/items/{it.id}", type="application/geo+json"
        )
        for it in items
    ]
    return attrs


def load(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str,
    variables: tuple[Variable, ...] = VARIABLES,
    orbit_state: OrbitState = OrbitState.any,
) -> xr.Dataset:
    """SMAP SPL3SMP_E for each day over *aoi*, lazily, on the native 9 km grid.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    variables : list of str
        Grids of each pass, by their AM names: soil moisture (soil_moisture,
        its error and the single-channel retrievals' _dca, _scah, _scav),
        brightness temperature (tb_v_corrected, tb_h_corrected: water bodies
        removed), the bit flags (retrieval_qual_flag, tb_qual_flag_v,
        surface_flag, ...) and the retrievals' inputs (surface_temperature,
        vegetation_opacity, ...).
    orbit_state : OrbitState
        Passes to read: any (both), descending (6 AM) or ascending (6 PM, the
        variables named with _pm).

    Returns
    -------
    xarray.Dataset
        One variable per grid and pass along ``(time, y, x)``, one time step
        per day, plus CF ``time_bounds``: float32 grids with -9999 as NaN, and
        uint16 bit flags with 65534 as their fill value.
    """
    orbit_state = OrbitState(orbit_state)
    passes = list(GROUPS) if orbit_state is OrbitState.any else [orbit_state]
    # Each variable's array in the files, e.g. soil_moisture_pm in the PM group
    paths = {
        f"{name}{suffix}": f"/{group}/{name}{suffix}"
        for group, suffix in (GROUPS[p] for p in passes)
        for name in variables
    }
    rows, cols = window(aoi)
    days = pd.date_range(start_date, end_date, freq="D", normalize=True)

    records = earthdata.cmr_search(CONCEPT_ID, aoi, start_date, end_date)
    found = earthdata.newest([earthdata.cmr_item(r) for r in records], key=file_day)
    if older := len(records) - len(found):
        log.info("SMAP: skipped %d older processing(s)", older)
    files = {file_day(it): it for it in found}
    items = [files[day] for day in days if day in files]
    if not items:
        raise ValueError(f"No SMAP SPL3SMP_E files for {start_date} to {end_date}")
    if missing := [day for day in days if day not in files]:
        log.warning(
            "SMAP: %d of %d day(s) are missing from the archive (first %s); "
            "filled with NaN",
            len(missing),
            len(days),
            f"{missing[0]:%Y-%m-%d}",
        )
    # Each file's URL (its .dmrpp is next to it) and the URL signed for its chunks
    urls = {file_day(it): it.assets["data"].href for it in items}
    earthdata.sign(items, ["data"])
    signed = {file_day(it): it.assets["data"].href for it in items}
    first = file_day(items[0])
    meta = describe(urls[first], signed[first], list(paths.values()))

    # One read of each day's file gives all its arrays
    reads = {
        day: dask.delayed(read)(
            urls[day], signed[day], list(paths.values()), rows, cols
        )
        for day in urls
    }
    shape = (rows.stop - rows.start, cols.stop - cols.start)
    parts = []
    for i, (name, path) in enumerate(paths.items()):
        dtype, fill, attrs = meta[path]
        steps = [
            dask.array.from_delayed(reads[day][i], shape, dtype=dtype)
            if day in reads
            else dask.array.full(shape, fill, dtype=dtype)
            for day in days
        ]
        da = xr.DataArray(
            dask.array.stack(steps),
            dims=("time", "y", "x"),
            coords={"time": days, **coords(rows, cols)},
            name=name,
        ).rio.write_crs(CRS)
        floats = np.issubdtype(dtype, np.floating)
        if floats:
            da = da.where(da != fill)
        # The clip fills the cells outside the AOI with the nodata value
        nodata = np.nan if floats else dtype.type(fill)
        da = io.clip_to_aoi(da.rio.write_nodata(nodata), aoi)
        # Replaced, not merged: rioxarray's nodata is a _FillValue attribute
        da.attrs = dict(attrs)
        if not floats:
            da.encoding["_FillValue"] = int(fill)
            # TODO: teach `gaia align` CF flag_masks (bit flags), as it does flag_values
            da.attrs[regrid.METHOD_ATTR] = "nearest"
        if standard_name := STANDARD_NAMES.get(name.removesuffix("_pm")):
            da.attrs["standard_name"] = standard_name
        # Each value is one snapshot within its UTC day
        parts.append(io.with_time_bounds(da, 0, 24, "time: point"))
    # The parts share their time_bounds and grid
    ds = xr.merge(parts, compat="override")
    ds = ds.chunk(chunks(ds.sizes["y"], ds.sizes["x"]))
    log.info(
        "SMAP: %d day(s) (%d file(s)), %d x %d cells of 9 km",
        len(days),
        len(items),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    ds.attrs = {**source_attrs(items), **provenance.time_extent(ds)}
    return ds
