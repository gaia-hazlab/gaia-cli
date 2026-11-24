"""SMAP L-band brightness temperature and soil moisture (9 km), clipped to an AOI.

Each day in the window is a time step of NSIDC's enhanced Level-3 radiometer
product (SPL3SMP_E version 006): the day's morning (descending) and evening
(ascending) passes on the 9 km EASE-Grid 2.0, the evening's variables with the
files' ``_pm`` suffix. Days missing from the archive are NaN, with a warning.
Read remotely from each day's HDF5 file through its DMR++ index (``dmrpp``).
Needs EARTHDATA_TOKEN.

Examples
--------
    $ gaia stage smap -i palouse.geojson -s 2025-09-01 -e 2025-09-10 -o palouse-smap
    $ gaia stage smap -i palouse.geojson -s 2025-09-01 -e 2025-09-10 -o palouse-smap --orbit-state descending --variables soil_moisture
"""

from __future__ import annotations

import logging
import math
from collections.abc import Callable, Mapping, Sequence
from typing import Literal, get_args

import dask
import dask.array
import geopandas as gpd
import numpy as np
import pandas as pd
import pystac
import xarray as xr
import zarr

from .. import dmrpp, earthdata, io, provenance, regrid
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
# Every array in a pass's group, by its AM name: the grids and the cells'
# coordinates and times (the reader skips the ones not asked for)
ARRAYS = (
    *get_args(Variable),
    "latitude",
    "longitude",
    "latitude_centroid",
    "longitude_centroid",
    "EASE_row_index",
    "EASE_column_index",
    "tb_time_seconds",
    "tb_time_utc",
)
GROUP_ARRAYS = {
    group: [f"{a}{suffix}" for a in ARRAYS] for group, suffix in GROUPS.values()
}

# The files' units in UDUNITS
UNITS = {"Kelvin": "K", "cm**3/cm**3": "m3 m-3", "dimensionless": "1"}
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
# Files signed together, in one task their reads wait on: one round of
# earthdata.signed_urls' 16 requests, read well within the signed URLs' hour
SIGN_BATCH = 16


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
    arrays: Mapping[str, zarr.Array],
) -> dict[str, tuple[np.dtype, float, dict]]:
    """The dtype, fill value and attributes of each of *arrays*."""
    out = {}
    for path, array in arrays.items():
        attrs = dict(array.attrs)
        fill = dmrpp.fill_value(attrs, array.dtype)
        out[path] = (array.dtype, fill, cf_attrs(attrs))
    return out


def find_files(
    concept_id: str,
    product: str,
    aoi: gpd.GeoDataFrame,
    start: str,
    end: str,
    steps: pd.DatetimeIndex,
    key: Callable[[pystac.Item], pd.Timestamp],
    suffix: str = "",
) -> dict[pd.Timestamp, pystac.Item]:
    """The newest processing of the file of each of *steps* (its time step:
    *key*) in the collection *concept_id*, searched from *start* to *end*.
    Steps missing from the archive get a warning; none at all is an error."""
    records = earthdata.cmr_search(concept_id, aoi, start, end)
    items = [earthdata.cmr_item(r, suffix) for r in records]
    found = earthdata.newest(items, key=key)
    if older := len(records) - len(found):
        log.info("SMAP %s: skipped %d older processing(s)", product, older)
    # The search also finds files that only touch the window
    files = {key(it): it for it in found if key(it) in steps}
    if not files:
        raise ValueError(f"No SMAP {product} files for {start} to {end}")
    if missing := [step for step in steps if step not in files]:
        log.warning(
            "SMAP %s: %d of %d time step(s) are missing from the archive "
            "(first %s); filled with NaN",
            product,
            len(missing),
            len(steps),
            missing[0].isoformat(),
        )
    return dict(sorted(files.items()))


def stack(
    files: Mapping[pd.Timestamp, pystac.Item],
    steps: pd.DatetimeIndex,
    paths: Mapping[str, str],
    arrays: Mapping[str, Sequence[str]],
    aoi: gpd.GeoDataFrame,
    standard_names: Mapping[str, str],
    time_bounds: tuple[float, float, str],
) -> xr.Dataset:
    """The arrays at *paths* (by variable name) of each of *files* over *aoi*,
    lazily, along ``(time, y, x)``: one time step per file, and a step of
    NaN (an integer array's fill value) for each of *steps* without one.
    *arrays* names every array in the files' groups (``dmrpp.open_arrays``),
    and *time_bounds* is ``io.with_time_bounds``' hours and cell methods."""
    rows, cols = window(aoi)
    urls = {step: it.assets["data"].href for step, it in files.items()}
    # The first file, opened once: its attributes, then its read
    first = next(iter(urls))
    first_arrays = dmrpp.open_arrays(urls[first], list(paths.values()), arrays)
    meta = describe(first_arrays)

    # The others are signed in batches, each just before its files' reads (the
    # first is signed when its first chunk is read)
    rest = [url for step, url in urls.items() if step != first]
    signed = {}
    for i in range(0, len(rest), SIGN_BATCH):
        batch = rest[i : i + SIGN_BATCH]
        hrefs = dask.delayed(earthdata.signed_urls)(batch)
        signed.update({url: hrefs[url] for url in batch})

    # One read of each file gives all its arrays
    reads = {
        step: dask.delayed(dmrpp.read)(first_arrays, list(paths.values()), rows, cols)
        if step == first
        else dask.delayed(dmrpp.read_file)(
            url, list(paths.values()), arrays, rows, cols, signed[url]
        )
        for step, url in urls.items()
    }
    shape = (rows.stop - rows.start, cols.stop - cols.start)
    parts = []
    for i, (name, path) in enumerate(paths.items()):
        dtype, fill, attrs = meta[path]
        values = [
            dask.array.from_delayed(reads[step][i], shape, dtype=dtype)
            if step in reads
            else dask.array.full(shape, fill, dtype=dtype)
            for step in steps
        ]
        da = xr.DataArray(
            dask.array.stack(values),
            dims=("time", "y", "x"),
            coords={"time": steps, **coords(rows, cols)},
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
        if standard_name := standard_names.get(name):
            da.attrs["standard_name"] = standard_name
        parts.append(io.with_time_bounds(da, *time_bounds))
    # The parts share their time_bounds and grid
    ds = xr.merge(parts, compat="override")
    return ds.chunk(chunks(ds.sizes["y"], ds.sizes["x"]))


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
    days = pd.date_range(start_date, end_date, freq="D", normalize=True)
    files = find_files(
        CONCEPT_ID, "SPL3SMP_E", aoi, start_date, end_date, days, file_day
    )
    standard_names = {
        name: sn
        for name in paths
        if (sn := STANDARD_NAMES.get(name.removesuffix("_pm")))
    }
    # Each value is one snapshot within its UTC day
    ds = stack(
        files, days, paths, GROUP_ARRAYS, aoi, standard_names, (0, 24, "time: point")
    )
    log.info(
        "SMAP: %d day(s) (%d file(s)), %d x %d cells of 9 km",
        len(days),
        len(files),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    ds.attrs = {**source_attrs(list(files.values())), **provenance.time_extent(ds)}
    return ds
