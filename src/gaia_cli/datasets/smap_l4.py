"""SMAP Level-4 surface and root-zone soil moisture (9 km, 3-hourly), clipped to an AOI.

Each time step is one of NSIDC's SPL4SMGP version 008 files: 3-hour means of the
L4_SM land model, into which SMAP's brightness temperatures are assimilated, on
``smap``'s 9 km EASE-Grid 2.0. Values are stamped at the end of their 3 hours.
Steps missing from the archive are NaN, with a warning. Each file stores each
variable as one whole-globe chunk, so any AOI reads all of it: about 4.3 MB per
soil moisture variable and step. Read remotely through each file's DMR++ index
(``dmrpp``). Needs EARTHDATA_TOKEN.

Examples
--------
    $ gaia stage smap-l4 -i palouse.geojson -s 2025-09-01 -e 2025-09-03 -o palouse-smap-l4
    $ gaia stage smap-l4 -i palouse.geojson -s 2025-09-01 -e 2025-09-03 -o palouse-smap-l4 --variables sm_surface snow_mass
"""

from __future__ import annotations

import logging
from typing import Literal, get_args

import geopandas as gpd
import pandas as pd
import pystac
import xarray as xr

from .. import earthdata, provenance
from . import smap

log = logging.getLogger(__name__)

# SPL4SMGP version 008 in NSIDC's cloud, pinned
CONCEPT_ID = "C3480440870-NSIDC_CPRD"
STAC = "https://cmr.earthdata.nasa.gov/stac/NSIDC_CPRD/collections/SPL4SMGP_008"
# The global grid's, as smap's
EXTENT = {
    "spatial": {"bbox": [[-180.0, -85.044, 180.0, 85.044]]},
    "temporal": {"interval": [["2015-03-31T00:00:00Z", None]]},
}
KIND = "raster"
# `gaia align` default
REGRID_METHOD = "average"
TOKENS = earthdata.TOKENS

GROUP = "Geophysical_Data"
STEP = pd.Timedelta(hours=3)

# The grids of the files' Geophysical_Data group (the root's cell coordinates,
# times and grid mapping aren't read)
Variable = Literal[
    "baseflow_flux",
    "depth_to_water_table_from_surface",
    "free_surface_water_on_peat_flux",
    "heat_flux_ground",
    "heat_flux_latent",
    "heat_flux_sensible",
    "height_lowatmmodlay",
    "land_evapotranspiration_flux",
    "land_fraction_saturated",
    "land_fraction_snow_covered",
    "land_fraction_unsaturated",
    "land_fraction_wilting",
    "leaf_area_index",
    "mwrtm_vegopacity",
    "net_downward_longwave_flux",
    "net_downward_shortwave_flux",
    "overland_runoff_flux",
    "precipitation_total_surface_flux",
    "radiation_longwave_absorbed_flux",
    "radiation_shortwave_downward_flux",
    "sm_profile",
    "sm_profile_pctl",
    "sm_profile_wetness",
    "sm_rootzone",
    "sm_rootzone_pctl",
    "sm_rootzone_wetness",
    "sm_surface",
    "sm_surface_wetness",
    "snow_depth",
    "snow_mass",
    "snow_melt_flux",
    "snowfall_surface_flux",
    "soil_temp_layer1",
    "soil_temp_layer2",
    "soil_temp_layer3",
    "soil_temp_layer4",
    "soil_temp_layer5",
    "soil_temp_layer6",
    "soil_water_infiltration_flux",
    "specific_humidity_lowatmmodlay",
    "surface_pressure",
    "surface_temp",
    "temp_lowatmmodlay",
    "vegetation_greenness_fraction",
    "windspeed_lowatmmodlay",
]
VARIABLES: tuple[Variable, ...] = ("sm_surface", "sm_rootzone", "sm_profile")
ARRAYS = {GROUP: list(get_args(Variable))}

# CF standard names carry no depth, so the layers sharing one (sm_rootzone and
# sm_profile, soil_temp_layer1 ... 6) have none: their long_name gives it.
# sm_surface's 0-5 cm is smap's soil_moisture's
STANDARD_NAMES = {
    "sm_surface": smap.SOIL_MOISTURE,
    "surface_temp": "surface_temperature",
    "snow_mass": "surface_snow_amount",
    "snow_depth": "surface_snow_thickness",
    "precipitation_total_surface_flux": "precipitation_flux",
    "snowfall_surface_flux": "snowfall_flux",
    "overland_runoff_flux": "surface_runoff_flux",
    "land_evapotranspiration_flux": "water_evapotranspiration_flux",
    "heat_flux_latent": "surface_upward_latent_heat_flux",
    "heat_flux_sensible": "surface_upward_sensible_heat_flux",
    "radiation_shortwave_downward_flux": "surface_downwelling_shortwave_flux_in_air",
    "surface_pressure": "surface_air_pressure",
    "leaf_area_index": "leaf_area_index",
    "temp_lowatmmodlay": "air_temperature",
    "specific_humidity_lowatmmodlay": "specific_humidity",
    "windspeed_lowatmmodlay": "wind_speed",
}


def file_step(item: pystac.Item) -> pd.Timestamp:
    """The end of *item*'s 3 hours, its time step (naive UTC)."""
    end = pd.Timestamp(item.properties["end_datetime"])
    assert isinstance(end, pd.Timestamp)  # not NaT
    return end.tz_convert(None)


def steps(start_date: str, end_date: str) -> pd.DatetimeIndex:
    """The 3-hourly time steps from *start_date* to *end_date*, inclusive; a
    date-only end includes its whole day, as ``sel(time=slice(...))`` does."""
    first = pd.Timestamp(start_date).floor("D")
    last = pd.Timestamp(end_date).floor("D") + pd.Timedelta(days=1)
    grid = pd.date_range(first, last, freq=STEP)
    return grid[grid.slice_indexer(start_date, end_date)]


def source_attrs(items: list[pystac.Item]) -> dict:
    """Provenance attributes: the collection's, plus a link to every file."""
    attrs = provenance.from_stac_collection(STAC)
    attrs["description"] = (
        "SMAP L4 global 3-hourly 9 km EASE-Grid surface and root zone soil "
        "moisture geophysical data (SPL4SMGP), from NSIDC DAAC."
    )
    attrs["version"] = "SPL4SMGP_008"
    attrs["processing:lineage"] = (
        "the newest processing of each 3-hour file; its Geophysical_Data "
        "group's cells clipped to the AOI on the native 9 km EASE-Grid 2.0; "
        "-9999 and steps missing from the archive → NaN"
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
) -> xr.Dataset:
    """SMAP SPL4SMGP for each 3 hours over *aoi*, lazily, on the native 9 km grid.

    Parameters
    ----------
    aoi : geopandas.GeoDataFrame
        Area of interest polygon.
    start_date : str
        Start date (inclusive).
    end_date : str
        End date (inclusive).
    variables : list of str
        Grids of the files' Geophysical_Data group: soil moisture (sm_surface,
        0-5 cm; sm_rootzone, 0-100 cm; sm_profile, to the model's bedrock;
        each also as _wetness, relative saturation, and the latter two as
        _pctl, percentiles), soil and surface
        temperature (soil_temp_layer1 ... 6, surface_temp), snow (snow_mass,
        snow_depth), water and energy fluxes (precipitation_total_surface_flux,
        land_evapotranspiration_flux, overland_runoff_flux, heat_flux_latent,
        ...) and the model's forcing (temp_lowatmmodlay, surface_pressure, ...).

    Returns
    -------
    xarray.Dataset
        One variable per grid along ``(time, y, x)``, one time step per 3
        hours, stamped at their end, plus CF ``time_bounds``: float32 means
        with -9999 as NaN.
    """
    paths: dict[str, str] = {name: f"/{GROUP}/{name}" for name in variables}
    times = steps(start_date, end_date)
    # The file of the first step covers the 3 hours before it
    files = smap.find_files(
        CONCEPT_ID,
        "SPL4SMGP",
        aoi,
        (times[0] - STEP).isoformat(),
        times[-1].isoformat(),
        times,
        file_step,
        suffix=".h5",
    )
    # Each value is the mean of the 3 hours ending at its time step
    ds = smap.stack(
        files, times, paths, ARRAYS, aoi, STANDARD_NAMES, (-3, 0, "time: mean")
    )
    log.info(
        "SMAP L4: %d step(s) (%d file(s)), %d x %d cells of 9 km",
        len(times),
        len(files),
        ds.sizes["y"],
        ds.sizes["x"],
    )
    ds.attrs = {**source_attrs(list(files.values())), **provenance.time_extent(ds)}
    return ds
