"""One module per data source, registered in :data:`DATASETS`.

Each module has ``load()``, ``EXTENT`` (a STAC Collection extent) and ``KIND``,
and optionally ``VARIANT``: a ``load()`` option named by the key
(``dem/usgs10m``), each value its own top-level group; and ``REGRID_METHOD``:
``gaia align``'s default for its raster variables (else nearest); and ``SCENES``:
it takes ``pre``/``post`` windows and ``select`` to choose scenes, and ``gaia
align`` keeps its acquisition times; and ``TOKENS``: environment variables it
needs, which ``gaia stage`` warns about up front. See
dev/adding-a-dataset.md.
"""

from types import ModuleType

from . import (
    burn_severity,
    dem,
    era5,
    hrrr,
    nlcd,
    prism,
    sentinel1_rtc,
    sentinel2,
    soil,
    stage4,
    synoptic,
)

# Name -> module. The name is the `gaia stage <name>` command, the stage-file
# key and the store group's last part
DATASETS: dict[str, ModuleType] = {
    "synoptic": synoptic,
    "prism": prism,
    "hrrr": hrrr,
    "stage4": stage4,
    "era5": era5,
    "nlcd": nlcd,
    "dem": dem,
    "soil": soil,
    "burn-severity": burn_severity,
    "sentinel1-rtc": sentinel1_rtc,
    "sentinel2": sentinel2,
}


def cf_criteria() -> dict[str, dict[str, str]]:
    """cf-xarray custom criteria selecting variables by their common name.

    Each dataset's ``COMMON_NAMES`` is stamped on its staged variables as
    ``gaia:common-name``; with these criteria ``ds.cf["precipitation"]`` gives PRISM
    ``ppt``, ERA5 ``tp`` or HRRR ``precipitation_surface``, staged or aligned::

        import cf_xarray as cfxr
        cfxr.set_options(custom_criteria=cf_criteria())
    """
    from ..provenance import COMMON_NAME_ATTR

    names = {
        c for m in DATASETS.values() for c in getattr(m, "COMMON_NAMES", {}).values()
    }
    return {c: {COMMON_NAME_ATTR: f"^{c}$"} for c in sorted(names)}
