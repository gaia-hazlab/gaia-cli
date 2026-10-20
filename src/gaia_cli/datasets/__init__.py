"""One module per data source, registered in :data:`DATASETS`.

Each module has ``load()``, ``EXTENT`` (a STAC Collection extent) and ``KIND``,
and optionally ``VARIANT``: a ``load()`` option named by the key
(``dem/usgs10m``), each value its own top-level group; and ``REGRID_METHOD``:
``gaia align``'s default for its raster variables (else nearest). See
dev/adding-a-dataset.md.
"""

from types import ModuleType

from . import burn_severity, dem, hrrr, nlcd, prism, soil, synoptic

# Name -> module. The name is the `gaia stage <name>` command, the stage-file
# key and the store group's last part
DATASETS: dict[str, ModuleType] = {
    "synoptic": synoptic,
    "prism": prism,
    "hrrr": hrrr,
    "nlcd": nlcd,
    "dem": dem,
    "soil": soil,
    "burn-severity": burn_severity,
}
