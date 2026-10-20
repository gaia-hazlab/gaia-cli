"""
A CLI following https://packaging.python.org/en/latest/guides/creating-command-line-tools/
"""

import logging
import os
import sys

import cyclopts

from .align import run as align
from .catalog import rebuild as catalog_rebuild
from .datasets import DATASETS
from .datasets.prism import fill as prism_fill
from .datasets.prism import rebuild_catalog as prism_rebuild_catalog
from .provenance import set_command
from .runner import dataset_command
from .runner import stage as stage_config

app = cyclopts.App(name="gaia")

# gaia stage --config {stage.yaml,PRESET} | gaia stage <dataset> ...
stage_app = cyclopts.App(
    name="stage", help="Stage datasets clipped to an AOI and time range."
)
app.command(stage_app)
stage_app.default(stage_config)

for name in DATASETS:
    stage_app.command(dataset_command(name), name=name)

# gaia align ...
app.command(align, name="align")

# gaia catalog rebuild ...
catalog_app = cyclopts.App(name="catalog", help="STAC catalog of staged stores.")
app.command(catalog_app)
catalog_app.command(catalog_rebuild, name="rebuild")

# gaia prism {fill,rebuild-catalog} ...
prism_app = cyclopts.App(name="prism", help="Manage the gaia PRISM cache.")
app.command(prism_app)
prism_app.command(prism_fill, name="fill")
prism_app.command(prism_rebuild_catalog, name="rebuild-catalog")


class _Formatter(logging.Formatter):
    """``12:03:04 [  12.3s] prism: message``: wall clock plus time since start."""

    def format(self, record):
        record.elapsed = record.relativeCreated / 1000
        record.short = record.name.removeprefix("gaia_cli.")
        if record.levelno != logging.INFO:
            record.short = f"{record.levelname} {record.short}"
        return super().format(record)


def configure_logging():
    """Log gaia_cli messages to stderr; level from ``GAIA_LOG_LEVEL`` (default INFO).

    Other libraries only log warnings and errors.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(
        _Formatter("%(asctime)s [%(elapsed)6.1fs] %(short)s: %(message)s", "%H:%M:%S")
    )
    logging.basicConfig(level=logging.WARNING, handlers=[handler], force=True)
    logging.getLogger("gaia_cli").setLevel(
        os.environ.get("GAIA_LOG_LEVEL", "INFO").upper()
    )
    logging.captureWarnings(True)


def main():
    configure_logging()
    set_command(sys.argv[1:])
    app()


if __name__ == "__main__":
    main()
