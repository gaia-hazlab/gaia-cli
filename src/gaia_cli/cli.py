"""
A CLI following https://packaging.python.org/en/latest/guides/creating-command-line-tools/
"""

import logging
import os
import sys

import cyclopts

from .align import run as align
from .all import stage as stage_all
from .catalog import rebuild as catalog_rebuild
from .hrrr import stage as stage_hrrr
from .nlcd import stage as stage_nlcd
from .prism import fill as prism_fill
from .prism import rebuild_catalog as prism_rebuild_catalog
from .prism import stage as stage_prism
from .provenance import set_command
from .stagefile import stage as stage_config
from .synoptic import stage as stage_synoptic

app = cyclopts.App()

# gaia stage --config stage.yaml | gaia stage {synoptic,prism,hrrr,all,nlcd} ...
stage_app = cyclopts.App(
    name="stage", help="Stage datasets clipped to an AOI and time range."
)
app.command(stage_app)
stage_app.default(stage_config)

stage_app.command(stage_synoptic, name="synoptic")
stage_app.command(stage_prism, name="prism")
stage_app.command(stage_hrrr, name="hrrr")
stage_app.command(stage_all, name="all")
stage_app.command(stage_nlcd, name="nlcd")

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
