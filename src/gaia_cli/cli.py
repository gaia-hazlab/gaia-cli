"""
A CLI following https://packaging.python.org/en/latest/guides/creating-command-line-tools/
"""

import cyclopts

from .synoptic import stage as stage_synoptic
from .prism import stage as stage_prism
from .hrrr import stage as stage_hrrr
from .all import stage as stage_all

app = cyclopts.App()

# gaia stage {synoptic,prism,hrrr,all} ...
stage_app = cyclopts.App(
    name="stage", help="Stage datasets clipped to an AOI and time range."
)
app.command(stage_app)

stage_app.command(stage_synoptic, name="synoptic")
stage_app.command(stage_prism, name="prism")
stage_app.command(stage_hrrr, name="hrrr")
stage_app.command(stage_all, name="all")

if __name__ == "__main__":
    app()
