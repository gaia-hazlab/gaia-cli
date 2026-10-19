from gaia_cli.cli import app, stage_app


def test_stage_subcommands_registered():
    names = set(stage_app._commands)
    assert {"synoptic", "prism", "hrrr", "nlcd", "dem"} <= names
    assert "all" not in names  # replaced by the conus-precip preset


def test_align_command_registered():
    assert "align" in app._commands
