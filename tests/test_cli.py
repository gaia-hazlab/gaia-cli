from gaia_cli.cli import stage_app


def test_stage_subcommands_registered():
    names = set(stage_app._commands)
    assert {"synoptic", "prism", "hrrr", "all", "nlcd"} <= names
