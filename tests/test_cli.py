import yaml

from gaia_cli import stagefile
from gaia_cli.cli import app, stage_app


def test_stage_subcommands_registered():
    assert set(stagefile.DATASETS) <= set(stage_app._commands)


def test_align_command_registered():
    assert "align" in app._commands


def test_dataset_commands_take_load_options_as_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    stagefile.load_aoi({"bbox": [-121.6, 47.6, -121.0, 47.9]}).to_file(aoi)
    common = ["-i", str(aoi), "-o", str(tmp_path / "out"), "--dry-run"]
    window = ["-s", "2025-12-01", "-e", "2025-12-02"]

    def dry_run(*args):
        stage_app([*args, *common], result_action="return_value")
        return yaml.safe_load(capsys.readouterr().out)["datasets"]

    prism = dry_run("prism", *window, "-v", "ppt", "tmin")["prism"]
    assert prism["variables"] == ["ppt", "tmin"]
    assert prism["window"] == "2025-12-01T00:00:00/2025-12-03T00:00:00"
    assert list(dry_run("dem", "-p", "usgs10m")) == ["dem/usgs10m"]
