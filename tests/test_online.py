"""Network smoke tests: stage each dataset for a small Bolt Creek AOI through the
stage-file runner. Run with ``pixi run test-online``; skipped by default."""

import os

import numpy as np
import pytest
import xarray as xr

from gaia_cli import stagefile, stations

pytestmark = pytest.mark.online

# Around the 2022 Bolt Creek fire and Stevens Pass, WA (~45 x 33 km).
# The padded read, 2025-12-09 to 2025-12-11, is already in the PRISM cache, so no write
# credentials are needed.
STAGE = """
name: bolt-creek-smoke
aoi: {{bbox: [-121.6, 47.6, -121.0, 47.9]}}
window: 2025-12-10/2025-12-10
datasets:
  {dataset}: {{}}
"""


# Each window is read one day beyond each end (stagefile.PAD)
def check_synoptic(ds):
    assert ds.sizes["geometry"] > 0
    assert int(ds.row_size.sum()) == ds.sizes["obs"]
    assert stations.to_grid(ds).sizes["time"] >= 72


def check_prism(ds):
    assert ds.sizes["time"] == 3
    assert float(ds.precipitation.min()) >= 0


def check_hrrr(ds):
    assert ds.sizes["time"] == 72
    assert float(ds.precipitation.min()) >= 0


def check_nlcd(ds):
    classes = np.unique(ds.landcover.values)
    assert {42, 52} & set(classes.tolist())  # evergreen forest, shrub/scrub


@pytest.mark.parametrize(
    "dataset, check",
    [
        pytest.param(
            "synoptic",
            check_synoptic,
            marks=pytest.mark.skipif(
                not os.environ.get("SYNOPTIC_TOKEN"), reason="needs SYNOPTIC_TOKEN"
            ),
        ),
        ("prism", check_prism),
        ("hrrr", check_hrrr),
        ("nlcd", check_nlcd),
    ],
)
def test_stage(tmp_path, dataset, check):
    config = tmp_path / "stage.yaml"
    config.write_text(STAGE.format(dataset=dataset))
    out = tmp_path / "out"
    stagefile.stage(config=str(config), output_path=str(out))

    tree = xr.open_datatree(out / "bolt-creek-smoke.zarr", engine="zarr")
    (path,) = [p for p in tree.groups if p.count("/") == 2]
    assert path.endswith(dataset)
    check(tree[path].to_dataset())
    assert tree[path].attrs["license"]
