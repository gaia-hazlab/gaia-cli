"""Network smoke tests: stage each dataset for a small Bolt Creek AOI through the
stage-file runner. Run with ``pixi run test-online``; skipped by default."""

import json
import os

import numpy as np
import pytest
import xarray as xr

from gaia_cli import align, stagefile, stations

pytestmark = pytest.mark.online

# Around the 2022 Bolt Creek fire and Stevens Pass, WA (~45 x 33 km).
# The padded read, 2025-12-09 to 2025-12-11, is already in the PRISM cache (ppt, tmin,
# tmax), so no write credentials are needed.
STAGE = """
name: bolt-creek-smoke
aoi: {{bbox: [-121.6, 47.6, -121.0, 47.9]}}
window: 2025-12-10/2025-12-10
datasets:
  {dataset}: {options}
"""


# Each window is read one day beyond each end (stagefile.PAD)
def check_synoptic(ds):
    assert ds.sizes["geometry"] > 0
    assert int(ds.row_size.sum()) == ds.sizes["obs"]
    assert stations.to_grid(ds).sizes["time"] >= 72


def check_prism(ds):
    assert ds.sizes["time"] == 3
    assert float(ds.precipitation.min()) >= 0


def check_prism_temperature(ds):
    assert ds.sizes["time"] == 3
    assert bool((ds.tmin <= ds.tmax).all())
    assert -40 < float(ds.tmin.min()) and float(ds.tmax.max()) < 30  # December, degC


def check_hrrr(ds):
    assert ds.sizes["time"] == 72
    assert float(ds.precipitation.min()) >= 0


def check_nlcd(ds):
    classes = np.unique(ds.landcover.values)
    assert {42, 52} & set(classes.tolist())  # evergreen forest, shrub/scrub


def check_dem(ds):
    # Skykomish valley (~150 m) to the Cascade crest (~2000 m)
    assert 50 < float(ds.elevation.min()) < 300
    assert 1800 < float(ds.elevation.max()) < 2500
    assert ds.elevation.attrs["vertical_datum"].startswith("EPSG:")


def check_soil(ds):
    # Decoded to physical units: acidic forest soils, bulk density near 1 g cm-3
    assert 3 < float(ds.ph1to1h2o.min()) < float(ds.ph1to1h2o.max()) < 9
    assert 0.3 < float(ds.dbovendry.min()) < float(ds.dbovendry.max()) < 2
    total = ds.claytotal + ds.sandtotal + ds.silttotal
    assert 90 < float(total.mean()) < 110  # percent
    assert ds.ph1to1h2o.attrs["gaia:regrid-method"] == "bilinear"


def check_burn_severity(ds):
    # The 2022 Bolt Creek fire: every class, NaN outside the mapped area
    classes = np.unique(ds.burn_severity.values[~np.isnan(ds.burn_severity.values)])
    assert classes.tolist() == [1, 2, 3, 4]
    assert ds.attrs["fire_name"] == "bolt_creek_sbs"


@pytest.mark.parametrize(
    "dataset, check, options",
    [
        pytest.param(
            "synoptic",
            check_synoptic,
            {},
            marks=pytest.mark.skipif(
                not os.environ.get("SYNOPTIC_TOKEN"), reason="needs SYNOPTIC_TOKEN"
            ),
        ),
        ("prism", check_prism, {}),
        pytest.param(
            "prism",
            check_prism_temperature,
            {"variables": ["ppt", "tmin", "tmax"]},
            id="prism-temperature",
        ),
        ("hrrr", check_hrrr, {}),
        ("nlcd", check_nlcd, {}),
        ("dem/cop30", check_dem, {}),
        ("dem/usgs10m", check_dem, {}),
        ("soil/solus100", check_soil, {}),
        ("burn_severity", check_burn_severity, {"fire": "bolt creek"}),
    ],
)
def test_stage(tmp_path, dataset, check, options):
    config = tmp_path / "stage.yaml"
    config.write_text(STAGE.format(dataset=dataset, options=json.dumps(options)))
    out = tmp_path / "out"
    stagefile.stage(config=str(config), output_path=str(out))

    tree = xr.open_datatree(out / "bolt-creek-smoke.zarr", engine="zarr")
    (path,) = [p for p in tree.groups if p.count("/") == 2]
    assert path.endswith(dataset)
    check(tree[path].to_dataset())
    assert tree[path].attrs["license"]


# The landlab_data_prep reference event: Pioneer Fire, WA (WA4818312053120240608),
# the perimeter's bbox. Landlab puts every layer on a 10 m UTM grid.
PIONEER = """
name: pioneer
aoi: {bbox: [-120.6807, 48.0978, -120.4009, 48.3895]}
datasets:
  dem/cop30: {}
  dem/usgs10m: {}
"""


def test_pioneer_dems_on_the_landlab_grid(tmp_path):
    config = tmp_path / "stage.yaml"
    config.write_text(PIONEER)
    out = tmp_path / "out"
    stagefile.stage(config=str(config), output_path=str(out))

    tree = xr.open_datatree(out / "pioneer.zarr", engine="zarr")
    # The DEM's default method (bilinear), as gaia align would use it
    aligned = align.align(tree, crs="EPSG:32610", res=10.0)
    cop30, usgs10m = (aligned[f"dem/{p}"].elevation for p in ("cop30", "usgs10m"))
    assert cop30.attrs["regrid_method"] == usgs10m.attrs["regrid_method"] == "bilinear"
    assert cop30.rio.crs.to_epsg() == usgs10m.rio.crs.to_epsg() == 32610
    assert cop30.shape == usgs10m.shape == (3311, 2185)
    # Lake Chelan's shore (~335 m) to the Chelan Mountains (~2680 m)
    for da in (cop30, usgs10m):
        assert 320 < float(da.min()) < 350
        assert 2600 < float(da.max()) < 2750
    # Different geoids and a surface vs a bare-earth model: close, not equal
    assert abs(float((cop30 - usgs10m).median())) < 5
