"""Network smoke tests: stage each dataset for a small Bolt Creek AOI through the
stage-file runner. Run with ``pixi run test-online``; skipped by default."""

import json
import os

import geopandas as gpd
import numpy as np
import pytest
import xarray as xr
from shapely.geometry import box

from gaia_cli import align, io, runner, stations, waveforms

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


# Each window is read one day beyond each end (PAD)
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


def check_statsgo(ds):
    # STATSGO K factors are 0.02 to 0.62; -0.1 (no K factor) is masked
    kf = ds.kffact.values[~np.isnan(ds.kffact.values)]
    assert kf.size and 0 < kf.min() and kf.max() < 0.7
    assert ds.kffact.attrs["gaia:regrid-method"] == "nearest"


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
        ("burn-severity", check_burn_severity, {"fire": "bolt creek"}),
        ("soil/statsgo", check_statsgo, {}),
    ],
)
def test_stage(tmp_path, dataset, check, options):
    config = tmp_path / "stage.yaml"
    config.write_text(STAGE.format(dataset=dataset, options=json.dumps(options)))
    out = tmp_path / "out"
    runner.stage(config=str(config), output_path=str(out))

    tree = xr.open_datatree(out / "bolt-creek-smoke.zarr", engine="zarr")
    (path,) = [p for p in tree.groups if p.count("/") == 2]
    assert path.endswith(dataset)
    check(tree[path].to_dataset())
    assert tree[path].attrs["license"]


# The landlab_data_prep reference event: Pioneer Fire, WA (WA4818312053120240608),
# the perimeter's bbox. Landlab puts every layer on a 10 m UTM grid.
PIONEER_BBOX = [-120.6807, 48.0978, -120.4009, 48.3895]
PIONEER = f"""
name: pioneer
aoi: {{bbox: {PIONEER_BBOX}}}
datasets:
  dem/cop30: {{}}
  dem/usgs10m: {{}}
"""


def test_pioneer_dems_on_the_landlab_grid(tmp_path):
    config = tmp_path / "stage.yaml"
    config.write_text(PIONEER)
    out = tmp_path / "out"
    runner.stage(config=str(config), output_path=str(out))

    tree = xr.open_datatree(out / "pioneer.zarr", engine="zarr")
    aoi = gpd.GeoDataFrame(geometry=[box(*PIONEER_BBOX)], crs="EPSG:4326")
    # The DEM's default method (bilinear), as gaia align would use it
    aligned = align.align(tree, crs="utm", res=10.0, aoi=aoi)
    cop30, usgs10m = (aligned[f"dem/{p}"].elevation for p in ("cop30", "usgs10m"))
    assert cop30.attrs["regrid_method"] == usgs10m.attrs["regrid_method"] == "bilinear"
    assert cop30.rio.crs.to_epsg() == usgs10m.rio.crs.to_epsg() == 32610
    # The contract's grid: AOI bounds snapped outward to 10 m
    assert cop30.shape == usgs10m.shape == (3309, 2182)
    assert cop30.rio.bounds() == (671_690, 5_329_770, 693_510, 5_362_860)
    # Lake Chelan's shore (~335 m) to the Chelan Mountains (~2680 m)
    for da in (cop30, usgs10m):
        assert 320 < float(da.min()) < 350
        assert 2600 < float(da.max()) < 2750
    # Different geoids and a surface vs a bare-earth model: close, not equal
    assert abs(float((cop30 - usgs10m).median())) < 5


# AK.BPAW, 107 km from Peters Dome, two minutes around the 2023-09-14 avalanche
BPAW = """
name: bpaw
time: 2023-09-14T02:06:49Z
origin: {latitude: 63.1455, longitude: -151.1904}
aoi: {latitude: 64.0997, longitude: -150.9873, radius_km: 5}
datasets:
  seismic: {window: {before: 49s, after: 71s}}
"""


def test_seismic_bpaw(tmp_path):
    config = tmp_path / "stage.yaml"
    config.write_text(BPAW)
    out = tmp_path / "out"
    runner.stage(config=str(config), output_path=str(out))

    tree = io.open_datatree(str(out / "bpaw.zarr"))
    # BPAW has BH channels only, at 50 Hz
    assert set(tree["waveforms/seismic"].children) == {"50hz"}
    ds = waveforms.decoded(tree["waveforms/seismic/50hz"].to_dataset())
    assert ds.channel.values.tolist() == [f"AK.BPAW..BH{c}" for c in "ENZ"]
    assert ds.sizes["time"] == 120 * 50
    assert not ds.counts.isnull().any()
    assert ds.sensitivity.values.tolist() == [4.8876e8] * 3
    assert ds.sensitivity_units.values.tolist() == ["m/s"] * 3
    # BPAW samples 8.4 ms after each step
    assert (ds.time_offset.values == np.timedelta64(8_400_000, "ns")).all()
    assert ds.distance_km.values == pytest.approx([106.8] * 3, abs=0.1)
    assert ds.attrs["networks"] == ["AK"]
