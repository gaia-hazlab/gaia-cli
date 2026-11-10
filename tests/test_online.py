"""Network smoke tests: stage each dataset for a small Bolt Creek AOI through the
stage-file runner. Run with ``pixi run test-online``; skipped by default."""

import json
import os

import geopandas as gpd
import numpy as np
import pandas as pd
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


def check_nrcs(ds):
    # Stevens Pass SNOTEL: hourly, so 72 reports over the padded three days
    assert "791:WA:SNTL" in ds.station.values.tolist()
    assert int(ds.row_size.sum()) == ds.sizes["obs"]
    assert ds.WTEQ.attrs["units"] == "mm" and float(ds.WTEQ.min()) >= 0
    assert -30 < float(ds.TOBS.min()) < float(ds.TOBS.max()) < 20  # degC
    assert stations.to_grid(ds).sizes["time"] >= 72


def check_prism(ds):
    assert ds.sizes["time"] == 3
    assert float(ds.ppt.min()) >= 0
    assert ds.ppt.attrs["gaia:common-name"] == "precipitation"


def check_prism_temperature(ds):
    assert ds.sizes["time"] == 3
    assert bool((ds.tmin <= ds.tmax).all())
    assert -40 < float(ds.tmin.min()) and float(ds.tmax.max()) < 30  # December, degC


def check_hrrr(ds):
    assert ds.sizes["time"] == 72
    assert float(ds.precipitation_surface.min()) >= 0


def check_era5(ds):
    # Hourly, hour-ending; three December days in the Cascades
    assert ds.sizes["time"] == 72
    assert float(ds.tp.min()) >= 0 and ds.tp.attrs["units"] == "mm"
    assert bool((ds.sf <= ds.tp + 1e-3).all())  # snowfall is part of the total
    assert -30 < float(ds.t2m.min()) < float(ds.t2m.max()) < 20  # degC
    assert 0 < float(ds.swvl1.min()) and float(ds.swvl2.max()) < 1
    assert 0 < float(ds.orography.min()) < float(ds.orography.max()) < 2500  # m
    assert ds.attrs["gaia:icechunk-snapshot"]


def check_stage4(ds):
    # Read from GDEX into a fresh cache: the 24 h totals ending 12 UTC on each day
    assert ds.sizes["time"] == 3
    assert float(ds.apcp.min()) >= 0
    assert ds.attrs["stage4:domain"] == "conus"


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
        ("nrcs", check_nrcs, {"variables": ["WTEQ", "TOBS"]}),
        ("prism", check_prism, {}),
        pytest.param(
            "prism",
            check_prism_temperature,
            {"variables": ["ppt", "tmin", "tmax"]},
            id="prism-temperature",
        ),
        ("hrrr", check_hrrr, {}),
        ("stage4/daily", check_stage4, {}),
        pytest.param(
            "era5",
            check_era5,
            {},
            marks=pytest.mark.skipif(
                not os.environ.get("ARRAYLAKE_TOKEN"), reason="needs ARRAYLAKE_TOKEN"
            ),
        ),
        ("nlcd", check_nlcd, {}),
        ("dem/cop30", check_dem, {}),
        ("dem/usgs10m", check_dem, {}),
        ("soil/solus100", check_soil, {}),
        ("burn-severity", check_burn_severity, {"fire": "bolt creek"}),
        ("soil/statsgo", check_statsgo, {}),
    ],
)
def test_stage(tmp_path, monkeypatch, dataset, check, options):
    monkeypatch.setenv("STAGE4_CACHE", str(tmp_path / "stage4"))
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


# Peters Dome, AK: the 2023-09-14 landslide lies between two passes of descending
# relative orbit 131 over the source area (only Sentinel-1A flew in 2023)
PETERS_S1 = """
name: peters-dome-s1
aoi: {latitude: 63.1455, longitude: -151.1904, radius_km: 5}
window: 2023-09-01/2023-09-30
datasets:
  sentinel1-rtc: {}
"""


def test_sentinel1_rtc_passes_around_peters_dome(tmp_path):
    config = tmp_path / "stage.yaml"
    config.write_text(PETERS_S1)
    runner.stage(config=str(config), output_path=str(tmp_path / "out"))

    tree = xr.open_datatree(tmp_path / "out" / "peters-dome-s1.zarr", engine="zarr")
    ds = tree["rasters/sentinel1-rtc"].to_dataset()
    assert [str(t)[:10] for t in ds.time.values] == ["2023-09-09", "2023-09-21"]
    assert ds.relative_orbit.values.tolist() == [131, 131]
    assert set(ds.orbit_state.values.tolist()) == {"descending"}
    assert ds.attrs["license"] == "CC-BY-4.0"
    # Linear gamma0: VH below VV, both mostly between -30 and +5 dB
    vv, vh = (float(ds[p].median()) for p in ("vv", "vh"))
    assert 10**-3 < vh < vv < 10**0.5
    assert ds.vv.attrs["gaia:regrid-method"] == "average"


def test_sentinel1_rtc_pre_and_post_passes_on_each_track(tmp_path):
    # The regional AOI sees tracks 102 and 131; the last pass before and the first
    # after the event on each. Loaded lazily, so nothing is read
    from gaia_cli.datasets import sentinel1_rtc
    from gaia_cli.window import parse_window

    windows = parse_window(
        {
            "pre": {"before": "30d", "select": "last"},
            "post": {"after": "30d", "select": "first"},
        },
        pd.Timestamp("2023-09-14T02:06:49"),
    )
    aoi = gpd.GeoDataFrame(
        geometry=[box(-152.2, 62.75, -150.2, 63.55)], crs="EPSG:4326"
    )
    ds = sentinel1_rtc.load(
        aoi,
        "2023-08-15T02:06:49",
        "2023-10-14T02:06:48",
        match=sentinel1_rtc.Match.relative_orbit,
        windows=windows,
    )
    assert [str(t)[:10] for t in ds.time.values] == [
        "2023-09-07",
        "2023-09-09",
        "2023-09-19",
        "2023-09-21",
    ]
    assert ds.window.values.tolist() == ["pre", "pre", "post", "post"]
    assert ds.relative_orbit.values.tolist() == [102, 131, 102, 131]


@pytest.mark.skipif(
    not os.environ.get("EARTHDATA_TOKEN"), reason="needs EARTHDATA_TOKEN"
)
def test_opera_rtc_passes_around_peters_dome(tmp_path):
    # The same passes as sentinel1-rtc, at 30 m, read lazily while the store is
    # written; the token must not reach the store
    config = tmp_path / "stage.yaml"
    config.write_text(PETERS_S1.replace("sentinel1-rtc", "opera-rtc"))
    runner.stage(config=str(config), output_path=str(tmp_path / "out"))

    tree = xr.open_datatree(tmp_path / "out" / "peters-dome-s1.zarr", engine="zarr")
    ds = tree["rasters/opera-rtc"].to_dataset()
    assert [str(t)[:10] for t in ds.time.values] == ["2023-09-09", "2023-09-21"]
    assert ds.relative_orbit.values.tolist() == [131, 131]
    assert ds.rio.resolution() == (30.0, -30.0)
    vv, vh = (float(ds[p].median()) for p in ("vv", "vh"))
    assert 10**-3 < vh < vv < 10**0.5
    assert {0, 2} <= set(np.unique(ds["mask"].values).tolist())  # valid, layover
    token = os.environ["EARTHDATA_TOKEN"]
    assert not any(
        token in p.read_text(errors="ignore")
        for p in (tmp_path / "out").rglob("*")
        if p.is_file()
    )


@pytest.mark.skipif(
    not os.environ.get("EARTHDATA_TOKEN"), reason="needs EARTHDATA_TOKEN"
)
def test_opera_rtc_static_layers_of_each_pass():
    # The pass on track 131 gets that track's static layers, read onto its grid
    from gaia_cli.aoi import load_aoi
    from gaia_cli.datasets import opera_rtc

    aoi = load_aoi({"latitude": 63.1455, "longitude": -151.1904, "radius_km": 2})
    static = ("mask", "local_incidence_angle", "incidence_angle")
    ds = opera_rtc.load(
        aoi, "2023-09-09T00:00:00", "2023-09-09T23:59:59", static=static
    ).load()
    assert ds.relative_orbit.values.tolist() == [131]
    # Steep terrain: local angles spread far around the ellipsoid's ~40 degrees
    assert 30 < float(ds.incidence_angle.median()) < 45
    lia = ds.local_incidence_angle
    assert float(lia.min()) < 20 and float(lia.max()) > 60
    # Layover on slopes facing the radar (small local angles), shadow on slopes
    # facing away (large ones)
    median = {code: float(lia.where(ds["mask"] == code).median()) for code in (0, 1, 2)}
    assert median[2] < median[0] < median[1]


@pytest.mark.skipif(
    not os.environ.get("EARTHDATA_TOKEN"), reason="needs EARTHDATA_TOKEN"
)
def test_nisar_gcov_provisional_pass_at_paradise(tmp_path):
    # One descending pass over Paradise, Mt Rainier, read remotely from the HDF5;
    # the token must not reach the store
    config = tmp_path / "stage.yaml"
    config.write_text(
        """
name: paradise-nisar
aoi: {latitude: 46.7827, longitude: -121.7477, radius_km: 1}
window: 2026-09-06/2026-09-06
datasets:
  nisar-gcov: {maturity: provisional, relative_orbit: 172}
"""
    )
    runner.stage(config=str(config), output_path=str(tmp_path / "out"))

    tree = xr.open_datatree(tmp_path / "out" / "paradise-nisar.zarr", engine="zarr")
    ds = tree["rasters/nisar-gcov"].to_dataset()
    assert [str(t)[:10] for t in ds.time.values] == ["2026-09-06"]
    assert ds.frames.values.tolist() == ["65"]
    assert ds.rio.resolution() == (10.0, -10.0)
    assert ds.attrs["version"] == "NISAR_L2_GCOV_PROVISIONAL_V1"
    # Linear gamma0 of forest and meadow at L-band: HV below HH
    hh, hv = (float(ds[t].median()) for t in ("HHHH", "HVHV"))
    assert 10**-2.5 < hv < hh < 10**0.5
    assert ds.HVHV.attrs["gaia:common-name"] == "gamma0_hv"
    assert 1 in np.unique(ds["mask"].values)
    token = os.environ["EARTHDATA_TOKEN"]
    assert not any(
        token in p.read_text(errors="ignore")
        for p in (tmp_path / "out").rglob("*")
        if p.is_file()
    )


def test_sentinel2_after_the_event_is_snow_covered():
    # The paper notes snowfall right after the event; the clear 2023-09-24 pass
    # over the source is mostly snow: bright in the visible, dark at 1.6 um
    from gaia_cli.aoi import load_aoi
    from gaia_cli.datasets import sentinel2

    aoi = load_aoi({"latitude": 63.1455, "longitude": -151.1904, "radius_km": 5})
    ds = sentinel2.load(aoi, "2023-09-24T00:00:00", "2023-09-24T23:59:59")
    assert ds.sizes["time"] == 1
    assert ds.mgrs_tiles.values[0].split() == [
        "05VNK",
        "05VNL",
    ]  # the AOI spans two tiles
    scene = ds.isel(time=0).load()
    inside = scene.scl.values > 0
    assert (scene.scl.values[inside] == 11).mean() > 0.5  # snow_or_ice
    assert float(scene.blue.median()) > 0.5 > float(scene.swir16.median())


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
    node = tree["waveforms/seismic"]
    assert set(node.children) == {"sites", "channels", "50hz"}
    ds = node["50hz"].to_dataset()
    assert ds.channel.values.tolist() == [f"AK.BPAW..BH{c}" for c in "ENZ"]
    assert ds.sizes["time"] == 120 * 50
    assert not ds.counts.isnull().any()
    # BPAW samples 8.4 ms after each step
    assert (ds.time_offset.values == np.timedelta64(8_400_000, "ns")).all()
    channels = node["channels"].to_dataset()
    assert channels.sensitivity.values.tolist() == [4.8876e8] * 3
    assert channels.sensitivity_units.values.tolist() == ["m/s"] * 3
    sites = node["sites"].to_dataset().xvec.decode_cf()
    assert sites.site.values.tolist() == ["AK.BPAW"]
    assert sites.distance_km.values == pytest.approx([106.8], abs=0.1)
    assert node.attrs["networks"] == ["AK"]


# The AV.KENI infrasound array, 279 km from Peters Dome: its six elements are
# location codes 01 to 06, each at its own position
KENI = """
name: keni
time: 2023-09-14T02:06:49Z
origin: {latitude: 63.1455, longitude: -151.1904}
aoi: {latitude: 60.6409, longitude: -151.0713, radius_km: 1}
datasets:
  infrasound: {window: {before: 49s, after: 71s}}
"""


def test_infrasound_keni_array(tmp_path):
    config = tmp_path / "stage.yaml"
    config.write_text(KENI)
    out = tmp_path / "out"
    runner.stage(config=str(config), output_path=str(out))

    node = io.open_datatree(str(out / "keni.zarr"))["waveforms/infrasound"]
    assert set(node.children) == {"sites", "channels", "100hz"}
    ds = node["100hz"].to_dataset()
    assert ds.channel.values.tolist() == [f"AV.KENI.0{i}.HDF" for i in range(1, 7)]
    assert ds.sizes["time"] == 120 * 100
    assert not ds.counts.isnull().any()
    channels = node["channels"].to_dataset()
    assert channels.sensitivity.values == pytest.approx([11541.508] * 6)
    assert channels.sensitivity_units.values.tolist() == ["Pa"] * 6
    sites = node["sites"].to_dataset().xvec.decode_cf()
    assert sites.distance_km.values == pytest.approx([279.2] * 6, abs=0.2)
    assert waveforms.to_physical(node, "100hz").name == "pressure"
