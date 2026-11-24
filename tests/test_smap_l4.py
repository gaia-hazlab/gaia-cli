import base64
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
import zarr
from shapely.geometry import box

from gaia_cli import dmrpp, earthdata, fetch, regrid
from gaia_cli.cli import stage_app
from gaia_cli.datasets import smap, smap_l4

# Cells (1, 1) to (2, 2) of the grid, near its upper-left corner, so the
# stand-in files need only its first few rows and columns (as test_smap's)
AOI = gpd.GeoDataFrame(
    geometry=[
        box(
            smap.X0 + 1.2 * smap.CELL,
            smap.Y0 - 2.8 * smap.CELL,
            smap.X0 + 2.8 * smap.CELL,
            smap.Y0 - 1.2 * smap.CELL,
        )
    ],
    crs=smap.CRS,
).to_crs(4326)
CELLS = 6  # rows and columns in a stand-in file
COLLECTION = {
    "title": "SMAP L4 Global 3-hourly 9 km EASE-Grid Surface and Root Zone Soil "
    "Moisture Geophysical Data V008",
    "license": "CC0-1.0",
    "links": [{"rel": "license", "href": "https://www.earthdata.nasa.gov/"}],
}
FILL = {"_FillValue": base64.b64encode(np.float64(-9999).tobytes()).decode()}
# name -> attributes of the stand-ins' arrays, as VirtualiZarr gives them
ARRAYS = {
    "sm_surface": {
        **FILL,
        "units": "m3 m-3",
        "long_name": "Top layer soil moisture (0-5 cm)",
    },
    "sm_rootzone": {**FILL, "units": "m3 m-3"},
    "sm_profile": {**FILL, "units": "m3 m-3"},
    "snow_mass": {**FILL, "units": "kg m-2"},
    "sm_surface_wetness": {**FILL, "units": "dimensionless"},
}


def step_file(value: float) -> dict[str, zarr.Array]:
    """A stand-in for one 3-hour file: its arrays by path, *value* but -9999
    at cell (1, 1)."""
    out = {}
    for name, attrs in ARRAYS.items():
        values = np.full((CELLS, CELLS), value, np.float32)
        values[1, 1] = -9999
        out[f"/{smap_l4.GROUP}/{name}"] = zarr.create_array(
            zarr.storage.MemoryStore(), data=values, attributes=attrs
        )
    return out


def record(end: str, produced: str = "2025-07-03T00:00:00.000Z") -> dict:
    """A CMR UMM-G record of the file of the 3 hours ending at *end*."""
    stop = pd.Timestamp(end)
    center = stop - pd.Timedelta(minutes=90)
    name = f"SMAP_L4_SM_gph_{center:%Y%m%dT%H%M%S}_Vv8012_001"
    path = f"nsidc-cumulus-prod-protected/SMAP/SPL4SMGP/008/{center:%Y/%m/%d}/{name}"
    https = "https://data.nsidc.earthdatacloud.nasa.gov"
    return {
        "meta": {"concept-id": "G2-NSIDC_CPRD"},
        "umm": {
            "GranuleUR": f"{name}.h5_u9loH5HD",
            "TemporalExtent": {
                "RangeDateTime": {
                    "BeginningDateTime": f"{stop - smap_l4.STEP:%Y-%m-%dT%H:%M:%S}.000Z",
                    "EndingDateTime": f"{stop:%Y-%m-%dT%H:%M:%S}.000Z",
                }
            },
            "SpatialExtent": {
                "HorizontalSpatialDomain": {
                    "Geometry": {
                        "BoundingRectangles": [
                            {
                                "WestBoundingCoordinate": -180,
                                "SouthBoundingCoordinate": -85.04456,
                                "EastBoundingCoordinate": 180,
                                "NorthBoundingCoordinate": 85.04456,
                            }
                        ]
                    }
                }
            },
            # The QA report is listed first
            "RelatedUrls": [
                {"Type": "GET DATA", "URL": f"{https}/{path}.qa"},
                {"Type": "GET DATA", "URL": f"{https}/{path}.h5"},
                {"Type": "GET DATA VIA DIRECT ACCESS", "URL": f"s3://{path}.h5"},
            ],
            "DataGranule": {"ProductionDateTime": produced},
        },
    }


def day_records(day: str) -> list[dict]:
    """Records of the 3 hours ending at each of *day*'s 8 time steps, from
    00:00 (the previous day's last file) to 21:00."""
    return [record(str(t)) for t in pd.date_range(day, periods=8, freq="3h")]


@pytest.fixture
def nsidc(monkeypatch):
    """NSIDC stand-ins: the search returns *records* (the file of each step
    of 2025-07-01 and those touching it), and each file holds its end hour /
    100 (its center's hour + 2), or the value *values* gives its URL."""
    stand_in = {
        "records": [
            *day_records("2025-07-01"),
            record("2025-07-02T00:00:00"),
        ],
        "values": {},
        "searches": [],
        "opened": [],
        "signed": [],
    }

    def search(concept_id, aoi, start=None, end=None):
        stand_in["searches"].append((concept_id, start, end))
        return stand_in["records"]

    def signed_urls(hrefs):
        stand_in["signed"].append(list(hrefs))
        return {h: f"{h}?Signature=batch" for h in hrefs}

    def open_arrays(url, paths, arrays, signed=None):
        assert url.endswith(".h5") and arrays == smap_l4.ARRAYS
        stand_in["opened"].append(url)
        assert signed in (None, f"{url}?Signature=batch")
        hour = int(url.rsplit("T", 1)[1][:2]) + 2
        arrays = step_file(stand_in["values"].get(url, hour / 100))
        return {p: arrays[p] for p in paths}

    def get(url, **kwargs):
        assert url.startswith(smap_l4.STAC)
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(earthdata, "cmr_search", search)
    monkeypatch.setattr(earthdata, "signed_urls", signed_urls)
    monkeypatch.setattr(dmrpp, "open_arrays", open_arrays)
    monkeypatch.setattr(fetch, "get", get)
    return stand_in


def test_load_reads_each_3_hours_of_the_window(nsidc):
    ds = smap_l4.load(AOI, "2025-07-01", "2025-07-01")
    # The steps' files, from the 3 hours before the first step
    assert nsidc["searches"] == [
        (smap_l4.CONCEPT_ID, "2025-06-30T21:00:00", "2025-07-01T21:00:00")
    ]
    assert set(ds.data_vars) == {*smap_l4.VARIABLES, "time_bounds"}
    assert list(ds.time.values) == list(
        pd.date_range("2025-07-01", periods=8, freq="3h")
    )
    ds = ds.compute()
    sm = ds.sm_surface.sel(time="2025-07-01T03:00")
    assert sm.shape == (2, 2) and sm.dtype == np.float32
    # The file of 00:00-03:00 holds 3 / 100; -9999 is NaN
    assert np.isnan(sm.values[0, 0]) and float(sm[1, 1]) == pytest.approx(0.03)
    # Each file opened once, the .h5 rather than its QA report; the file
    # ending 2025-07-02 00:00 only touches the window
    assert len(nsidc["opened"]) == len(set(nsidc["opened"])) == 8


def test_load_is_on_smaps_ease_grid(nsidc):
    ds = smap_l4.load(AOI, "2025-07-01", "2025-07-01")
    assert ds.rio.crs.to_epsg() == 6933
    assert ds.x.values == pytest.approx(smap.X0 + smap.CELL * np.array([1.5, 2.5]))
    assert ds.y.values == pytest.approx(smap.Y0 - smap.CELL * np.array([1.5, 2.5]))


def test_each_value_is_the_mean_of_the_3_hours_before_it(nsidc):
    ds = smap_l4.load(AOI, "2025-07-01", "2025-07-01")
    assert ds.time.attrs["bounds"] == "time_bounds"
    assert list(pd.to_datetime(ds.time_bounds.values[1])) == [
        pd.Timestamp("2025-07-01T00:00"),
        pd.Timestamp("2025-07-01T03:00"),
    ]
    assert ds.sm_surface.attrs["cell_methods"] == "time: mean"


def test_load_attributes(nsidc):
    variables = ("sm_surface", "sm_rootzone", "snow_mass", "sm_surface_wetness")
    ds = smap_l4.load(AOI, "2025-07-01", "2025-07-01", variables)
    sm = ds.sm_surface
    assert sm.attrs["units"] == "m3 m-3"
    assert sm.attrs["standard_name"] == "volume_fraction_of_condensed_water_in_soil"
    assert sm.attrs["long_name"] == "Top layer soil moisture (0-5 cm)"
    # The root zone isn't the surface layer: no shared standard name
    assert "standard_name" not in ds.sm_rootzone.attrs
    assert ds.snow_mass.attrs["standard_name"] == "surface_snow_amount"
    assert ds.sm_surface_wetness.attrs["units"] == "1"
    assert "standard_name" not in ds.sm_surface_wetness.attrs
    assert regrid.METHOD_ATTR not in sm.attrs  # the runner sets REGRID_METHOD
    assert ds.attrs["version"] == "SPL4SMGP_008"
    links = [x["href"] for x in ds.attrs["links"] if "/items/" in x["href"]]
    assert len(links) == 8
    assert links[0] == (
        f"{smap_l4.STAC}/items/SMAP_L4_SM_gph_20250630T223000_Vv8012_001.h5_u9loH5HD"
    )


def test_missing_steps_are_filled_with_a_warning(nsidc, caplog):
    nsidc["records"] = [
        r for r in nsidc["records"] if "20250701T043000" not in r["umm"]["GranuleUR"]
    ]
    ds = smap_l4.load(AOI, "2025-07-01", "2025-07-01").compute()
    assert (
        "1 of 8 time step(s) are missing from the archive (first 2025-07-01T06:00:00)"
        in caplog.text
    )
    assert ds.sizes["time"] == 8
    assert ds.sm_surface.sel(time="2025-07-01T06:00").isnull().all()


def test_no_files_is_an_error(nsidc):
    nsidc["records"] = [record("2025-07-02T00:00:00")]  # only touches the window
    with pytest.raises(ValueError, match="No SMAP SPL4SMGP files"):
        smap_l4.load(AOI, "2025-07-01", "2025-07-01")


def test_load_reads_the_newest_processing_of_a_step(nsidc, caplog):
    caplog.set_level(logging.INFO)
    newer = record("2025-07-01T03:00:00", produced="2025-10-01T00:00:00Z")
    newer["umm"]["RelatedUrls"][1]["URL"] = newer["umm"]["RelatedUrls"][1][
        "URL"
    ].replace(".h5", "_002.h5")
    nsidc["records"].append(newer)
    nsidc["values"][newer["umm"]["RelatedUrls"][1]["URL"]] = 0.4
    ds = smap_l4.load(AOI, "2025-07-01", "2025-07-01").compute()
    assert "skipped 1 older processing(s)" in caplog.text
    assert float(ds.sm_surface.sel(time="2025-07-01T03:00")[1, 1]) == pytest.approx(0.4)


def test_steps_of_the_runners_window():
    # The runner's inclusive bounds of the window 2025-07-01/2025-07-01,
    # padded by a day
    steps = smap_l4.steps("2025-06-30T00:00:00", "2025-07-02T23:59:59.999999")
    assert steps[0] == pd.Timestamp("2025-06-30T00:00")
    assert steps[-1] == pd.Timestamp("2025-07-02T21:00")
    assert len(steps) == 24
    # Within a day
    assert list(smap_l4.steps("2025-07-01T02:00", "2025-07-01T07:00")) == [
        pd.Timestamp("2025-07-01T03:00"),
        pd.Timestamp("2025-07-01T06:00"),
    ]


def test_cli_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["smap-l4", "-i", str(aoi), "-s", "2025-07-01", "-e", "2025-07-01"]
        + ["-o", str(tmp_path / "o"), "--variables", "sm_surface", "snow_mass"]
        + ["--dry-run"],
        result_action="return_value",
    )
    assert "variables: [sm_surface, snow_mass]" in capsys.readouterr().out


def test_an_unknown_variable_is_refused_before_reading(nsidc, tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    with pytest.raises(SystemExit):
        stage_app(
            ["smap-l4", "-i", str(aoi), "-s", "2025-07-01", "-e", "2025-07-01"]
            + ["-o", str(tmp_path / "o"), "--variables", "soil_moisture"],
            result_action="return_value",
        )
    out = capsys.readouterr()
    assert "soil_moisture" in out.out + out.err
    assert nsidc["searches"] == []


def test_the_runner_writes_it(nsidc, tmp_path):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["smap-l4", "-i", str(aoi), "-s", "2025-07-01", "-e", "2025-07-01"]
        + ["-o", str(tmp_path / "o")],
        result_action="return_value",
    )
    (store,) = (tmp_path / "o").glob("*.zarr")
    ds = xr.open_datatree(store, engine="zarr")["rasters/smap-l4"].to_dataset()
    assert ds.sm_surface.attrs[regrid.METHOD_ATTR] == "average"
    # The window and a day beyond each end, of which the stand-in has 2025-07-01
    assert ds.sizes["time"] == 24
    assert float(ds.sm_rootzone.sel(time="2025-07-01T21:00").max()) == pytest.approx(
        0.21
    )
