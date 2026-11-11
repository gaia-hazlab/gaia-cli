import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from osgeo import gdal
from shapely.geometry import box

from gaia_cli import earthdata, fetch, io, regrid
from gaia_cli.cli import stage_app
from gaia_cli.datasets import smap
from gaia_cli.datasets.sentinel1_rtc import OrbitState

# Cells (1, 1) to (2, 2) of the grid, near its upper-left corner, so the
# stand-in files need only its first few rows and columns
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
    "title": "SMAP Enhanced L3 Radiometer Global and Polar Grid Daily 9 km",
    "license": "CC0-1.0",
    "links": [{"rel": "license", "href": "https://www.earthdata.nasa.gov/"}],
}
FLAGS = {"flag_masks": "1s, 2s, 4s, 8s", "flag_meanings": "a b c d"}
# AM name -> (type, unit, attributes) of the stand-ins' arrays
ARRAYS = {
    "soil_moisture": (gdal.GDT_Float32, "cm**3/cm**3", {"long_name": "soil moisture"}),
    "soil_moisture_error": (gdal.GDT_Float32, "cm**3/cm**3", {}),
    "tb_v_corrected": (gdal.GDT_Float32, "Kelvin", {"valid_min": 0, "valid_max": 330}),
    "tb_h_corrected": (gdal.GDT_Float32, "Kelvin", {}),
    "retrieval_qual_flag": (gdal.GDT_UInt16, "", FLAGS),
    "tb_qual_flag_v": (gdal.GDT_UInt16, "", FLAGS),
    "tb_qual_flag_h": (gdal.GDT_UInt16, "", FLAGS),
    "surface_flag": (gdal.GDT_UInt16, "", FLAGS),
}


def day_file(value: float) -> gdal.Dataset:
    """A stand-in for one day's file: every array of both passes, floats
    *value* (PM: *value* + 0.5) but -9999 at cell (1, 1), flags 1 (PM: 2)."""
    ds = gdal.GetDriverByName("MEM").CreateMultiDimensional("smap")
    root = ds.GetRootGroup()
    for (group, suffix), offset in zip(smap.GROUPS.values(), (0, 1), strict=True):
        g = root.CreateGroup(group)
        dims = [g.CreateDimension(f"dim{i}", None, None, CELLS) for i in (0, 1)]
        for name, (kind, unit, attrs) in ARRAYS.items():
            array = g.CreateMDArray(
                f"{name}{suffix}", dims, gdal.ExtendedDataType.Create(kind)
            )
            if kind == gdal.GDT_Float32:
                values = np.full((CELLS, CELLS), value + offset / 2, np.float32)
                values[1, 1] = -9999
                array.SetNoDataValueDouble(-9999)
            else:
                values = np.full((CELLS, CELLS), 1 + offset, np.uint16)
                array.SetNoDataValueDouble(65534)
            array.Write(values)
            if unit:
                array.SetUnit(unit)
            for key, v in attrs.items():
                kind = (
                    gdal.ExtendedDataType.CreateString()
                    if isinstance(v, str)
                    else gdal.ExtendedDataType.Create(gdal.GDT_Float64)
                )
                array.CreateAttribute(key, [], kind).Write(v)
    return ds


def record(day: str, produced: str = "2025-09-02T11:33:38.529Z") -> dict:
    """A CMR UMM-G record of *day*'s file."""
    name = f"SMAP_L3_SM_P_E_{day.replace('-', '')}_R19240_001.h5"
    path = f"nsidc-cumulus-prod-protected/SMAP/SPL3SMP_E/006/{name}"
    return {
        "meta": {"concept-id": "G1-NSIDC_CPRD"},
        "umm": {
            "GranuleUR": name,
            "TemporalExtent": {
                "RangeDateTime": {
                    "BeginningDateTime": f"{day}T00:00:00.000Z",
                    "EndingDateTime": f"{day}T23:59:59.999Z",
                }
            },
            "SpatialExtent": {
                "HorizontalSpatialDomain": {
                    "Geometry": {
                        "BoundingRectangles": [
                            {
                                "WestBoundingCoordinate": -180,
                                "SouthBoundingCoordinate": -85.0445,
                                "EastBoundingCoordinate": 180,
                                "NorthBoundingCoordinate": 85.0445,
                            }
                        ]
                    }
                }
            },
            "RelatedUrls": [
                {"Type": "GET DATA VIA DIRECT ACCESS", "URL": f"s3://{path}"},
                {
                    "Type": "GET DATA",
                    "URL": f"https://data.nsidc.earthdatacloud.nasa.gov/{path}",
                },
                {
                    "Type": "EXTENDED METADATA",
                    "URL": f"https://data.nsidc.earthdatacloud.nasa.gov/{path}.iso.xml",
                },
            ],
            "DataGranule": {"ProductionDateTime": produced},
        },
    }


@pytest.fixture
def nsidc(monkeypatch):
    """NSIDC stand-ins: the search returns *records*, and each file holds the
    value *values* gives its URL. Returns both, and the searches made."""
    stand_in = {
        "records": [record("2025-09-01"), record("2025-09-02"), record("2025-09-03")],
        "values": {},
        "searches": [],
    }

    def search(concept_id, aoi, start=None, end=None):
        stand_in["searches"].append((concept_id, start, end))
        return stand_in["records"]

    def open_file(href):
        day = href.rsplit("_E_", 1)[1][:8]
        return day_file(stand_in["values"].get(href, int(day[-2:]) / 10))

    def get(url, **kwargs):
        assert url.startswith(smap.STAC)
        return type(
            "R",
            (),
            {"raise_for_status": lambda self: None, "json": lambda self: COLLECTION},
        )()

    monkeypatch.setattr(earthdata, "cmr_search", search)
    monkeypatch.setattr(earthdata, "sign", lambda items, assets: None)
    monkeypatch.setattr(smap, "open_file", open_file)
    monkeypatch.setattr(fetch, "get", get)
    return stand_in


def test_load_reads_both_passes_of_each_day(nsidc):
    ds = smap.load(AOI, "2025-09-01", "2025-09-03")
    assert nsidc["searches"] == [(smap.CONCEPT_ID, "2025-09-01", "2025-09-03")]
    assert set(ds.data_vars) == {
        *smap.VARIABLES,
        *(f"{v}_pm" for v in smap.VARIABLES),
        "time_bounds",
    }
    assert [str(t)[:10] for t in ds.time.values] == [
        "2025-09-01",
        "2025-09-02",
        "2025-09-03",
    ]
    ds = ds.compute()
    # The AOI's 2 x 2 cells, each day's value; -9999 is NaN
    sm = ds.soil_moisture.sel(time="2025-09-02")
    assert sm.shape == (2, 2)
    assert np.isnan(sm.values[0, 0]) and float(sm.values[1, 1]) == pytest.approx(0.2)
    assert float(ds.soil_moisture_pm.sel(time="2025-09-02")[1, 1]) == pytest.approx(0.7)
    assert ds.soil_moisture.dtype == np.float32
    flags = ds.surface_flag_pm.sel(time="2025-09-02")
    assert flags.dtype == np.uint16 and set(np.unique(flags.values)) == {2}


def test_load_is_on_the_ease_grid(nsidc):
    ds = smap.load(AOI, "2025-09-01", "2025-09-01")
    assert ds.rio.crs.to_epsg() == 6933
    assert ds.rio.resolution() == pytest.approx((smap.CELL, -smap.CELL))
    # Cell centers of rows and columns 1 and 2
    assert ds.x.values == pytest.approx(smap.X0 + smap.CELL * np.array([1.5, 2.5]))
    assert ds.y.values == pytest.approx(smap.Y0 - smap.CELL * np.array([1.5, 2.5]))


def test_load_attributes(nsidc):
    ds = smap.load(AOI, "2025-09-01", "2025-09-01")
    sm, tb, flag = ds.soil_moisture, ds.tb_v_corrected_pm, ds.retrieval_qual_flag
    assert sm.attrs["units"] == "m3 m-3"
    assert sm.attrs["standard_name"] == "volume_fraction_of_condensed_water_in_soil"
    assert sm.attrs["long_name"] == "soil moisture"
    assert (
        tb.attrs["units"] == "K"
        and tb.attrs["standard_name"] == "brightness_temperature"
    )
    assert (tb.attrs["valid_min"], tb.attrs["valid_max"]) == (0, 330)
    # Bit flags keep their fill value and masks, and regrid with nearest
    assert flag.encoding["_FillValue"] == 65534
    assert flag.attrs["flag_masks"] == [1, 2, 4, 8]
    assert flag.attrs["flag_meanings"] == "a b c d"
    assert flag.attrs[regrid.METHOD_ATTR] == "nearest"
    assert regrid.METHOD_ATTR not in sm.attrs  # the runner sets REGRID_METHOD
    assert ds.attrs["version"] == "SPL3SMP_E_006"
    links = [x["href"] for x in ds.attrs["links"] if "/items/" in x["href"]]
    assert links == [f"{smap.STAC}/items/SMAP_L3_SM_P_E_20250901_R19240_001.h5"]


def test_each_day_is_bounded_by_its_utc_day(nsidc):
    ds = smap.load(AOI, "2025-09-01", "2025-09-02")
    assert ds.time.attrs["bounds"] == "time_bounds"
    assert list(pd.to_datetime(ds.time_bounds.values[1])) == [
        pd.Timestamp("2025-09-02"),
        pd.Timestamp("2025-09-03"),
    ]
    assert ds.soil_moisture.attrs["cell_methods"] == "time: point"


def test_missing_days_are_filled_with_a_warning(nsidc, caplog):
    nsidc["records"] = [record("2025-09-01"), record("2025-09-03")]
    ds = smap.load(AOI, "2025-09-01", "2025-09-03").compute()
    assert (
        "1 of 3 day(s) are missing from the archive (first 2025-09-02)" in caplog.text
    )
    assert ds.sizes["time"] == 3
    assert ds.soil_moisture.sel(time="2025-09-02").isnull().all()
    assert ds.surface_flag.dtype == np.uint16
    assert (ds.surface_flag.sel(time="2025-09-02") == 65534).all()


def test_no_files_is_an_error(nsidc):
    nsidc["records"] = []
    with pytest.raises(ValueError, match="No SMAP SPL3SMP_E files"):
        smap.load(AOI, "2025-09-01", "2025-09-03")


def test_load_reads_the_newest_processing_of_a_day(nsidc, caplog):
    caplog.set_level(logging.INFO)
    newer = record("2025-09-01", produced="2025-10-01T00:00:00Z")
    newer["umm"]["RelatedUrls"][1]["URL"] += "?newer"
    nsidc["records"] = [record("2025-09-01"), newer]
    nsidc["values"][newer["umm"]["RelatedUrls"][1]["URL"]] = 0.4
    ds = smap.load(AOI, "2025-09-01", "2025-09-01").compute()
    assert "skipped 1 older processing(s)" in caplog.text
    assert float(ds.soil_moisture[0, 1, 1]) == pytest.approx(0.4)


@pytest.mark.parametrize(
    "orbit_state, names",
    [
        (OrbitState.descending, {"soil_moisture", "surface_flag"}),
        (OrbitState.ascending, {"soil_moisture_pm", "surface_flag_pm"}),
    ],
)
def test_orbit_state_picks_the_passes(nsidc, orbit_state, names):
    variables = ("soil_moisture", "surface_flag")
    ds = smap.load(AOI, "2025-09-01", "2025-09-01", variables, orbit_state)
    assert set(ds.data_vars) == {*names, "time_bounds"}


def test_large_datasets_are_chunked_by_time_steps(nsidc, monkeypatch):
    # Two of the AOI's 2 x 2 float32 time steps (16 bytes each) to a chunk
    monkeypatch.setattr(smap, "CHUNK_BYTES", 32)
    ds = smap.load(AOI, "2025-09-01", "2025-09-03", ("soil_moisture",))
    assert ds.soil_moisture.chunks == ((2, 1), (2,), (2,))


def test_chunks():
    assert smap.chunks(283, 623) == {"time": 14, "y": -1, "x": -1}  # CONUS
    # A global time step is over 10 MB: bands of rows
    assert smap.chunks(1624, 3856) == {"time": 1, "y": 679, "x": -1}


def test_window_is_the_cells_the_aoi_touches_and_one_more():
    paradise = gpd.GeoDataFrame(
        geometry=[box(-121.9, 46.7, -121.6, 46.9)], crs="EPSG:4326"
    )
    assert smap.window(paradise) == (slice(217, 222), slice(621, 627))
    # Within the grid at its edges
    assert smap.window(AOI) == (slice(0, 4), slice(0, 4))


def test_an_aoi_beyond_the_grid_is_an_error():
    pole = gpd.GeoDataFrame(geometry=[box(0, 86, 10, 88)], crs="EPSG:4326")
    with pytest.raises(ValueError, match="outside SMAP's grid"):
        smap.window(pole)


def test_cmr_search_pages_through_the_results(monkeypatch):
    monkeypatch.setattr(earthdata, "CMR_PAGE_SIZE", 2)
    pages = [[record("2025-09-01"), record("2025-09-02")], [record("2025-09-03")]]
    calls = []

    class Response:
        def __init__(self, items, after):
            self.items = items
            self.headers = {"CMR-Search-After": after} if after else {}

        def raise_for_status(self):
            pass

        def json(self):
            return {"items": self.items}

    def get(url, params, headers):
        calls.append((url, dict(params), headers))
        page = pages[len(calls) - 1]
        return Response(page, "token" if len(calls) == 1 else None)

    monkeypatch.setattr(fetch, "get", get)
    records = earthdata.cmr_search(
        smap.CONCEPT_ID, AOI, "2025-09-01", "2025-09-03T23:59:59.999999"
    )
    assert len(records) == 3
    assert [c[2] for c in calls] == [{}, {"CMR-Search-After": "token"}]
    url, params, _ = calls[0]
    assert url == earthdata.CMR
    assert params["concept_id"] == smap.CONCEPT_ID
    assert params["temporal"] == "2025-09-01T00:00:00Z,2025-09-03T23:59:59Z"
    west, south, east, north = map(float, params["bounding_box"].split(","))
    assert west < east and south < north


def test_cmr_item_is_the_granules_https_file():
    item = earthdata.cmr_item(record("2025-09-01"))
    assert item.id == "SMAP_L3_SM_P_E_20250901_R19240_001.h5"
    assert str(item.datetime) == "2025-09-01 00:00:00+00:00"
    assert item.properties["end_datetime"] == "2025-09-01T23:59:59.999Z"
    assert item.properties["processing:datetime"] == "2025-09-02T11:33:38.529Z"
    assert item.assets["data"].href == (
        "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected/"
        "SMAP/SPL3SMP_E/006/SMAP_L3_SM_P_E_20250901_R19240_001.h5"
    )
    assert item.bbox == [-180, -85.0445, 180, 85.0445]


def test_cli_flags(tmp_path, capsys):
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["smap", "-i", str(aoi), "-s", "2025-09-01", "-e", "2025-09-03"]
        + ["-o", str(tmp_path / "o"), "--orbit-state", "descending"]
        + ["--variables", "soil_moisture", "--dry-run"],
        result_action="return_value",
    )
    out = capsys.readouterr().out
    assert "orbit_state: descending" in out and "variables: [soil_moisture]" in out


def test_the_runner_reads_it_in_processes_and_writes_it(nsidc, tmp_path, monkeypatch):
    calls = []

    def compute(tree, env):
        calls.append(env)
        return tree

    monkeypatch.setattr(io, "compute_in_processes", compute)
    aoi = tmp_path / "aoi.geojson"
    AOI.to_file(aoi)
    stage_app(
        ["smap", "-i", str(aoi), "-s", "2025-09-01", "-e", "2025-09-03"]
        + ["-o", str(tmp_path / "o")],
        result_action="return_value",
    )
    assert calls == [smap.PROCESS_ENV]
    assert calls[0]["CPL_VSIL_CURL_CHUNK_SIZE"] == "16384"
    (store,) = (tmp_path / "o").glob("*.zarr")
    ds = xr.open_datatree(store, engine="zarr")["rasters/smap"].to_dataset()
    # Stored as uint16 with its fill value (read back as float with NaN)
    assert ds.surface_flag.encoding["dtype"] == np.uint16
    assert ds.surface_flag.encoding["_FillValue"] == 65534
    assert ds.soil_moisture.attrs[regrid.METHOD_ATTR] == "average"
    assert ds.surface_flag.attrs[regrid.METHOD_ATTR] == "nearest"
    assert float(ds.tb_v_corrected.max()) == pytest.approx(0.3)
