import numpy as np
import rioxarray  # noqa: F401  registers the .rio accessor
import xarray as xr
from affine import Affine

from gaia_cli import regrid


def _grid(x0, y0, res, nx, ny, crs, data=None):
    x = x0 + res * (np.arange(nx) + 0.5)
    y = y0 - res * (np.arange(ny) + 0.5)
    if data is None:
        data = np.zeros((ny, nx))
    da = xr.DataArray(data, coords={"y": y, "x": x}, dims=("y", "x"), name="v")
    return da.rio.write_crs(crs)


def test_conservative_equals_block_mean_on_nested_grids():
    rng = np.random.default_rng(0)
    src = _grid(500_000, 5_300_000, 1_000, 8, 6, "EPSG:32610", rng.random((6, 8)))
    tgt = _grid(500_000, 5_300_000, 2_000, 4, 3, "EPSG:32610")

    out = regrid.regrid(src, tgt, "conservative")

    expected = src.coarsen(x=2, y=2).mean()
    np.testing.assert_allclose(
        out.values, expected.values, rtol=1e-6
    )  # UTM → local LAEA area distortion
    assert out.rio.crs == tgt.rio.crs


def test_conservative_skips_nan_and_renormalizes():
    data = np.array([[1.0, np.nan], [3.0, 5.0]])
    src = _grid(500_000, 5_300_000, 1_000, 2, 2, "EPSG:32610", data)
    tgt = _grid(500_000, 5_300_000, 2_000, 1, 1, "EPSG:32610")
    # a single cell has no coordinate spacing; stage output carries a GeoTransform
    tgt = tgt.rio.write_transform(Affine(2_000, 0, 500_000, 0, -2_000, 5_300_000))

    out = regrid.regrid(src, tgt, "conservative")

    np.testing.assert_allclose(out.values, [[3.0]])


def test_conservative_across_crs_preserves_constant_and_extra_dims():
    src = _grid(-122.5, 48.5, 0.01, 60, 40, "EPSG:4326", np.full((40, 60), 7.0))
    src = src.expand_dims(time=3)
    # 3 km UTM grid well inside the source extent
    tgt = _grid(560_000, 5_365_000, 3_000, 5, 4, "EPSG:32610")

    out = regrid.regrid(src, tgt, "conservative")

    assert out.dims == ("time", "y", "x")
    np.testing.assert_allclose(out.values, 7.0, rtol=1e-9)


def test_cell_areas_geographic():
    da = _grid(-122.0, 48.0, 0.01, 2, 2, "EPSG:4326")
    # ~0.01° × 0.01° near 48°N ≈ 1.11 km × 0.74 km
    np.testing.assert_allclose(regrid.cell_areas(da), 1.11e3 * 0.744e3, rtol=0.01)


def test_regrid_rio_methods_match_template_coords():
    src = _grid(500_000, 5_300_000, 1_000, 8, 6, "EPSG:32610", np.ones((6, 8)))
    tgt = _grid(500_000, 5_300_000, 2_000, 4, 3, "EPSG:32610")
    for method in ["average", "bilinear", "nearest"]:
        out = regrid.regrid(src, tgt, method)
        np.testing.assert_array_equal(out.x, tgt.x)
        np.testing.assert_allclose(out.values, 1.0)
