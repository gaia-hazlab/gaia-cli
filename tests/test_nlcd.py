import numpy as np
import xarray as xr

from gaia_cli.datasets import nlcd


def test_nlcd_classes_are_cf_flags(monkeypatch):
    da = xr.DataArray(np.full((1, 2, 2), 42, "uint8"), dims=("band", "y", "x"))
    monkeypatch.setattr(nlcd.rioxarray, "open_rasterio", lambda *a, **k: da)
    attrs = nlcd.open_nlcd(2021).attrs
    n = len(attrs["flag_values"])
    assert (
        n
        == len(attrs["flag_meanings"].split())
        == len(attrs["flag_colors"].split())
        == 16
    )
    assert (
        attrs["flag_meanings"].split()[attrs["flag_values"].index(42)]
        == "Evergreen_Forest"
    )
