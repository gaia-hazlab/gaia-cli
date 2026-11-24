import asyncio
import base64

import numpy as np
import pytest
import zarr

from gaia_cli import dmrpp, earthdata, fetch


def response(**attrs):
    return type("R", (), {"raise_for_status": lambda self: None, **attrs})()


def test_fill_values_from_virtualizarr():
    fill = base64.b64encode(np.float64(-9999).tobytes()).decode()
    assert dmrpp.fill_value({"_FillValue": fill}, np.dtype(np.float32)) == -9999
    assert dmrpp.fill_value({"_FillValue": 65534}, np.dtype(np.uint16)) == 65534


def test_granule_fetches_its_dmrpp_once_with_the_token(monkeypatch):
    calls = []

    def get(url, headers):
        calls.append((url, headers))
        return response(status_code=200, content=b"<xml/>")

    monkeypatch.setenv(earthdata.TOKEN, "t")
    monkeypatch.setattr(fetch, "get", get)
    store = dmrpp.Granule("https://nsidc/f.h5")
    assert store.head("f.h5.dmrpp") == {"size": 6}
    assert bytes(store.get("f.h5.dmrpp").buffer()) == b"<xml/>"
    assert store.get_ranges("f.h5.dmrpp", starts=[0, 4], lengths=[2, 2]) == [
        b"<x",
        b"/>",
    ]
    assert calls == [("https://nsidc/f.h5.dmrpp", {"Authorization": "Bearer t"})]


def test_a_file_without_a_dmrpp_is_a_clear_error(monkeypatch):
    monkeypatch.setenv(earthdata.TOKEN, "t")
    monkeypatch.setattr(fetch, "get", lambda url, headers: response(status_code=404))
    store = dmrpp.Granule("https://nsidc/SMAP_E_20250901.h5")
    with pytest.raises(RuntimeError, match=r"no DMR\+\+ index for SMAP_E_20250901\.h5"):
        store.head("SMAP_E_20250901.h5.dmrpp")


def test_granule_signs_its_file_once_and_fetches_each_block_once(monkeypatch):
    data = bytes(range(256)) * 64  # 16 kB
    signs, calls = [], []

    def signed_url(href):
        signs.append(href)
        return "https://cdn/f.h5?Signature=x"

    def get(url, headers):
        start, end = map(int, headers["Range"].removeprefix("bytes=").split("-"))
        calls.append((url, start))
        return response(content=data[start : end + 1])

    monkeypatch.setattr(dmrpp, "BLOCK", 1000)
    monkeypatch.setattr(earthdata, "signed_url", signed_url)
    monkeypatch.setattr(fetch, "get", get)
    store = dmrpp.Granule("https://nsidc/f.h5")
    assert signs == []  # not until a chunk is read

    async def reads():
        return await asyncio.gather(
            store.get_range_async("ignored", start=10, end=20),
            store.get_range_async("ignored", start=990, length=20),  # 2 blocks
            store.get_range_async("ignored", start=1500, end=1600),
        )

    a, b, c = asyncio.run(reads())
    assert (a, b, c) == (data[10:20], data[990:1010], data[1500:1600])
    assert signs == ["https://nsidc/f.h5"]
    # Blocks 0 and 1, each once, from the signed URL with its query
    assert sorted(calls) == [
        ("https://cdn/f.h5?Signature=x", 0),
        ("https://cdn/f.h5?Signature=x", 1000),
    ]


def test_granule_given_a_signed_url_reads_it_without_signing(monkeypatch):
    calls = []

    def signed_url(href):
        raise AssertionError("signed again")

    def get(url, headers):
        calls.append(url)
        return response(content=b"x" * 1000)

    monkeypatch.setattr(dmrpp, "BLOCK", 1000)
    monkeypatch.setattr(earthdata, "signed_url", signed_url)
    monkeypatch.setattr(fetch, "get", get)
    store = dmrpp.Granule("https://nsidc/f.h5", "https://cdn/f.h5?Signature=batch")
    assert asyncio.run(store.get_range_async("ignored", start=0, end=10)) == b"x" * 10
    assert calls == ["https://cdn/f.h5?Signature=batch"]


def test_open_arrays_skips_the_other_arrays_of_each_group(monkeypatch):
    import virtualizarr.parsers

    seen = []

    class Parser:
        def __init__(self, group, skip_variables):
            seen.append((group, skip_variables))

        def __call__(self, url, registry):
            store = zarr.storage.MemoryStore()
            zarr.create_array(store, name="sm_surface", shape=(2, 2), dtype="f4")
            return store

    monkeypatch.setattr(virtualizarr.parsers, "DMRPPParser", Parser)
    out = dmrpp.open_arrays(
        "https://nsidc/L4.h5",
        ["/Geophysical_Data/sm_surface"],
        {"Geophysical_Data": ["sm_surface", "snow_mass", "sm_rootzone"]},
    )
    assert seen == [("Geophysical_Data", ["snow_mass", "sm_rootzone"])]
    assert out["/Geophysical_Data/sm_surface"].shape == (2, 2)


def test_signed_url_is_the_daacs_redirect(monkeypatch):
    calls = []

    def get(url, headers, allow_redirects):
        calls.append((url, headers, allow_redirects))
        return response(status_code=303, headers={"location": "https://cdn/f?S=x"})

    monkeypatch.setenv(earthdata.TOKEN, "t")
    monkeypatch.setattr(fetch, "get", get)
    assert earthdata.signed_url("https://nsidc/f") == "https://cdn/f?S=x"
    assert calls == [("https://nsidc/f", {"Authorization": "Bearer t"}, False)]


def test_signed_url_needs_a_current_token(monkeypatch):
    def get(url, headers, allow_redirects):
        return response(
            status_code=302, headers={"location": "https://urs.earthdata.nasa.gov/x"}
        )

    monkeypatch.setenv(earthdata.TOKEN, "t")
    monkeypatch.setattr(fetch, "get", get)
    with pytest.raises(RuntimeError, match="could not sign"):
        earthdata.signed_url("https://nsidc/f")
