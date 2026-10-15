from gaia_cli import hrrr


def test_open_hrrr_pins_aws_endpoint(monkeypatch):
    captured = {}

    def fake_s3_storage(**kwargs):
        captured.update(kwargs)
        raise RuntimeError("stop")

    monkeypatch.setattr(hrrr.icechunk, "s3_storage", fake_s3_storage)
    try:
        hrrr.open_hrrr(region="us-west-2")
    except RuntimeError:
        pass
    assert captured["endpoint_url"] == "https://s3.us-west-2.amazonaws.com"
    assert captured["anonymous"] is True
