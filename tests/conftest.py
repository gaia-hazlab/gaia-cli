import numpy as np
import pandas as pd
import pymseed
import pytest
import requests

from gaia_cli import waveforms


@pytest.fixture
def no_aws_credentials(monkeypatch, tmp_path):
    """Hide every AWS credential source from boto3."""
    for var in ("AWS_PROFILE", "AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("AWS_CONFIG_FILE", str(tmp_path / "config"))
    monkeypatch.setenv("AWS_SHARED_CREDENTIALS_FILE", str(tmp_path / "credentials"))
    monkeypatch.setenv("AWS_EC2_METADATA_DISABLED", "true")


STATION_HEADER = (
    "#Network | Station | Location | Channel | Latitude | Longitude | Elevation | "
    "Depth | Azimuth | Dip | SensorDescription | Scale | ScaleFreq | ScaleUnits | "
    "SampleRate | StartTime | EndTime"
)


def mseed(*traces, sample_type="i") -> bytes:
    """miniSEED records for ``(seed id, rate, start, samples)`` *traces*."""
    tl = pymseed.MS3TraceList()
    for sid, rate, start, data in traces:
        dtype = np.int32 if sample_type == "i" else np.float32
        tl.add_data(
            pymseed.nslc2sourceid(*sid.split(".")),
            np.asarray(data, dtype=dtype),
            sample_type,
            rate,
            starttime=pd.Timestamp(start).value,
        )
    encoding = (
        pymseed.DataEncoding.STEIM2
        if sample_type == "i"
        else pymseed.DataEncoding.FLOAT32
    )
    return b"".join(tl.generate(encoding=encoding, max_record_length=512))


def response(status: int, content: bytes = b"") -> requests.Response:
    r = requests.Response()
    r.status_code, r._content, r.reason = status, content, "fake"
    return r


class FakeFDSN:
    """Station and dataselect services answering from added channels and traces.

    ``statuses`` queues the status of the next POSTs (default 200), and every
    request is recorded in ``gets`` and ``posts``.
    """

    def __init__(self):
        self.rows: list[str] = []
        self.traces: list[tuple] = []
        self.statuses: list[int] = []
        self.gets: list[dict] = []
        self.posts: list[str] = []

    def channel(
        self,
        sid,
        lat=64.0,
        lon=-151.0,
        rate=50.0,
        scale=4.8876e8,
        units="m/s",
        azimuth=None,
        start="2018-06-07T00:00:00.0000",
        end="",
    ):
        net, sta, loc, cha = sid.split(".")
        if azimuth is None:
            azimuth = {"E": 90.0}.get(cha[-1], 0.0)
        dip = -90.0 if cha.endswith("Z") else 0.0
        self.rows.append(
            f"{net}|{sta}|{loc}|{cha}|{lat}|{lon}|619.0|0.0|{azimuth}|{dip}|"
            f"Trillium|{scale}|0.2|{units}|{rate}|{start}|{end}"
        )

    def trace(self, sid, rate, start, data):
        self.traces.append((sid, rate, start, data))

    def get(self, url, params, timeout):
        self.gets.append(params)
        if not self.rows:
            return response(404)
        text = "\n".join([STATION_HEADER, *self.rows]) + "\n"
        return response(200, text.encode())

    def post(self, url, data, timeout):
        self.posts.append(data)
        status = self.statuses.pop(0) if self.statuses else 200
        if status != 200:
            return response(status)
        wanted = {
            ".".join(line.split()[:4]).replace("--", "") for line in data.splitlines()
        }
        traces = [t for t in self.traces if t[0] in wanted]
        return response(200, mseed(*traces)) if traces else response(204)


@pytest.fixture
def fdsn(monkeypatch):
    """A :class:`FakeFDSN` in place of the FDSN web services."""
    fake = FakeFDSN()
    monkeypatch.setattr(waveforms.requests, "get", fake.get)
    monkeypatch.setattr(waveforms.requests, "post", fake.post)
    monkeypatch.setattr(waveforms.clock, "sleep", lambda s: None)
    return fake
