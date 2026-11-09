"""HTTP requests to upstream services, with one timeout and retry policy.

Every request has a timeout, so a stalled server can't hang a run, and is
retried with backoff on connection errors, timeouts and transient statuses
(429, 5xx). After the last retry the final response is returned, so
``raise_for_status`` reports it. Imports nothing else from gaia_cli, so any
module can use it.
"""

import threading

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

TIMEOUT = (10, 120)  # seconds to connect, and between bytes read
# Every method: the POSTs gaia makes (FDSN dataselect) are queries
RETRY = Retry(
    total=4,
    backoff_factor=2,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=None,
    raise_on_status=False,
)

_local = threading.local()


def session() -> requests.Session:
    """A new session that retries (see :data:`RETRY`); give it a timeout."""
    s = requests.Session()
    adapter = HTTPAdapter(max_retries=RETRY)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


def _shared() -> requests.Session:
    # One per thread: earthdata.sign resolves URLs from a thread pool
    if not hasattr(_local, "session"):
        _local.session = session()
    return _local.session


def get(url: str, **kwargs) -> requests.Response:
    """``requests.get`` with retries, and :data:`TIMEOUT` unless given."""
    kwargs.setdefault("timeout", TIMEOUT)
    return _shared().get(url, **kwargs)


def post(url: str, **kwargs) -> requests.Response:
    """``requests.post`` with retries, and :data:`TIMEOUT` unless given."""
    kwargs.setdefault("timeout", TIMEOUT)
    return _shared().post(url, **kwargs)
