import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from gaia_cli import fetch


@pytest.fixture
def server():
    """A local server answering each request with the next of ``statuses``
    (then 200), recording each request's method in ``seen``."""
    statuses: list[int] = []
    seen: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def answer(self):
            seen.append(self.command)
            n = int(self.headers.get("Content-Length") or 0)
            self.rfile.read(n)
            self.send_response(statuses.pop(0) if statuses else 200)
            self.send_header("Content-Length", "2")
            self.end_headers()
            self.wfile.write(b"ok")

        do_GET = do_POST = answer

        def log_message(self, *args):
            pass

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{httpd.server_port}/"
    yield url, statuses, seen
    httpd.shutdown()


def test_get_and_post_retry_transient_errors(server):
    url, statuses, seen = server
    statuses += [503]
    assert fetch.get(url).status_code == 200
    statuses += [502]
    assert fetch.post(url, data="query").status_code == 200  # dataselect POSTs
    assert seen == ["GET", "GET", "POST", "POST"]


def test_client_errors_are_returned_without_retrying(server):
    url, statuses, seen = server
    statuses += [404]
    assert fetch.get(url).status_code == 404
    assert seen == ["GET"]


def test_requests_have_a_timeout_by_default(monkeypatch):
    sent = {}
    monkeypatch.setattr(
        fetch.requests.Session, "request", lambda self, *a, **k: sent.update(k)
    )
    fetch.get("https://example.com")
    assert sent["timeout"] == fetch.TIMEOUT
    fetch.post("https://example.com", timeout=5)
    assert sent["timeout"] == 5
