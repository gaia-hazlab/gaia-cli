# HTTP retries and timeouts

gaia makes remote requests in two ways, and they retry differently:

- **`fetch.py`** sends every request gaia makes itself (APIs, searches, signing). It
  uses requests and urllib3.
- **GDAL** reads rasters (COGs, VRTs, HDF5) through `/vsicurl/`. odc-stac, rioxarray
  and rasterio all read this way.

GDAL's settings are in `io.GDAL_HTTP`. They're set as environment variables when
`gaia_cli.io` is imported, unless a variable is already set, so you can override any
of them. Because they're in the environment, they cover all three readers, in the
main process and in worker processes (`io.compute_in_processes`). odc-stac's
`cloud_defaults` isn't used, because its 10 retries would replace these.

As of 2026-11-10, with GDAL 3.13.3:

| | `fetch.py` | GDAL default | GDAL in gaia |
|---|---|---|---|
| Retries | 4 | 0 | 5 |
| Delays | 0, 4, 8, 16 s (~28 s in all) | 30 s, then ×2–2.5 each time | 0.5 s, then ×2–2.5 each time (~20 s in all) |
| Statuses retried | 429, 500, 502–504 | 429, 500, 502–504, and a 400 with "RequestTimeout" | same as default |
| Connection failures retried | any connection, read or timeout error | only these curl error messages: "Connection timed out", "Operation timed out", "Connection reset by peer", "Connection was reset", "SSL connection timeout" | same as default |
| `Retry-After` | honored (429, 503) | ignored | ignored |
| Timeouts | 10 s to connect, 120 s between bytes | none | 10 s to connect, 120 s for each request |

GDAL's delays grow by a random factor of 2–2.5 each time, with no cap. At odc-stac's
10 retries, a dead server takes ~20 minutes to fail. A 503 test URL failed after
21 s with gaia's settings.

GDAL's `GDAL_HTTP_TIMEOUT` limits a whole request, unlike urllib3's read timeout.
gaia's GDAL reads are range requests of a few MB at most, so 120 s is ample. GDAL
retries its timeouts, so a stalled server becomes a retry instead of hanging the run.

## Failures GDAL doesn't retry

Retries are chosen in `CPLHTTPGetNewRetryDelay()` (`port/cpl_http.cpp`). GDAL matches
curl's error text, not its error code. These fail on the first attempt:

- an empty reply (response code 0)
- a refused connection
- a DNS failure
- an SSL error other than a timeout
- a 403 or 404

The HTTP/1.1 setting in `io.py` works around one of these, where CloudFront drops an
HTTP/2 stream and GDAL sees code 0.

`GDAL_HTTP_RETRY_CODES` (GDAL ≥ 3.10) **replaces** that list, it doesn't add to it.
With `ALL`, GDAL retries every non-200 response, 403 and 404 included. It also stops
retrying the connection timeouts and resets above. Don't use it as a catch-all.

## Reported errors

rasterio raises `RasterioIOError: Read failed. See previous exception for details.`
GDAL's own error is attached as the exception's cause. `gaia stage` reports each
failed dataset with its causes, so the actual error is visible: an HTTP status or a
curl message.

On 2026-11-10, CI's online test of PRISM failed once with only the rasterio message.
The file was `prism/us/an/800m/tmax/daily/2025/prism_tmax_us_30s_20251211.tif` on
Kopah, and it read fine afterwards. With odc-stac's 10 retries in effect, the error
was probably one GDAL doesn't retry.
