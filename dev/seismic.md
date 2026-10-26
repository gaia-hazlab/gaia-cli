# Seismic waveforms (`seismic`)

`gaia stage seismic` stages raw seismic waveforms from EarthScope's FDSN web services
for the stations inside an AOI. Waveforms are stored as **raw digitizer counts at each
channel's native sample rate**, with the channel metadata needed to convert them to
physical units. Nothing is demeaned, filtered or resampled at stage time. A common grid
and physical units are read-side views (`gaia_cli.waveforms.to_grid`, `to_physical`),
like `stations.to_grid` for station reports.

Code: `src/gaia_cli/waveforms.py` (shared with `infrasound`, roadmap step 2.2) and
`src/gaia_cli/datasets/seismic.py`.

```bash
gaia stage --config examples/peters-dome-seismic.yaml -o peters-dome-seismic
gaia stage seismic -i aoi.geojson -s 2023-09-14T02:06 -e 2023-09-14T02:08 -n AK --station BPAW -o bpaw
```

## Sources

| service | used for |
|---|---|
| `https://service.earthscope.org/fdsnws/station/1/query` | Channels in the AOI: `level=channel&format=text` gives coordinates, elevation, depth, azimuth, dip, the overall sensitivity (`Scale`, `ScaleFreq`, `ScaleUnits`), the sample rate and each epoch's start and end. |
| `https://service.earthscope.org/fdsnws/dataselect/1/query` | miniSEED, by POST (one line per channel), decoded with [pymseed](https://github.com/EarthScope/pymseed) (libmseed). |

EarthScope serves networks worldwide, so `EXTENT` is global with an open time span.
Data are open for most networks, but each network has its own terms and citation
(network DOIs at <https://www.fdsn.org/networks/>). Each group lists its `networks`.
Restricted channels aren't requested (`includerestricted=false`).

## Selecting channels

Options are FDSN codes, with `?` and `*` wildcards and commas between alternatives:
`network`, `station`, `location` (`--` is the empty code) and `channel` (default
`HH?,BH?`). The AOI does the spatial selection: a 250 km circle
(`{latitude, longitude, radius_km: 250}`) replaces peters-dome's `max_distance_km`.

1. The station service is queried with the AOI's WGS84 bounds and the window. Channels
   whose station point is outside the AOI polygon are dropped.
2. **Epochs.** A channel can have several epochs (a sensor swap, a new sensitivity).
   gaia keeps the epoch covering the window's start, else the first one that overlaps
   the window, and logs a warning when the sensitivity changes inside the window.
3. **Band preference.** `channel` is a list in order of preference. Each station keeps
   every channel of the **first pattern it has**, of one band (the channel code's first
   two letters) and of its **lowest location code**. With `HH?,BH?`, a station with both
   keeps HHZ, HHN and HHE only, and a BH-only station keeps its BH channels. This mirrors
   peters-dome's `_pick_channels`, and **differs from FDSN's own semantics**, where
   `channel=HH?,BH?` is the union of both. To keep both bands, stage them as two runs
   (`channel: "HH?"` and `channel: "BH?"`) into separate stores, or ask for a plain-union
   option.
4. **Size guard.** Before downloading, gaia estimates the counts' size, channels ×
   rate × duration × 4 bytes, and refuses requests over 2 GB
   (`waveforms.MAX_BYTES`). For scale: Peters Dome's 157 channels over 40 min are 82 MB.

Dataselect requests go 50 stations at a time, each retried up to 3 times with backoff
on 5xx errors and timeouts. No data (204 or 404) is fine, and a request that still
fails is logged and skipped, so the other stations are kept (as in peters-dome).

## Store layout

```
waveforms/seismic/            provenance attributes
  sites/      one row per position (NET.STA.LOC), along an xvec geometry dimension
  channels/   one row per channel requested (NET.STA.LOC.CHA), with or without data
  40hz/       counts (channel, time) int32, site and time_offset along channel
  50hz/
  100hz/
```

The metadata is normalized into two small tables, so each fact is stored once, and the
counts refer to them by id:

- **`sites`** has one row per position, `NET.STA.LOC` (`NET.STA` for an empty
  location, e.g. `AK.BPAW`). It is laid out like a `stations/*` group, along an xvec
  `geometry` dimension, so xvec CF-encodes and decodes it with its index, and the STAC
  Item's bbox comes from it. A station's three components share one site; an
  infrasound array's elements (2.2) are separate sites, by location code.
- **`channels`** has one row per channel requested, including channels that returned
  no data (their `group` is `""`), so the request is recorded in full.
- **One group per sample rate**, named `<rate>hz` (`%g`: `50hz`, `0.1hz`, `12.5hz`),
  holds only the counts and what is needed to place them in time. Rates come from the
  **data**, not the metadata, with float jitter rounded (within 1e-4 of a whole number
  of Hz, else to 1e-6 Hz). A channel whose data mixes rates keeps the rate with most
  samples, with a warning.

`sites` (dimension `geometry`):

| name | notes |
|---|---|
| `geometry` | shapely Point (lon, lat), EPSG:4326, xvec-encoded in the store; a site's first channel (by SEED id) gives its position, logged if its channels differ |
| `site` | `NET.STA.LOC` |
| `network`, `station`, `location` | the SEED codes; `location` is often `""` |
| `elevation` | m |
| `distance_km`, `back_azimuth` | from the origin, on WGS84 (pyproj `Geod`); back azimuth is from the site to the origin, clockwise from north |

`channels` (dimension `channel`):

| name | notes |
|---|---|
| `channel` | SEED id `NET.STA.LOC.CHA`, e.g. `AK.BPAW..BHZ` |
| `site` | the channel's row in `sites` |
| `code` | channel code, e.g. `BHZ` |
| `depth` | m below the site |
| `azimuth`, `dip` | degrees |
| `sensitivity`, `sensitivity_frequency`, `sensitivity_units` | overall sensitivity in counts per `sensitivity_units` (lowercased, `m/s` for seismometers, `m/s**2` for accelerometers) at a frequency in Hz |
| `sample_rate` | Hz, of the data (of the metadata for a channel without data) |
| `group` | the rate group holding its counts, `""` if it returned no data |

Each rate group (dimensions `channel`, `time`):

| name | notes |
|---|---|
| `counts` | `(channel, time)` int32 raw counts; gaps hold the fill value −2³¹, which xarray masks to NaN (as float64) on read. `xr.open_datatree(..., mask_and_scale=False)` keeps int32. Chunks: one channel by up to an hour of samples. |
| `site` | `(channel)`, the channel's row in `sites` |
| `time_offset` | `(channel)` timedelta, see below |

Attributes on `waveforms/seismic` itself: the provenance fields (`title`, `license`,
`providers`, `links` to both services, `start_datetime`/`end_datetime`), `processing`,
`networks` (those with data, to cite), `gaia:origin` (`latitude`, `longitude` and
`source`: the stage file's origin or the AOI centroid) and `gaia:query` (the service
URL and FDSN codes).

To look up a channel's metadata, select the tables by id:

```python
node = tree["waveforms/seismic"]
counts = node["50hz"].counts
meta = node["channels"].to_dataset().sel(channel=counts.channel)
sites = node["sites"].to_dataset().xvec.decode_cf()
```

The helpers below do these joins.

## Time

A seismic window is read **exactly**: `PAD = 0`, unlike the one-day padding of other
timed datasets. `load()` gets the window's start and an inclusive end one microsecond
before the window's end, so `[start, end)` keeps every sample, the last second
included. Data are fetched with a second's margin on each side and trimmed.

Each rate's `time` axis starts at the window's start, one step (1/rate) per sample, in
integer nanoseconds. Digitizer clocks aren't aligned to the second: AK.BPAW's samples
fall 8.4 ms after each 20 ms step. Each channel's samples go at the **nearest** step,
and `time_offset` records the shift, under half a step, so **sample *i* was recorded at
`time[i] + time_offset`**. It is the offset of the channel's longest segment; a segment
after a gap may differ slightly. Gaps hold the fill value; where segments overlap, the
later one wins (logged).

## Units

Physical units are `counts / sensitivity` (`waveforms.to_physical(node, "50hz")`), velocity in m/s for
the HH and BH broadband channels. This assumes a flat response at the sensitivity's
frequency, which is not a full deconvolution (see
[roadmap/obspy-dependency.md](roadmap/obspy-dependency.md)). Compared with obspy
`remove_response(output="VEL")` for six Peters Dome channels, the counts are identical,
the RMS amplitude agrees within 3% from 0.067 to 5 Hz, and the waveforms differ by 3–15%
RMS in that band. The rest is phase: at 15–100 s periods the amplitude still agrees
(within 0.4%) but the waveforms differ by 13–28%, and at 5–8 Hz counts / sensitivity
overestimates the amplitude by 3–13%. For response-corrected data, run obspy on the
staged counts in your own environment.

## Common grid

`waveforms.to_grid(tree["waveforms/seismic"], rate=20)` merges every rate group onto one
`(channel, time)` grid: a zero-phase 4th-order Butterworth low-pass at 0.4 × rate
(`scipy.signal.sosfiltfilt`, per contiguous run of samples, only for channels sampled
faster than `rate`), then linear interpolation from each sample's true time. Gaps stay
NaN. Channels without data are left out, and the columns of `channels` and `sites`
(but not `geometry`) become coordinates along `channel`. `waveforms.by_component(grid,
node["sites"])` reshapes it to `(geometry, component, time)` with Z/N/E (and 1/2, read
as N/E within 5° of 0°/90°), the layout peters-dome's report plots, with the sites'
variables along `geometry` as `stations.to_grid` gives for station reports:

```python
from gaia_cli import io, waveforms

tree = io.open_datatree("peters-dome-seismic/peters-dome-seismic.zarr")
node = tree["waveforms/seismic"]
grid = waveforms.to_grid(node, rate=20)  # velocity, m/s
z = waveforms.by_component(grid, node["sites"]).velocity.sel(component="Z")
z = z.sortby("distance_km")
```

## Why not seisfetch or obspy

obspy: see [roadmap/obspy-dependency.md](roadmap/obspy-dependency.md) (no Python 3.15
build, a pandas constraint, and large for what is used).

[seisfetch](https://github.com/Denolle-Lab/seisfetch) wraps pymseed and FDSN too, but:

- it pins `pymseed<0.10`, and only pymseed 1.0.x has Python 3.15 wheels (gaia-cli locks
  3.15); seisfetch isn't on conda-forge;
- its FDSN station query takes no latitude/longitude/bbox, and drops azimuth and dip;
- it has no dataselect POST bulk requests, chunking or retries;
- its `bundle_to_xarray` outer-joins each channel's own time axis, so mixed rates and
  sub-sample offsets blow up into a mostly-NaN union axis.

Its real strength is bulk access to EarthScope's S3 archive, worth revisiting for long
windows.

## Caveats

- Channels with non-integer samples (rare for seismic data) are skipped with a warning:
  counts are int32.
- One location per station: arrays whose elements are location codes (common for
  infrasound) need another rule in step 2.2.
- A whole station's request is skipped if it keeps failing after retries; the run logs
  it but doesn't fail.
- The STAC Item's bbox covers the stations, not the AOI.
