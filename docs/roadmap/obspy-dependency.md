# Sidecar · Can we stage seismic data without obspy?

Status: **decided**, 2026-10-15. gaia-cli does not depend on obspy. Supports
[design goal 2](README.md#design-goals).

## Why it matters

obspy would have been the only dependency forcing a Python pin. It has no 3.15 build
yet, and gaia-cli locks 3.15.0. It also constrains pandas (`<3` at the time of writing,
while gaia-cli resolves pandas 3.0), and it is large compared with what `seismic` and
`infrasound` actually use.

## What gaia-peters-dome uses obspy for

From `src/peters_dome/waveforms.py` (`main` @ `a9401c6`):

| obspy call | purpose | replacement | new dependency? |
|---|---|---|---|
| `Client.get_stations(level="response")` | Station and channel metadata inside the AOI | FDSN station service with `format=text&level=channel`. It returns lat/lon, elevation, azimuth, dip, sample rate and the overall sensitivity (scale) per channel. | none (`requests`) |
| `Client.get_waveforms_bulk` | Waveforms for all picked channels | POST to FDSN dataselect, which returns miniSEED bytes | none |
| `Stream` / miniSEED decoding | Records → arrays | `pymseed` (EarthScope's binding to libmseed; already used in `seisfetch`) | **pymseed** |
| `merge(method=1)`, `split()`, `detrend("demean")` | Join records, keep gaps, remove the mean | numpy | none |
| `filter("lowpass", zerophase=True)` | Anti-alias before interpolation to the common grid | `scipy.signal.butter(..., output="sos")` + `sosfiltfilt` | none (scipy is already in the lockfile, pulled in by other packages) |
| `gps2dist_azimuth` | Distance and back-azimuth from the event | `pyproj.Geod(ellps="WGS84").inv` | none (pyproj is already in the lockfile) |
| `remove_sensitivity` | Infrasound counts → Pa | Divide by the scale from the text station table | none |
| `remove_response(output="VEL", pre_filt=...)` | Seismic counts → m/s | See the options below | — |

Everything except full response removal is a thin layer over HTTP, numpy and scipy.

## Options for response removal

1. **Divide by sensitivity only.** HH/BH channels are broadband sensors whose response
   is flat in velocity from roughly 0.01 Hz to near Nyquist. Peters Dome low-passes to
   8 Hz (0.4 × 20 Hz) over a 40-minute window, so inside that band this is very close to
   `remove_response(output="VEL")`. It drifts at very long periods, and it cannot produce
   `DISP` or `ACC`. This option needs no extra code.
2. **Deconvolve poles and zeros in gaia-cli.** Fetch StationXML once and read the first
   stage's poles and zeros plus the total gain with the standard-library XML parser. Then
   do an FFT deconvolution with scipy, using the same cosine `pre_filt`. That's about
   60–80 lines. It skips the digital FIR stages, which matter little below 0.4·fs. It
   would need tests against obspy output.
3. **Server-side correction.** The IRIS `irisws/timeseries` service accepted
   `correct=true`, `units=VEL`, frequency limits and filtering, and returned corrected
   data. It has not been verified since EarthScope's service migration.

## Decision

Use option 1 and store the raw data:

- Stage **raw counts** with a per-channel `sensitivity` and `sensitivity_units`, and
  derive velocity (m/s) or pressure (Pa) as counts ÷ sensitivity. Record this in the
  `processing` attr.
- Anyone who needs full deconvolution, or displacement or acceleration, can run obspy on
  the staged counts in their own environment, outside gaia-cli. gaia-cli has no optional
  obspy import or `to_stream` helper either.
- This replaces obspy with pymseed, drops the `python <3.15` pin and the pandas
  constraint, and follows the fidelity argument in [03](03-storage.md): keep the raw
  data and derive views from it.

Move to option 2 only if a workflow needs response-corrected data outside the flat band.

## Checks

Done on 2026-10-15:

- [x] **pymseed has Python 3.15 wheels.** pymseed 1.0.1 on PyPI ships `cp315` wheels
  (`requires_python >=3.10`). Add it as a PyPI dependency in pixi; conda-forge is not
  needed.
- [x] **EarthScope's text station output includes sensitivity.** A
  `format=text&level=channel` query around Peters Dome returns `Scale`, `ScaleFreq` and
  `ScaleUnits` (e.g. `AK.BPAW..BHZ`: 4.8876E8 at 0.2 Hz, m/s), plus azimuth, dip and
  sample rate.

Left for PR 2.1:

- [ ] Compare counts ÷ sensitivity with obspy `remove_response(output="VEL")` for a few
  Peters Dome stations in the 0.01–8 Hz band. Do this once, in a separate environment,
  not as a gaia-cli test.
