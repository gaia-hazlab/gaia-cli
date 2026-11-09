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
The comparison under [Checks](#checks) shows the flat band is narrower than assumed
here: amplitudes agree within 3% from 0.067 to 5 Hz, but phase differs at long periods
and amplitude is 3–13% high at 5–8 Hz.

### Not seisfetch either

[seisfetch](https://github.com/Denolle-Lab/seisfetch) (Denolle Lab) also decodes with
pymseed, but PR 2.1 doesn't depend on it (decided 2026-10-21):

- It pins `pymseed<0.10`, and only pymseed 1.0.x has Python 3.15 wheels. It isn't on
  conda-forge.
- Its FDSN station query takes no latitude, longitude or bbox, and drops azimuth and
  dip.
- It has no dataselect POST bulk requests, chunking or retries.
- Its `bundle_to_xarray` outer-joins each channel's own time axis, so mixed rates and
  sub-sample offsets blow up into a mostly-NaN union axis.

Its strength is bulk access to EarthScope's S3 archive, worth revisiting for long
windows. gaia-cli's own layer is `gaia_cli.waveforms` ([../seismic.md](../seismic.md)).

## Checks

Done on 2026-10-15:

- [x] **pymseed has Python 3.15 wheels.** pymseed 1.0.1 on PyPI ships `cp315` wheels
  (`requires_python >=3.10`). conda-forge has `py315` builds of 1.0.1 too, so pixi
  installs it from conda-forge, like gaia-cli's other dependencies (PR 2.1).
- [x] **EarthScope's text station output includes sensitivity.** A
  `format=text&level=channel` query around Peters Dome returns `Scale`, `ScaleFreq` and
  `ScaleUnits` (e.g. `AK.BPAW..BHZ`: 4.8876E8 at 0.2 Hz, m/s), plus azimuth, dip and
  sample rate.

Done in PR 2.1, on 2026-10-21:

- [x] Compare counts ÷ sensitivity with obspy `remove_response(output="VEL")` for a few
  Peters Dome stations in the 0.01–8 Hz band. Done once, in a separate environment
  (`pixi exec --spec obspy`), not as a gaia-cli test. Six BH channels at 50 Hz,
  2023-09-14 01:50–02:30 UTC, staged with `gaia stage seismic`. obspy: demean,
  `remove_response(output="VEL", pre_filt=(0.005, 0.01, 8, 9))`; gaia: (counts −
  mean) ÷ sensitivity; both then band-passed (4th-order Butterworth, zero phase), the
  first and last minute dropped. Each cell is the RMS of the difference ÷ the RMS of
  obspy's velocity, and in brackets the ratio of their RMS amplitudes:

  | channel | 0.01–8 Hz | 15–100 s | 0.067–1 Hz | 1–5 Hz | 5–8 Hz |
  |---|---|---|---|---|---|
  | AK.KTH..BHZ | 0.054 (1.001) | 0.147 (1.000) | 0.044 (1.001) | 0.080 (1.006) | 0.283 (1.082) |
  | AK.CAST..BHZ | 0.099 (1.003) | 0.284 (0.996) | 0.082 (1.000) | 0.116 (1.017) | 0.297 (1.117) |
  | AK.PPLA..BHZ | 0.065 (1.001) | 0.230 (0.999) | 0.061 (1.000) | 0.031 (1.005) | 0.061 (1.028) |
  | AK.BPAW..BHZ | 0.071 (1.003) | 0.132 (1.000) | 0.047 (1.000) | 0.116 (1.014) | 0.288 (1.085) |
  | AK.BPAW..BHN | 0.071 (1.003) | 0.206 (1.001) | 0.044 (1.000) | 0.100 (1.010) | 0.293 (1.088) |
  | AK.CUT..BHZ | 0.153 (1.023) | 0.246 (0.997) | 0.071 (1.000) | 0.151 (1.029) | 0.314 (1.129) |

  The staged counts equal obspy's raw samples exactly for all six, and the text
  station sensitivity equals StationXML's. Correlation over 0.01–8 Hz is 0.989–0.999.
  Amplitudes agree within 3% up to 5 Hz; the differences are mostly phase, largest at
  15–100 s (the sensor's long-period corner), and at 5–8 Hz counts ÷ sensitivity is
  3–13% too large. Good enough for record sections normalized per trace (peters-dome's
  report); not for waveform modeling or amplitudes above 5 Hz, which need option 2 or
  obspy on the staged counts.
