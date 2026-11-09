# Infrasound waveforms (`infrasound`)

`gaia stage infrasound` stages pressure waveforms (microbarometers and infrasound
sensors) from the same EarthScope FDSN services as [`seismic`](seismic.md), stored the
same way: raw counts at each channel's native rate, one group per rate, next to the
`sites` and `channels` tables. Everything in [seismic.md](seismic.md) applies (sources,
store layout, time, `to_grid`), except how channels are chosen and the units.

Code: `src/gaia_cli/datasets/infrasound.py`, a thin wrapper of `waveforms.load()` with
`arrays=True`.

```bash
gaia stage --config examples/peters-dome.yaml --only seismic,infrasound -o peters-dome
gaia stage infrasound -i aoi.geojson -s 2023-09-14T02:10 -e 2023-09-14T02:20 -n IM --station "I53H?" -o i53
```

## Selecting channels

`channel` (default `BDF,HDF`) is a list in order of preference, as for seismic, but
the choice is made **per sensor position** rather than per station.

### Why not per station

A channel is named `NETWORK.STATION.LOCATION.CHANNEL`, e.g. `AK.CUT.EP.BDF`. Around
Peters Dome the station service returns three kinds of infrasound station, which use
the location code differently:

| kind | example | the location code means |
|---|---|---|
| 1. array, one element per station | IM.I53H1 to I53H8 (IMS array I53US) | nothing special: each element is its own station |
| 2. array, one element per location | AV.KENI, WHTR, BAEI, BAIS: `01` to `06` | a **different sensor**, tens of metres away |
| 3. one sensor, two streams | most AK and some AV stations: BDF at `20` (20 Hz) and `EP` (40 Hz) | **the same sensor**, recorded twice, at the same coordinates |

So "another location" is sometimes another sensor (kind 2) and sometimes another copy
of the same one (kind 3). The metadata doesn't say which, but the coordinates do:
kind 2 elements are at different positions, kind 3 streams at the same one.

peters-dome's `load_infrasound` kept one channel per station, at its lowest location
(the rule `seismic` uses):

- kind 1: all 8 I53 elements kept, which is why its I53 figure is fine;
- kind 2: only element `01` of each array kept, so 18 elements are dropped;
- kind 3: `20` chosen ("20" sorts before "EP"), but EarthScope has the `20` streams'
  metadata and **no data** for them (all 26 were empty for 2023-09-14), so 26
  sensors are skipped as "no data".

### The rule

1. Group the channels by station **and position** (latitude and longitude rounded to
   1e-5°, about a metre).
2. In each group keep the first pattern of `channel` the sensor has (BDF before HDF),
   then the **fastest** stream, then the lowest location.

Kinds 1 and 2 are then one group per element, so every element is kept. Kind 3's two
streams share a group, and `EP` at 40 Hz wins: the native stream rather than the
decimated copy, and the one with data.

On the same station metadata (400 km, 2023-09-14 01:50 to 03:00 UTC, 79 positions at
61 stations), peters-dome's rule gets data from **31 sensors**, this rule from **74**.
The other five have no data under either rule (see [Caveats](#caveats)).

The code is `waveforms.prefer_bands(..., arrays=True)`, and
`test_prefer_bands_with_arrays_keeps_one_stream_per_position` builds one station of
each kind.

### Limits

- Two array elements within about a metre of each other would be merged into one.
- "Fastest" is a heuristic that suits EarthScope's archive today. A sturdier rule would
  fall back to the next stream when the first returns no data, at the cost of a second
  round of requests; not worth it yet.
- `seismic` is unchanged: one band and location per station, so a strong-motion sensor
  at location `10` doesn't join the broadband at `00`.

## Units

Pressure in Pa is counts ÷ sensitivity (`waveforms.to_physical(node, "40hz")`, named
`pressure`). This is exactly obspy's `remove_sensitivity`, which peters-dome used, so
unlike seismic it gives up nothing relative to the reference. Infrasound sensors are
flat over their passband, but their low-frequency corners vary (Hyperion ~0.01 Hz,
Chaparral 64 ~0.03 Hz), so amplitudes below about 0.05 Hz are not reliable.

## Checks

Done on 2026-11-09, staging `examples/peters-dome.yaml`'s infrasound
(400 km, 01:50 to 03:00 UTC):

- **Size and speed.** 79 channels at 61 stations were requested and 74 returned data
  (8 at 20 Hz, 37 at 40 Hz, 6 at 50 Hz, 23 at 100 Hz): 75 MB in memory, 41 MB on disk,
  about 10 s.
- **Against obspy.** For IM.I53H1..BDF, AK.CUT.EP.BDF, AV.KENI.03.HDF and
  AV.SPU..BDF (one per rate), the staged counts equal obspy's raw samples exactly, the
  text station sensitivity equals StationXML's, and counts ÷ sensitivity equals
  `remove_sensitivity`.
- **The event.** On `to_grid(node, rate=20)`, band-passed 0.5 to 5 Hz, each I53
  element's largest pressure between 300 and 1100 s after the origin is 0.61 to
  0.74 Pa, at about 870 s (the main pulse, ~0.29 km/s over 252 km) or 405 s (the
  precursor), as in the peters-dome report.

## Caveats

- Five channels returned no data for this window (AK.L17K, L26K, SSN, AT.DORN,
  AV.IVE); they are listed in `channels` with an empty `group`.
- No array processing (beamforming, back azimuth, trace velocity). The sites' positions
  and `time_offset` are what it needs.
- The low-rate `LDF` (1 Hz) channels aren't requested by default; add them as their own
  run (`channel: LDF`) if wanted.
