# Soil moisture station timeseries

Status: step 1 (Synoptic) implemented on branch `add-soil-moisture`; step 2 (NRCS) proposed.
Survey of Washington stations made 2026-10-30.

Goal: for any AOI in Washington, stage soil moisture timeseries from stations in the
same way `gaia stage synoptic` stages precipitation: every report, as a CF ragged array
(see [stations.md](stations.md)).

## Summary

- Synoptic has 28 WA stations whose metadata lists `soil_moisture` (24 active, 4
  inactive) and 13 with `soil_moisture_tension`. Only **8 actually report soil moisture
  now**, and only 4 of them (CW3E) say at which depth. Synoptic carries **no SNOTEL
  soil moisture**.
- NRCS's AWDB web service needs **no credentials**. It has soil moisture at **24 active
  WA SNOTEL stations and 2 SCAN stations**, hourly, at labelled depths (2, 4, 8, 20 and
  sometimes 40 in), with records back to 2000–2012. That is the better source for
  Washington, and step 2 below adds it.
- Free Synoptic accounts can only reach about the last year of history. Our token was
  granted longer access on 2026-11-02 (see [Caveats](#caveats)).

## What Synoptic provides in Washington

Found with `/stations/metadata?state=WA&vars=soil_moisture&sensorvars=1` (`status`
active and inactive), then a timeseries request for 2026-10-15 to 2026-10-30.

### Stations reporting now

| station | network | sensors (depth) | report interval | soil record from |
|---|---|---|---|---|
| C3MYD Maywood | CW3E | 6: 5, 10, 15, 20, 50, 100 cm | hourly | 2025-08 |
| C3STP Stampede | CW3E | 6: 5, 10, 15, 20, 50, 100 cm | hourly | 2026-10 |
| PVC54 Paradise | NWAC | 2, depth not given | hourly | 2023-07 |
| SUN64 Sunrise Base | NWAC | 2, depth not given | hourly | 2020-11 |
| OHAW1 Ohanapecosh | RAWS | 1, depth not given | hourly (:41) | 2016-11 |
| HZCW1 Hozomeen | RAWS | 1, depth not given | hourly (:23) | 2025-09 |
| COOW1 Cook Farm Field D | SCAN | 1 of SCAN's 5, depth not given | hourly, about half the hours missing | 2017-08 |
| LDDW1 Lind #1 | SCAN | 1 of SCAN's 5, depth not given | hourly, gaps; reads 0.0 for the whole window | 2017-08 |

Synoptic's soil sensors report hourly in Washington. Just over the border in Oregon,
some Columbia Plateau Conservation Research Center (CPC) stations report every 30
minutes, and SKYPN (Portland) about every 10 minutes. CW3E logs every 2 minutes, but
Synoptic serves hourly values.

### Listed but not usable

- **HADS (9) and CRN (2)**: Synoptic's station status is "active", but the soil sensors
  only reported in 2006–2008. Station status describes the station, not each sensor:
  check `SENSOR_VARIABLES.<var>.<set>.PERIOD_OF_RECORD` instead.
- **RAWS `soil_moisture_tension` (13)**: one per station. Twelve ended between 2007 and
  2024. BMFW1's period of record runs to today, but the timeseries has no values.
- **RAWS `soil_moisture` at FTFW1, PMFW1, TSOW1**: ended 2011–2024.
- **C3HDW Headworks, C3YLI YLine (CW3E)**: listed with 6 depths, but timeseries
  requests return "No stations found … or your account does not have access". Probably
  restricted by the provider.
- **USGS Hydro (3, inactive)**: a few hours in 2017.

### Instruments: how they differ

Two kinds of quantity are reported:

- **`soil_moisture`** (`%`): volumetric water content, the volume of water per volume
  of soil × 100. Dielectric sensors measure it indirectly. Water's dielectric
  permittivity (~80) is far higher than soil minerals' (~4) or air's (1), so the
  sensor's electrical response tracks water content. A calibration curve, often
  generic (loam), turns that response into water content.
  - **Stevens HydraProbe** (coaxial impedance): CW3E, at 5, 10, 15, 20, 50 and 100 cm
    ([CW3E surface met](https://cw3e.ucsd.edu/cw3e_observations_surfacemet/)), and
    NRCS SCAN and SNOTEL. It also reports soil temperature and electrical
    conductivity.
  - **Campbell CS616/CS650** (water content reflectometers): some older CW3E stations,
    and offered in fire-weather (RAWS) station packages.
  - Synoptic does not say which sensor NWAC or RAWS stations use.

  Different sensors, calibrations and soils mean absolute values are **not directly
  comparable between networks**, or even between stations. Wetting and drying events,
  anomalies and percentiles compare better than raw percentages. Saturation in a
  forest loam may read ~40 %, in a sand ~25 %.
- **`soil_moisture_tension`** (`centibars` = kPa): matric potential, how tightly the soil
  holds its water, which is what plants and drainage "feel". Granular matrix sensors
  such as the Irrometer **Watermark** measure it as electrical resistance through a
  porous block in equilibrium with the soil: 0 cb is saturated, ~10–30 cb field
  capacity, and 200 cb the sensor's dry limit. They are temperature-sensitive. Turning
  tension into water content needs the soil's **water retention curve** (for example
  van Genuchten parameters from texture), which is site-specific and approximate.
  People do it, but it is a model, not a measurement. In Washington this is academic:
  no tension sensor on Synoptic reports now.

## Step 1: Synoptic, every sensor set (this branch)

`gaia stage synoptic -v soil_moisture` already ran, but kept only the first sensor set
(`soil_moisture_set_1`) and dropped depth. Now:

- The timeseries request asks for `sensorvars=1`, which gives each set's `position`
  (metres, negative below ground). The GeoJSON metadata output does not include it.
- Every sensor set of each requested variable is kept, as its own `obs` variable:
  - a set below ground is named for its depth: `soil_moisture_5cm`, …,
    `soil_moisture_100cm`;
  - otherwise set 1 keeps the variable's name and set *n* is `<var>_set_n` (NWAC's
    second sensor is `soil_moisture_set_2`). A second set at an already-named depth also
    keeps its set number.
- This applies to every variable, so a station with two `air_temp` sets now also gets
  `air_temp_set_2`. Derived sets (`_set_1d`) are used only when there is no measured
  set, as before.

No new dimension: one `obs` variable per depth keeps `stations.to_grid` and
`to_dataframe` unchanged, and stations without a sensor at that depth have NaN. Depths
are free-form, so a large AOI collects many variables (a WA bounding box gave
`_5cm`, `_7cm`, `_10cm`, `_15cm`, `_20cm`, `_30cm`, `_45cm`, `_50cm`, `_75cm`,
`_100cm`, `_105cm`, mostly from Oregon stations).

```bash
gaia stage synoptic -i wa.geojson -s 2026-10-22 -e 2026-10-29 -v soil_moisture -o wa-soil
```

```python
ds = xr.open_dataset("wa-soil/wa-soil.zarr", engine="zarr", group="stations/synoptic")
daily = stations.to_grid(ds, freq="1D")       # last report of each day
daily.soil_moisture_10cm
```

## Step 2: NRCS AWDB, SNOTEL and SCAN (proposed, separate PR)

The [AWDB REST API](https://wcc.sc.egov.usda.gov/awdbRestApi/swagger-ui/index.html)
is public: no token, JSON.

```
GET https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/stations
      ?stationTriplets=*:WA:SNTL,*:WA:SCAN&elements=SMS:*&returnStationElements=true&activeOnly=true
GET https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data
      ?stationTriplets=679:WA:SNTL&elements=SMS:*&duration=HOURLY&beginDate=2026-10-28 00:00&endDate=2026-10-29 00:00
```

- `SMS` is soil moisture (`pct`, volumetric). Each station element has a `heightDepth`
  in inches (negative is below ground) and a `beginDate`. 24 of the 79 active WA SNOTEL
  stations have it, at depths of {-2, -4, -8, -20} in plus -40 at some, from 2000–2020.
  Both WA SCAN stations (Lind #1, Cook Farm Field D) have -2, -4, -8, -20 and -40 in.
  Durations: HOURLY, DAILY, SEMIMONTHLY, MONTHLY.
- Proposed: `datasets/nrcs.py`, `KIND = "stations"`, group `stations/nrcs`, same ragged
  layout and helpers. Variables by AWDB element code (`SMS`, also `STO` soil temperature,
  `WTEQ` snow water equivalent, `PREC`), with depth naming as in step 1 (see
  [Decisions](#decisions)). AWDB stores inches and °F; convert to mm and °C. Query stations by
  state(s) covering the AOI bounds, then clip to the AOI as `synoptic` does.
- AWDB's hourly times are local standard time (PST, UTC−8, all year), so convert to
  UTC. Check this against a station before relying on it.

USCRN (Marblemount, Quinault) publishes hourly soil moisture at 5, 10, 20, 50 and 100
cm as public files from NCEI. That is a possible step 3: only 2 WA stations.

## Caveats

- **History**: Synoptic refuses any request entirely more than about a year old ("does
  not have access to the requested history") and silently truncates one that straddles
  the limit. A 2025-01…2025-12 request returned only Nov–Dec. Synoptic's FAQ gives free
  accounts one year; the academic Open Access program can grant more on request, and
  granted it to our token on 2026-11-02: both requests below now start where asked.
  Other free tokens still hit the limit, which affects precipitation staging too, so
  we should warn when the first report is much later than the requested start
  (follow-up). On 2026-10-30, before the grant, with OHAW1 (soil moisture since
  2016-11):

  ```bash
  B="https://api.synopticdata.com/v2/stations/timeseries?token=$SYNOPTIC_TOKEN&stid=OHAW1&vars=soil_moisture&obtimezone=UTC"
  # Entirely older than a year: refused
  curl -s "$B&start=202401010000&end=202401020000" | jq -c '.SUMMARY | {RESPONSE_CODE, RESPONSE_MESSAGE}'
  # {"RESPONSE_CODE":403,"RESPONSE_MESSAGE":"Account associated with this token does not have access to the requested history. ..."}
  # Straddling the limit: "OK", but starts exactly one year before now, not 2025-10-01
  curl -s "$B&start=202510010000&end=202511300000" | jq -c '{summary: .SUMMARY.RESPONSE_MESSAGE, first: .STATION[0].OBSERVATIONS.date_time[0]}'
  # {"summary":"OK","first":"2025-10-30T17:41:00Z"}
  ```
- **Active stations only**: `get_station_metadata` asks for `status=active`, so a
  historical window misses stations that have since closed (existing behaviour).
- **No QC**: Synoptic QC flags are not stored. LDDW1 reads 0.0 for weeks.
- **Earlier 502s** (2026-10-30) were Synoptic gateway errors that recovered within
  minutes, not throttling. Synoptic reports limits in its JSON `SUMMARY`.

## Decisions

- SI units throughout. Soil moisture is staged as a volume fraction, m3 m-3 (Synoptic's
  and AWDB's percent / 100), and tension in kPa. AWDB depths are inches, so step 2
  converts them to centimetres:
  2, 4, 8, 20 and 40 in become `soil_moisture_5cm`, `_10cm`, `_20cm`, `_51cm` and
  `_102cm`, rounded to the centimetre like Synoptic's names, with the exact depth in
  metres as an attribute.

## Open questions

- Which sensor sits at which depth at NWAC, RAWS and Synoptic's SCAN feed? Ask NWAC and
  compare SCAN's single Synoptic series with AWDB's depths.
