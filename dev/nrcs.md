# NRCS SNOTEL and SCAN

`gaia stage nrcs` stages hourly station reports from the NRCS
[Air and Water Database (AWDB) REST API](https://wcc.sc.egov.usda.gov/awdbRestApi/swagger-ui/index.html),
public and without a token. Its main use is soil moisture: see
[soil-moisture.md](soil-moisture.md) for why NRCS rather than Synoptic, and for the
other sources.

```bash
gaia stage nrcs -i wa.geojson -s 2026-10-01 -e 2026-10-31 -o wa-soil             # SMS
gaia stage nrcs -i wa.geojson -s 2026-01-01 -e 2026-03-31 -v SMS STO WTEQ -o wa-snow
```

The group `stations/nrcs` has the same CF ragged layout as `stations/synoptic` (see
[stations.md](stations.md)), so `stations.to_grid` and `to_dataframe` read it.

## Coverage

Networks: SNOTEL (`SNTL`), SNOTEL-lite (`SNTLT`), SCAN and Tribal SCAN (`CSCAN`).
Stations with soil moisture (`SMS`) on 2026-11-09, active only:

| network | stations | where |
|---|---|---|
| SNOTEL | 523 | western mountains: UT 141, NV 58, CO 53, ID 49, AK 41, WY 38, OR 32, CA 32, MT 25, WA 24, NM 18, AZ 11, SD 1 |
| SCAN | 210 | 40 states plus HI, PR, VI; mostly agricultural sites. AK 3, WA 2, OR 0 |
| SNOTEL-lite | 6 | CO, NV, UT |

Counting closed stations too, 774 have had soil moisture. Closed stations are
staged when their record overlaps the window. Soil sensors at SNOTEL sites were mostly
added from 2000 to 2020, so a site's snow record is often much longer than its soil
record.

Each station element has a `heightDepth` in inches, negative below ground. SNOTEL
and SCAN soil sensors are at 2, 4, 8 and 20 in, and 40 in at some sites. SCAN
also has two soil temperature sensors at each depth (ordinals 1 and 2).

## What is staged

- One request lists every station of the four networks with any of the
  `--variables` (`elements=SMS:*`, `activeOnly=false`), about 2 s. Stations are
  clipped to the AOI polygon and to those open at some time in the window.
- Then `/data` with `duration=HOURLY`, 10 stations per request. Data AWDB marks
  suspect are left out (AWDB's default, `returnSuspectData=false`). QC and QA flags
  aren't stored.
- **Names**: the AWDB element code, plus the depth in centimetres, rounded, for a
  sensor below ground: 2, 4, 8, 20 and 40 in become `SMS_5cm`, `SMS_10cm`, `SMS_20cm`,
  `SMS_51cm` and `SMS_102cm`. The exact depth is the variable's `depth` attribute
  (`"5.08 cm"`). A second sensor at the same depth gets its ordinal (`STO_5cm_2`), and
  one above ground gets `_above` (`<code>_305cm_above` at 120 in). Elements without a depth keep
  the code (`WTEQ`). One `obs` variable per depth keeps the stations helpers
  unchanged; a station without a sensor at that depth has NaN.
- **Units**: SI. Soil moisture (`SMS`, and `SMV`/`SMN`/`SMX`) is AWDB's percent / 100,
  m3 m-3. Inches become mm, degF degC and mph m s-1. Other units are AWDB's, in the
  `units` attribute. `standard_name` is set for `SMS`, `STO`, `TOBS`, `WTEQ`, `SNWD`
  and `PREC`.
- Station variables: `station` (the AWDB triplet, `679:WA:SNTL`), `name`, `elevation`
  (m, from ft), `network_id`, `utc_offset`.

### Time

AWDB's hourly dates are the station's **local standard time all year**, a fixed offset
from UTC given by the station's `dataTimeZone` (−8 in WA and OR, −9 in AK). `gaia`
converts them to UTC and keeps the offset as `utc_offset`. The request starts 12 hours
before the window so that every UTC hour in it is covered, and the reports are trimmed
to the window afterwards.

Checked on 2026-11-09 at Cook Farm Field D (2198:WA:SCAN), which Synoptic also serves
(COOW1, in UTC): AWDB's `TOBS` + 8 h matches Synoptic's `air_temp` to 0.09 °C on
average over 83 hours, against 1.2 °C or more for an offset one hour either way.

`SMS` and `TOBS` are instantaneous readings (AWDB function code `C`). `PREC` is the
running total since the start of the water year (October 1), so it drops to zero once
a year. Use `to_grid(how="sum")` on it only within one water year.

## Caveats

- **Sensor and calibration**: SNOTEL and SCAN use Stevens HydraProbes, with a generic
  calibration. Absolute values differ between soils; wetting and drying, anomalies and
  percentiles compare better than raw values. Lind #1 (2021:WA:SCAN) reads about 0.00
  at 5 cm through October 2026.
- **Daily values aren't staged**: AWDB's daily `SMS` is derived. `to_grid(freq="1D")`
  gives the last hourly report of each day instead.
- **Request size**: AWDB refuses a station query that asks for every station's element
  list at once ("number of station elements exceeds the amount allowed"), so `nrcs.py`
  doesn't ask for it. Long windows over many stations make large `/data` responses:
  about 5 depths × 8,760 hours per station per year.

## Follow-ups

- Store QC flags (`returnFlags=true`), and optionally suspect data
  (`returnSuspectData=true`), which the soil-moisture-fusion validation keeps.
- A common name for volumetric soil moisture, shared with ERA5's `swvl1` and `swvl2`, once
  `gaia align` handles station depths.
