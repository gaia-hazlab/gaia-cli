# Soil moisture sources

Where to get soil moisture for an AOI, and what gaia stages. Surveyed for Washington
on 2026-10-30 and extended to Oregon and Alaska on 2026-11-09.

**`gaia stage nrcs` is the station source** ([nrcs.md](nrcs.md)): NRCS SNOTEL and SCAN,
hourly, at labelled depths, no token, records back to 2000–2012 in Washington, and
coverage across the western US and Alaska. ERA5 (`gaia stage era5`) is the gridded one.
The other sources below are documented, not staged.

## Summary

| source | what | WA / OR / AK | depths | access | in gaia |
|---|---|---|---|---|---|
| **NRCS SNOTEL, SCAN** (AWDB) | in situ, hourly | 26 / 32 / 44 active stations | 5, 10, 20, 51, (102) cm | public REST API | **`nrcs`** |
| Synoptic | in situ, hourly, many networks | 8 reporting in WA | known at CW3E only | `SYNOPTIC_TOKEN`; about one year of history on free accounts | `synoptic -v soil_moisture`: first sensor only, no depth |
| NOAA USCRN | in situ, hourly | 3 / 4 / 24 stations | 5, 10, 20, 50, 100 cm | public files (NCEI) | no |
| USBR AgriMet | in situ, agricultural, Pacific Northwest | WA, OR, ID, MT | labelled, in inches | public | no |
| WSU AgWeatherNet | in situ, agricultural, mostly eastern WA | WA | 8 in (20 cm) | public API | no |
| AmeriFlux | flux towers, soil water content (SWC) | a few | site-specific | free account, CC-BY-4.0 | no |
| Cook Farm (USDA ARS) | dense research probe network, Pullman, WA | WA | 30–150 cm | figshare archive | no |
| ISMN | harmonized archive of many networks (SCAN, SNOTEL, USCRN, …) | global | per network | free account | no |
| ERA5 | reanalysis, 0.25°, hourly | global | layers 0–7, 7–28 cm (and 28–100, 100–289 cm, not staged) | `ARRAYLAKE_TOKEN` | **`era5`** (`swvl1`, `swvl2`) |
| SMAP L3/L4 | satellite, 9–36 km; L4 adds root zone | global | ~5 cm; L4 0–100 cm | NASA Earthdata | no (roadmap, later) |
| NISAR SME2 | satellite, 200 m, provisional | global | ~5 cm | NASA Earthdata | no |

### How soil-moisture-fusion gets its data

The [soil-moisture-fusion](https://github.com/meganfrisella/soil-moisture-fusion)
project (validating NISAR SME2 in eastern Washington) calls the providers directly
(`src/insitu/coverage.py`), caching and checksumming each response:

- AWDB `/data`, `SMS:-2:1` (the 2 in, 5.08 cm sensor) hourly at 2 SCAN and 2 SNOTEL
  stations, percent / 100, times from each station's fixed offset to UTC (as
  `nrcs.py` does), keeping NRCS suspect values;
- NOAA USCRN hourly02 files (`SOIL_MOISTURE_5`, already m3 m-3), Spokane 17 SSW;
- USBR AgriMet (`XSM<n>` channels, percent at labelled depths);
- for its station map only: AgWeatherNet, AmeriFlux and the Cook Farm archive.

`gaia stage nrcs -v SMS` covers its NRCS part, at every depth.

## NRCS vs Synoptic in Washington

Synoptic lists 28 WA stations with `soil_moisture` (24 active, 4 inactive) and 13 with
`soil_moisture_tension`, but only 8 reported in October 2026, and it carries no
SNOTEL soil moisture:

| station | network | sensors (depth) | soil record from | also in NRCS? |
|---|---|---|---|---|
| C3MYD Maywood | CW3E | 6: 5, 10, 15, 20, 50, 100 cm (HydraProbe) | 2025-08 | no |
| C3STP Stampede | CW3E | 6: 5, 10, 15, 20, 50, 100 cm (HydraProbe) | 2026-10 | no |
| PVC54 Paradise | NWAC | 2, depth and sensor not given | 2023-07 | no |
| SUN64 Sunrise Base | NWAC | 2, depth and sensor not given | 2020-11 | no |
| OHAW1 Ohanapecosh | RAWS | 1, depth not given | 2016-11 | no |
| HZCW1 Hozomeen | RAWS | 1, depth not given | 2025-09 | no |
| COOW1 Cook Farm Field D | SCAN | 1 of 5: the 5 cm sensor | 2017-08 | **yes**, all 5 depths |
| LDDW1 Lind #1 | SCAN | 1 of 5, reads 0.0 | 2017-08 | **yes** |

The Cook Farm series matches AWDB's 2 in sensor exactly (checked 2026-11-09), and AWDB
has every hour where Synoptic's copy misses about half. So Synoptic adds two
well-documented but young CW3E profiles and four sensors of unknown depth. Not
usable: HADS (9) and CRN (2) soil sensors that stopped in 2006–2008 at stations still
"active" (station status is not sensor status: see
`SENSOR_VARIABLES.<var>.<set>.PERIOD_OF_RECORD`), RAWS tension sensors (none reporting),
and CW3E Headworks and YLine (listed, but access denied).

`gaia stage synoptic -v soil_moisture` works, but keeps only the first sensor set and
no depth. Branch `synoptic-soil-sensors` (commit `5a5163a`, not merged) kept every set,
named by depth from `sensorvars=1` positions (`soil_moisture_10cm`) or `_set_n`, and
converted to m3 m-3. Revisit it if the CW3E sites, or Synoptic's denser coverage
elsewhere (Oregon's CPC and SKYPN), are needed. Since main now stages source names, a
lighter version would keep `_set_n` and put the depth in an attribute.

Synoptic history: free accounts get about a year. Older requests are refused with
"does not have access to the requested history", and a request straddling the limit
is silently truncated to it. Our token was granted longer access through the academic
Open Access program on 2026-11-02.

## Other sources

- **NOAA USCRN**: [hourly02 files](https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/),
  one per station and year, with soil moisture (m3 m-3) and temperature at 5, 10, 20, 50
  and 100 cm. WA: Darrington 21 NNE, Quinault 4 NE, Spokane 17 SSW. OR: Coos Bay, Corvallis,
  John Day, Riley. AK: 24 stations (not every Alaska station has soil sensors: check
  before relying on one). Few stations, but high-quality, triple-redundant sensors. The
  natural next station source.
- **USBR AgriMet**: Pacific Northwest agricultural weather stations;
  [station parameters](https://www.usbr.gov/pn/agrimet/aginfo/station_params.html) list
  soil moisture channels (`XSM<n>`, percent at <n> in).
- **WSU AgWeatherNet**: eastern Washington agricultural stations, soil moisture at
  8 in where installed ([API](https://weather.wsu.edu/api/v1/network/stations)).
- **AmeriFlux**: flux towers with `SWC` (percent by volume); depths in site metadata.
- **Cook Farm**: USDA ARS R.J. Cook Agronomy Farm sensor network, 30–150 cm,
  [doi:10.15482/USDA.ADC/1349683](https://doi.org/10.15482/USDA.ADC/1349683).
- **ISMN** (International Soil Moisture Network): harmonized, quality-flagged in-situ
  data from many networks, worldwide; useful beyond the US, but it duplicates NRCS
  and USCRN here.
- **Gridded**: ERA5 is staged. SMAP and NISAR need NASA Earthdata (see
  [roadmap/04-datasets.md](roadmap/04-datasets.md)); soil-moisture-fusion processes
  NISAR SME2 itself.

## Instruments and quantities

- **Volumetric water content** (`SMS`, Synoptic `soil_moisture`, percent; staged as
  m3 m-3): water volume per soil volume, measured indirectly by dielectric sensors.
  Water's permittivity (~80) is far higher than minerals' (~4) or air's (1), and a
  calibration curve, often generic, turns the sensor's response into water content.
  Stevens HydraProbe (coaxial impedance) at SNOTEL, SCAN and CW3E; Campbell CS616/CS650
  (reflectometers) at some older CW3E and RAWS stations. Different sensors,
  calibrations and soils make absolute values hard to compare between stations:
  saturation reads ~40 % in a forest loam, ~25 % in a sand.
- **Matric potential** (Synoptic `soil_moisture_tension`, centibars = kPa): how
  tightly the soil holds its water, from granular matrix sensors such as the Irrometer
  Watermark. 0 is saturated, 10–30 kPa field capacity, 200 kPa the sensor's dry limit.
  Turning it into water content needs the soil's water retention curve, a model.
