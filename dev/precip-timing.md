# Synoptic precipitation timing

How station precipitation is put on an hourly time axis, what was wrong with the first
version, and why Synoptic's precipitation service does not solve the problem.

`gaia stage synoptic` now stores every report ([stations.md](stations.md)), and
`gaia_cli.stations.to_grid(ds)` builds the hourly axis described here when you read
the store. Its `precipitation_report_time` is the `report_time` below.

## Summary

- Synoptic stations do not all report on the hour. In the Skagit test case (19 stations,
  2025-12-08 to 2025-12-11) four report hourly at :12, :23 or :28. Several RAWS stations
  report on the hour but skip hours (gaps of 4 h and 11 h).
- **Old behavior (fixed):** `resample("h").bfill(limit=1)` labeled each hour with the
  *next* report. Hourly values therefore included precipitation that fell after their
  timestamp, by up to 48 minutes. The hour just before a gap ended received the
  accumulation reported after the gap.
- **New behavior:** each hour takes the accumulation at the *latest report at or before*
  it, as long as that report is less than an hour old. Otherwise the hour is NaN. A new
  `report_time` coordinate `(station, time)` records which report each value came from.
- **Effect:** stations that report on the hour have identical values. For the four
  off-hour stations, two-day totals over 2025-12-08T12 → 12-10T12 drop by 3.6–5.8 mm
  (about 4–5 %), because the old values included up to 48 minutes of rain after the
  window end.
- **The precipitation service (`/stations/precip?pmode=intervals`) is not a fix.**
  Its intervals are bounded by actual reports, which can start before the requested
  interval. Intervals that span a reporting gap are dropped, and with them the rain
  that fell during the gap.

## How Synoptic reports precipitation

`gaia stage synoptic` calls `/stations/timeseries` with `vars=precip_accum&precip=1`.
Synoptic returns its derived `precip_accumulated_set_1d`, a running total since the
start of the request, at each report time. Report times vary by network:

| cadence | example stations (Skagit) | report minutes |
|---|---|---|
| hourly, on the hour | RAIW1, THBW1, SWSW1, BVPW1, HZOW1, MRTW1, DCKW1, GODW1, UBDW1 | :00 (RAWS here skip some hours) |
| hourly, off the hour | FIFW1, GHFW1 (:28), MBMW1 (:12), HZCW1 (:23) | one fixed minute |
| sub-hourly | UBKW1, DARW1, MARW1, MMCW1, COOP… | every 5–15 min |

An hourly report at :28 cannot tell how much of an hour's rain fell before :00. **No
resampling can produce clock-hour totals for these stations.** The best we can do is
label honestly and record when each value was observed.

## Old vs new labeling

For a station reporting at :28:

| hourly label | old (`bfill`) | new (as-of) |
|---|---|---|
| 01:00 | accumulation at 01:28 (28 min *after* the label) | accumulation at 00:28 (32 min *before*) |
| 02:00 | accumulation at 02:28 | accumulation at 01:28 |

Neither label equals the clock hour. The new one never includes rain that had not
fallen by the timestamp, which is what code like
`precip.sel(time=t1) - precip.sel(time=t0)` assumes (see `precip-example.ipynb`). The
offset is now visible in `report_time` rather than hidden.

For gaps, consider a station that reports at 01:00 and then next at 04:00:

| hourly label | old | new |
|---|---|---|
| 02:00 | NaN | NaN |
| 03:00 | accumulation at 04:00 (all gap rain assigned before it fell) | NaN |
| 04:00 | accumulation at 04:00 | accumulation at 04:00 |

### Measured on the Skagit test case

The table shows totals for the two PRISM days 2025-12-08T12:00 → 2025-12-10T12:00, in mm.
Stations not listed are unchanged.

| station | report minute | lag (new) | total old | total new | largest change in an hourly increment |
|---|---|---|---|---|---|
| FIFW1 | :28 | 32 min | 127.51 | 121.67 | 4.07 |
| GHFW1 | :28 | 32 min | 92.45 | 88.90 | 3.54 |
| MBMW1 | :12 | 48 min | 118.11 | 112.77 | 4.57 |
| HZCW1 | :23 | 37 min | 65.79 | 61.72 | 2.54 |

Hourly maxima at these stations are 6–10 mm. A 2.5–4.6 mm change in individual
hourly increments therefore matters for any hourly comparison with HRRR. For hourly
work, use only values where `report_time == time`, or account for the lag.

## Why not Synoptic's precipitation service?

The service returns one total per requested interval, with `count`, `first_report`,
`last_report` and `report_type`. landlab_data_prep uses it. Requesting 1-hour
intervals (`interval=1`, `interval_window=0`) for the same stations and dates shows:

1. **Intervals follow the reports, not the clock.** For FIFW1, the first interval of a
   request starting at 00:00 is `first_report 23:28 → last_report 00:28`. That is
   shifted 32 minutes earlier, and it starts before the requested start. The timing
   problem is the same as above, just in the other direction.
2. **Intervals that span a reporting gap are dropped.** RAIW1, THBW1 and HZOW1 have 4 h
   and 11 h gaps. The service returns 45 contiguous rows instead of 61. The rain that
   fell across the gaps (7.6, 22.9 and 25.4 mm) is missing from the sum of intervals,
   while the running total keeps it. Two-day totals from summed intervals are 23–36 %
   low at those stations.
3. **`interval` in the response is an index, not a length.** A single 6-hour request
   returns `"interval": 1`.

The service is useful for its completeness metadata. The running total plus
`report_time` carries the same information: an interval total `cum(t1) - cum(t0)`
is exact when both `report_time` values equal their labels, and you know the offset
when they don't.

### Notes for landlab_data_prep

landlab requests one interval at a time and accepts a total only if
`float(record["interval"]) == interval_hours` and `start <= first_report`. From the
behavior above:

- With `interval_hours` other than 1, every record has `interval == 1`, so every
  record would be rejected.
- Off-hour stations always have `first_report` before the requested start, so they are
  always rejected.
- Gap-spanning intervals are never returned.

All three are worth checking with @Almehedi06.

## Is the staged format enough for Landlab?

Yes. One staged format serves both uses:

- **gaia comparisons** use the running total on the hourly axis, plus `report_time`.
- **Landlab forcing** needs N-hour totals per interval, assigned to a grid cell from
  the nearest station. An N-hour total is `cum(t) - cum(t - N h)`, valid when both
  ends have a `report_time`. landlab's stricter rule (both reports inside the interval)
  is `report_time == time` at both ends. Nearest-station gridding and the ASCII files
  are model-specific, so they belong in `gaia convert landlab`
  ([issue #5](https://github.com/gaia-hazlab/gaia-cli/issues/5)), not in a second
  staged format.

## Reproducing

This comparison predates the ragged layout: `construct_uniform_datarray` has since
been replaced by `synoptic.construct_ragged` followed by `stations.to_grid`, which
gives the same hourly values.


The comparison used `gaia_cli.synoptic.get_station_metadata` and
`get_station_timeseries` for the Skagit AOI, 2025-12-08 → 2025-12-11. It then ran
`/stations/precip` with `pmode=intervals&interval=1&interval_window=0&units=precip|mm`
for the same station IDs, and compared the old and new `construct_uniform_datarray`
outputs. Synoptic responses for past dates can change as QC is applied, so expect
small differences on a rerun.
