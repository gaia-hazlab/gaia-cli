# Station data

`gaia stage synoptic` keeps **every report**, not a dense hourly grid. The
`stations/synoptic` group is a CF Discrete Sampling Geometries `timeSeries` stored as
a [contiguous ragged array](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#_contiguous_ragged_array_representation).
`gaia_cli.stations` turns it into a table or a regular grid when you need one. See
[roadmap/03-storage.md](roadmap/03-storage.md#a-one-zarr-store-with-cf-ragged-arrays-for-observations-recommended-to-start)
for why.

## Layout

| dimension | variables |
|---|---|
| `geometry` (one per station) | `station` (id, `cf_role: timeseries_id`), `name`, `elevation` (m), `network_id` (Synoptic `MNET_ID`), `precip_source`, `row_size`, and the point geometry (xvec CF-encoded) |
| `obs` (one per report) | `time` (UTC), `precipitation`, and any other requested variables |

The group has `featureType: timeSeries`. `row_size` has `sample_dimension: obs`: the
first `row_size[0]` reports belong to the first station, the next `row_size[1]` to the
second, and so on. Each station's reports are in time order.

- `precipitation` is the running total since the station's first report in the window
  (mm), as in [precip-timing.md](precip-timing.md). It is Synoptic's derived
  `precip_accumulated_set_1d`. If the account cannot use derived precipitation, it is
  the raw `precip_accum` with gauge resets removed. `precip_source` names the reported
  variable it was derived from (`precip_accum`, `precip_accum_one_hour`,
  `precip_accum_since_local_midnight`, …).
- `--variables` adds other Synoptic variables, such as `air_temp` or `wind_speed`,
  under their Synoptic names. Each takes the station's first sensor set (`_set_1`), as
  a number, with Synoptic's `units`. A report missing a variable has NaN for it.
- A report with none of the variables is dropped, and so is a station without any
  such report.
- QC flags are not stored yet.

```bash
gaia stage synoptic -i $AOI -s 2025-12-08 -e 2025-12-10 -v precip_accum air_temp -o skagit-stations
```

or, in a stage file:

```yaml
datasets:
  synoptic: {variables: [precip_accum, air_temp]}
```

## Reading

Both helpers take the group as read from the store (CF-encoded) or as returned by
`synoptic.load()`:

```python
import xarray as xr
from gaia_cli import stations

ds = xr.open_dataset("skagit-stations/skagit-stations.zarr", engine="zarr",
                     group="stations/synoptic")

df = stations.to_dataframe(ds)   # one row per report: station, time, precipitation, …
hourly = stations.to_grid(ds)    # dense (geometry, time) grid, hourly by default
totals = stations.to_grid(ds, freq="1D", how="sum", variables=["precipitation"])
```

`to_grid` puts the reports on a regular grid:

- **`how="last"`** (default): each grid time *t* holds the latest report at or before
  *t* that is less than one `freq` old, taken separately for each variable. A station
  that reports temperature every 5 minutes and precipitation hourly keeps both. Times
  without such a report are NaN. `<variable>_report_time` records when each value was
  reported. With `freq="1h"` and only `precipitation`, this is exactly the hourly array
  that `gaia stage synoptic` used to write.
- **`how="sum"`**, for running totals: the amount in each interval `(t - freq, t]`,
  the difference of the running total at both ends. `time_bounds` gives the
  intervals. It is NaN where either end has no report.

Station variables become coordinates along `geometry`, which keeps its xvec index, so
`hourly.max("time").xvec.to_geopandas()` works. Nothing derived is written to the
store.
