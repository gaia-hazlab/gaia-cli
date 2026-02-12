# Data Staging CLI for GAIA Analysis

This is currently just a prototype. It illustrates a basic CLI to get preciptation time series for a given AOI and do a basic comparison against a raster dataset.

The basic approach is to:

1. Get heteregenous data from many stations from an API.
1. Normalize / munge to a uniform Xarray Vector Data Cube with [XVec](https://xvec.readthedocs.io)
1. Save to [Zarr](https://zarr.readthedocs.io) in a temporary S3 bucket
1. Now any group can easily load the timeseries, compare to rasters, use different libraries if they want, etc.

## Usage

You need [Pixi](https://pixi.prefix.dev/latest/installation/) to run this CLI currently:

```bash
gh repo clone gaia-hazlab/gaia-cli
pixi shell
gaia stage --help
```

```
Usage: gaia stage [ARGS]

Stage data by reading an AOI, fetching precipitation data, and writing Zarr. Precipitation timeseries is resampled to construct a uniform hourly DataArray.
Requires SYNOPTIC_TOKEN environment variable (https://docs.synopticdata.com). Requires AWS_PROFILE environment variable if writing to S3.

╭─ Parameters ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ INPUT --input -i    Path to a vector file (e.g. GeoJSON, Shapefile) that contains the area-of-interest polygon. The first geometry in the file is used.       │
│ START --start -s    Start date/time for the timeseries (parseable by pandas.Timestamp).                                                                       │
│ END --end -e        End date/time for the timeseries (parseable by pandas.Timestamp).                                                                         │
│ OUTPUT --output -o  Path to write the resulting Zarr store.                                                                                                   │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

Example workflow (Skagit Basin 2025-12 Atmospheric River Event)
```bash
AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
gaia stage -i $AOI -s 2025-12-01 -e 2025-12-15 -o s3://cresst/test.zarr
```

## Analysis

See [precip-example.ipynb](./precip-example.ipynb)
