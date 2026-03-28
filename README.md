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
gaia stage -i $AOI -s 2025-12-01 -e 2025-12-20 -o s3://cresst/scratch/skagit-test.zarr
```

## Analysis

See [precip-example.ipynb](./precip-example.ipynb)


## Design ideas:

- We want a single entrypoint for people to open a catalog for a given AOI and time range. We put everything into a [xr.DataTree](https://xarray.pydata.org/en/stable/data-structures.html#data-tree) structure to preserve the original datasets dimensions and coordinates, and have a second step "e.g. create datacube" that merges into a common grid and CRS.

- We'll use Zarr to serialize to disk. For now, we'll use default chunking and compression.

- Various datasets use different CRS & grids (e.g. PRISM uses lon/lat NAD83, HRRR uses Lambert Conformal Conic NAD83, COP30DEM is WGS84 (G1150)). For now we reproject everything to  WGS84 (G1150) for simplicity.

- Coordinate names vary (longitude vs lon vs x, latitude vs lat vs y). We convert to 'x', 'y', 'time' before saving. If multiple 'time' dimensions exist (e.g. hourly station sampling vs daily reanalysis or satellite products, we'll rename to the dataset name + '_time' (e.g. 'prism_time') to avoid conflicts.

- Standardize variable names and units:
    - Precipitation: 'precipitation' in mm
    - Temperature: 'temperature' in degC
    - Wind Speed: 'wind_speed' in m/s
    - Soil Moisture: 'soil_moisture' in m3/m3


## GitHub Actions Workflow

Run the data stager as a github action, saving resulting zarr as a build artifact.

```
gh workflow run stage.yml --field aoi=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json --field start=2025-12-01 --field end=2025-12-20
```
