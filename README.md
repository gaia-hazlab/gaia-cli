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
Usage: gaia stage COMMAND

Stage datasets clipped to an AOI and time range.

╭─ Commands ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ hrrr      Stage HRRR analysis data clipped to an AOI. The HRRR precipitation rate (kg m-2 s-1) is converted to hourly accumulation (mm).  Requires AWS_PROFILE environment    │
│           variable if writing to S3.                                                                                                                                          │
│ prism     Stage PRISM daily precipitation clipped to an AOI.                                                                                                                  │
│ synoptic  Stage in-situ precipitation from the Synoptic API. Precipitation timeseries is resampled to construct a uniform hourly DataArray. Requires SYNOPTIC_TOKEN           │
│           environment variable (https://docs.synopticdata.com). Requires AWS_PROFILE environment variable if writing to S3.                                                   │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```


```
Usage: gaia stage synoptic [ARGS]

Stage in-situ precipitation from the Synoptic API. Precipitation timeseries is resampled to construct a uniform hourly DataArray. Requires SYNOPTIC_TOKEN environment variable
(https://docs.synopticdata.com). Requires AWS_PROFILE environment variable if writing to S3.

╭─ Parameters ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ INPUT --input -i    Path to a vector file (e.g. GeoJSON, Shapefile) that contains the area-of-interest polygon. The first geometry in the file is used.                             │
│ START --start -s    Start date/time for the timeseries (parseable by pandas.Timestamp).                                                                                             │
│ END --end -e        End date/time for the timeseries (parseable by pandas.Timestamp).                                                                                               │
│ OUTPUT --output -o  Path to write the resulting Zarr store.                                                                                                                         │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

### Examples


Example workflow (Skagit Basin 2025-12 Atmospheric River Event)
```bash
AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
gaia stage --config conus-precip -i $AOI -s 2025-12-01 -e 2025-12-20 -o /tmp/skagit-test.zarr
```

`conus-precip` is a built-in preset: Synoptic stations, PRISM and HRRR precipitation. Every dataset declares its coverage, so datasets that the AOI or window misses (for example PRISM and HRRR in Alaska) are skipped with a warning. A dataset that fails doesn't stop the others; rerun with `--skip-existing` to retry only the missing groups.

`stage` keeps every raster on its native grid and time step. To compare products, align them onto one grid (default: the coarsest, HRRR 3 km) and one time step (default: the coarsest, PRISM days) in a new store. See [docs/reprojecting.md](docs/reprojecting.md):
```bash
gaia align -i /tmp/skagit-test.zarr -o /tmp/skagit-test-aligned.zarr --aoi $AOI
```

Example stage file: one YAML file gives the AOIs, time windows and options for several datasets (see [docs/stage-file.md](docs/stage-file.md)). `--dry-run` validates it and prints the resolved configuration without fetching data:
```bash
gaia stage --config examples/skagit-2025.yaml --dry-run
gaia stage --config examples/skagit-2025.yaml -o /tmp/skagit-config.zarr
gaia stage --config examples/skagit-2025.yaml -o /tmp/skagit-config.zarr --only prism,hrrr
```

Example static landcover staging (NLCD 2021)
```bash
AOI=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json
gaia stage nlcd -i $AOI -o s3://gaia/scratch/skagit-nlcd-2021.zarr --year 2021
```

### Output to S3

Zarr outputs are well-suited to S3-compatible object stores (e.g. AWS S3, MinIO, etc.). The examples below write to either a University of Washington object store called [Kopah](https://hyak.uw.edu/docs/systems/kopah/overview/) or [backblaze B2](https://www.backblaze.com/docs). We output to a `scratch` prefix that has a 30-day lifecycle policy (any uploaded file is deleted after 30 days). You must have credentials to write to this bucket, credentials can be configured by setting environment variables `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` or by using the AWS config files and settings AWS_PROFILE to the profile name. For example, to use the `kopah` profile, set `AWS_PROFILE=kopah` in your environment:

`~/.aws/config`:
```
[profile kopah]
endpoint_url = https://s3.kopah.uw.edu
request_checksum_calculation = when_required
response_checksum_validation = when_required
s3 =
    payload_signing_enabled = false
    checksum_algorithm = crc32

[profile b2]
region = us-west-004
endpoint_url = https://s3.us-west-004.backblazeb2.com
```

`~/.aws/credentials`:
```
[kopah]
aws_access_key_id = YOUR_ACCESS_KEY
aws_secret_access_key = YOUR_SECRET_KEY

[b2]
aws_access_key_id = YOUR_ACCESS_KEY
aws_secret_access_key = YOUR_SECRET_KEY
```

## Analysis

See [precip-example.ipynb](./precip-example.ipynb)

## Provenance

Every store records where its data came from, using STAC field names. The
root group says when the store was created, by which gaia-cli commit
(`gaia:cli-version`) and with what command. Each product group (or the root of
a single-product store) carries its upstream `title`, `license`, `providers`,
`version` and `links`:

```python
dt = xr.open_datatree("skagit_all.zarr", engine="zarr")
dt.attrs["gaia:cli-version"]
dt["rasters/hrrr"].attrs["license"]  # 'CC-BY-4.0'
```

Products have different licenses, so check each group's `license` and `links`
before sharing a store. See [docs/provenance-metadata.md](docs/provenance-metadata.md).

Each store also gets a STAC Item next to it (`out.zarr` → `out.json`).
`gaia catalog rebuild -p s3://bucket/prefix` links every Item under a prefix
into a `catalog.json` for browsing in stac-browser; see
[docs/stac-catalog.md](docs/stac-catalog.md).


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

```bash
gh workflow run stage.yml \
  --field workflow_name=skagit_test \
  --field subcommand=all \
  --field aoi=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json \
  --field start_date=2025-12-01 \
  --field end_date=2025-12-02 \
  --field copy_to_s3=false
```


```bash
gh workflow run stage.yml \
  --field workflow_name=skagit \
  --field subcommand=all \
  --field aoi=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json \
  --field start_date=2025-12-01 \
  --field end_date=2025-12-30 \
  --field copy_to_s3=true
```
