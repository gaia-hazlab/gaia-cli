# Developer notes

Notes for people working on gaia-cli. User documentation is in [../docs/](../docs/),
published at <https://gaia-hazlab.github.io/gaia-cli/>.

| note | about |
|---|---|
| [roadmap/](roadmap/README.md) | plan: more datasets, regions and workflows |
| [adding-a-dataset.md](adding-a-dataset.md) | the module contract for a new dataset |
| [stage-file.md](stage-file.md) | full stage file grammar: AOIs, windows, coverage, output folder |
| [provenance-metadata.md](provenance-metadata.md) | provenance attributes in each store |
| [stac-catalog.md](stac-catalog.md) | STAC Items and `gaia catalog rebuild` |
| [prism.md](prism.md), [hrrr.md](hrrr.md), [stations.md](stations.md), [precip-timing.md](precip-timing.md), [dem.md](dem.md), [soil.md](soil.md), [burn-severity.md](burn-severity.md), [seismic.md](seismic.md) | per-dataset details: sources, grids, time conventions, caveats |
| [backblaze.md](backblaze.md) | Backblaze B2 setup |
| [huggingface.md](huggingface.md) | Hugging Face bucket setup (S3 gateway, monthly scratch wipe) |
| [releasing.md](releasing.md) | versions from git tags; cutting a release |

## Object store credentials

Zarr outputs are well-suited to S3-compatible object stores (e.g. AWS S3, MinIO, etc.). The examples below write to either a University of Washington object store called [Kopah](https://hyak.uw.edu/docs/systems/kopah/overview/), [backblaze B2](https://www.backblaze.com/docs), or a [Hugging Face bucket](huggingface.md). We output to a `scratch` prefix that has a 30-day lifecycle policy (any uploaded file is deleted after 30 days). You must have credentials to write to this bucket, credentials can be configured by setting environment variables `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` or by using the AWS config files and settings AWS_PROFILE to the profile name. For example, to use the `kopah` profile, set `AWS_PROFILE=kopah` in your environment:

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

[profile hf]
region = us-east-1
endpoint_url = https://s3.hf.co/gaia-hazlab
request_checksum_calculation = when_required
response_checksum_validation = when_required
s3 =
    addressing_style = path
```

`~/.aws/credentials`:
```
[kopah]
aws_access_key_id = YOUR_ACCESS_KEY
aws_secret_access_key = YOUR_SECRET_KEY

[b2]
aws_access_key_id = YOUR_ACCESS_KEY
aws_secret_access_key = YOUR_SECRET_KEY

[hf]
aws_access_key_id = HFAK...
aws_secret_access_key = YOUR_SECRET_KEY
```


## Staging with GitHub Actions

Run the data stager as a github action, saving the output folder (store, STAC Item and `stage.yaml`) as a build artifact. `subcommand` is a dataset (`synoptic`, `prism`, `hrrr`) or the `conus-precip` preset.

```bash
gh workflow run stage.yml \
  --field workflow_name=skagit_test \
  --field subcommand=conus-precip \
  --field aoi=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json \
  --field start_date=2025-12-01 \
  --field end_date=2025-12-02 \
  --field copy_to_s3=false
```


```bash
gh workflow run stage.yml \
  --field workflow_name=skagit \
  --field subcommand=conus-precip \
  --field aoi=https://raw.githubusercontent.com/DSHydro/skagit-met/refs/heads/main/data/GIS/SkagitBoundary.json \
  --field start_date=2025-12-01 \
  --field end_date=2025-12-30 \
  --field copy_to_s3=true
```


`stage-event.yml` stages a whole stage file: a path in the repository or an http(s)
URL (whose AOIs are then URLs or inline). It can then align the store and write
Landlab inputs from it. The stores go to `s3://gaia/scratch/github-workflows/<name>`
and `<name>-aligned`; every output folder is also kept as an artifact for 7 days.

```bash
gh workflow run stage-event.yml \
  --field config=examples/pioneer.yaml \
  --field name=pioneer \
  --field align="--crs utm --res 10 --aoi examples/pioneer-aoi.geojson -m dem/usgs10m=cubic -m soil/solus100=bilinear -m rasters/nlcd=nearest -m rasters/burn-severity=nearest" \
  --field landlab_inputs=true \
  --field copy_to_s3=true
```
