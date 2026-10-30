# Notes on Hugging Face Storage Buckets for Gaia

[Storage Buckets](https://huggingface.co/docs/hub/storage-buckets) are mutable, non-versioned object storage on the Hugging Face Hub, owned here by the [gaia-hazlab](https://huggingface.co/gaia-hazlab) org.

Buckets:
1. `gaia-hazlab/gaia-scratch`: scratch for Gaia. Public read. Emptied on the 1st of each month by [wipe-hf-scratch.yml](../.github/workflows/wipe-hf-scratch.yml).
2. `gaia-hazlab/gaia-public`: public test data (a USGS 3DEP COG).

There are no lifecycle rules, ACLs, bucket policies or versioning; the only access setting is public or private, and deletes are permanent. Manage buckets with the `hf` CLI (`pixi global install huggingface_hub`, upgrade with `pixi global update huggingface_hub`):

```bash
hf buckets list gaia-hazlab                      # buckets in the org
hf buckets list gaia-hazlab/gaia-scratch -h -R   # files in a bucket
hf buckets create gaia-hazlab/<name>             # add --private for a private bucket
```

## Writing with gaia: the S3 gateway

`gaia` writes through the [S3-compatible gateway](https://huggingface.co/docs/hub/storage-buckets-s3), so HF is just another AWS profile and no code is HF-specific. Put the org in the endpoint and use the bare bucket name in `s3://` URLs:

```ini
# ~/.aws/config
[profile hf]
region = us-east-1
endpoint_url = https://s3.hf.co/gaia-hazlab
request_checksum_calculation = when_required
response_checksum_validation = when_required
s3 =
    addressing_style = path

# ~/.aws/credentials
[hf]
aws_access_key_id = HFAK...
aws_secret_access_key = ...
```

```bash
AWS_PROFILE=hf aws s3 ls s3://gaia-scratch/
AWS_PROFILE=hf pixi run gaia stage prism ... -o s3://gaia-scratch/skagit-prism
```

S3 credentials are generated from a Hub access token: at [settings/tokens](https://huggingface.co/settings/tokens), make a fine-grained token with write access to `gaia-hazlab` repos, then choose **Generate S3 credentials** from its menu. The secret is shown once. The credentials have the token's permissions.

Compared with B2 (see [backblaze.md](backblaze.md)), the gateway honors conditional PUTs (`If-None-Match`), lists with ListObjectsV2 only, and rejects keys with `//` or a leading `/`.

## Reading anonymously

The gateway always requires a signature, even for a public bucket. Anonymous readers use HTTPS instead:

```
https://huggingface.co/buckets/gaia-hazlab/<bucket>/resolve/<key>     # 302 to a signed CDN URL, valid 1 hour
https://huggingface.co/api/buckets/gaia-hazlab/<bucket>/tree/<prefix>  # listing (JSON)
```

## GitHub Actions

- **Staging** ([stage.yml](../.github/workflows/stage.yml), [stage-event.yml](../.github/workflows/stage-event.yml)): pick `storage: huggingface`. It reads the `huggingface` GitHub Environment:
  - secrets `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` (the `HFAK...` S3 credentials)
  - vars `AWS_ENDPOINT_URL=https://s3.hf.co/gaia-hazlab`, `AWS_DEFAULT_REGION=us-east-1`, `S3_BUCKET=gaia-scratch`
- **Monthly wipe** ([wipe-hf-scratch.yml](../.github/workflows/wipe-hf-scratch.yml)): authenticates with a [Trusted Publisher](https://huggingface.co/docs/hub/trusted-publishers) instead of a secret. On the bucket's **Settings → Trusted Publishers**, add GitHub Actions with `repository=gaia-hazlab/gaia-cli`, `branch=main`, `workflow=wipe-hf-scratch.yml`. A manual run defaults to a dry run.
