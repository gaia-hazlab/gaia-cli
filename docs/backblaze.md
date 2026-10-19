# Notes on configuring Backblaze B2 for Gaia

UW IT Kopah Object store is running an older version of Ceph that does not support logging and it is on average down 3-4 days per month wich is a bit frustrating. We'll try to use Backblaze B2 for our scratch object store as a replacement. It is S3-compatible, the first 10 GB of storage is always free, and paid storage is $6.95/TB/month (compared to Kopah's $5/TB/month).

Docs: https://www.backblaze.com/docs

Buckets:
1. gaia-scratch : scratch bucket for Gaia. 30-day lifecycle policy (delete after 30 days). Public read. CORS policy set to allow [STACBrowser](https://github.com/radiantearth/stac-browser).
2. gaia-data : long-term storage for Gaia. No lifecycle policy. Credentials required to read/write.
3. gaia-logging : logging bucket for Gaia. 30-day lifecycle policy.

## Free tier vs. paid

The free tier is the same B2 account with a payment method missing. Most features work; a few are gated.

| Feature | Free account? | Notes |
|---|---|---|
| Private buckets, application keys, S3 API | Yes | |
| Lifecycle rules | Yes | |
| CORS rules | Yes | |
| Public buckets | **Needs verification** | Before your *first* public bucket you must verify your email **and** either have a payment history or make a small card payment (credited to your balance). |
| Bucket access logs | **No** | "Available to all paying customers." Docs also say it may still need to be enabled by Backblaze sales; confirm before relying on it. |
| Storage beyond 10 GB | No | Uploads past the cap fail until you add a card. |

[Pricing](https://www.backblaze.com/cloud-storage/pricing), as of Oct 2026:
- Egress is free up to 3x your average monthly storage, then $0.01/GB.
- Class A/B/C API calls are free; Class D is $0.004 per 10k after 2,500/day free.
- There is no minimum storage duration and no minimum file size fee.

Watch out for these:
- **Bucket names are global across all B2 customers** (6–63 chars, letters/numbers/hyphens). `gaia-scratch` may already be taken; if so, use something like `gaia-hazlab-scratch`.
- **Buckets keep every version by default.** This differs from AWS S3, where new buckets have versioning off and an overwrite replaces the object. In B2, buckets are always versioned and the console default is "Keep All Versions", forever. Overwriting or deleting a file by name only hides it, and the old version stays stored ([S3 API bucket versions](https://www.backblaze.com/docs/cloud-storage-s3-compatible-api-bucket-versions), [lifecycle rules](https://www.backblaze.com/docs/cloud-storage-lifecycle-rules)). The docs don't explicitly say hidden versions are billed, but they are stored data, so assume they are. Re-writing Zarr stores (many small chunk overwrites) could quietly pile up storage. Every bucket below therefore gets at least a "delete hidden versions after 1 day" rule (equivalent to the console's "Keep Only Last Version"), including `gaia-data`.
- The region (e.g. `us-west-004`) is fixed when the account is created. Access logs must go to a bucket in the same region.

## S3 compatibility gaps

"S3-compatible" stores don't all implement the same features. These results come from testing against `gaia-scratch` (Oct 2026; zarr 3.4.1, obstore, icechunk 2.3.0). Neither feature can be turned on with a bucket setting.

| Feature | AWS S3 | Kopah (Ceph) | B2 |
|---|---|---|---|
| Anonymous LIST on a public bucket | Yes (bucket policy) | Yes (`public/` prefix) | **No.** 403 `Unauthenticated requests are not allowed for this api` |
| Conditional create (`If-None-Match: *`) | Yes | Yes | **No.** 501 `NotImplemented` |

How `gaia` copes:
- **Anonymous reads:** zarr v3 stores written without consolidated metadata must be *listed* to discover their arrays. That fails anonymously on B2. `gaia` writes consolidated metadata (`consolidated=True`, inline in the root `zarr.json`), so readers need only GETs. zarr warns that consolidated metadata isn't in the v3 spec yet; `io.py` silences that warning. Non-Python zarr readers may not understand it.
- **Writes:** zarr creates new keys with a conditional PUT. `io._UnconditionalObjectStore` swaps that for check-then-write, which is safe because each output has one writer and `check_output` refuses existing paths.

**Icechunk** relies on conditional writes for safe commits:
- **Default config:** `Repository.create` fails with the same 501.
- **With conditional writes disabled:** create, commits, and anonymous reads all work, but only if you also relax the AWS SDK's checksum handling. Otherwise writes fail with `body checksum mismatch`. Anonymous reads needed no listing.

```python
config = icechunk.RepositoryConfig(
    storage=icechunk.StorageSettings(
        unsafe_use_conditional_create=False, unsafe_use_conditional_update=False
    )
)
# and run with:
#   AWS_REQUEST_CHECKSUM_CALCULATION=when_required AWS_RESPONSE_CHECKSUM_VALIDATION=when_required
```

The `unsafe_` flags remove Icechunk's protection against concurrent commits, so two writers committing at once can lose data. Only use them for single-writer repos; for multi-writer Icechunk, use AWS S3 or Kopah.

## Configuring Backblaze B2 for Gaia

There are two APIs:
- The **B2 native API** (via the `b2` CLI) can create a bucket with type, lifecycle, and CORS in one command. It's the easiest option for admin tasks.
- The **S3-compatible API** (via `aws` CLI, `boto3`, `obstore`) is what `gaia` uses for reading and writing data. Bucket access logs can only be configured through this API.

### 1. Install the CLI and authorize

The `b2` CLI is on PyPI, not conda-forge, so run it ad hoc with `uvx`. The first run builds an environment in uv's cache (`uv cache dir`, usually `~/.cache/uv`). It is not deleted after running, so later runs reuse it. Remove it with `uv cache clean b2` ([uv tools docs](https://docs.astral.sh/uv/concepts/tools/)).

```bash
uvx b2 version
```

Make an admin key in the web console (**Application Keys → Add a New Application Key**, "All buckets", Read and Write). Its capabilities should include `writeBuckets`, `writeBucketLogging`, and `readBucketLogging`. The **master key cannot be used with the S3 API**, so make a separate key even if you only use the `b2` CLI for setup.

```bash
export B2_APPLICATION_KEY_ID=<keyID>
export B2_APPLICATION_KEY=<applicationKey>
uvx b2 account get   # confirms auth; shows the S3 endpoint/region
```

### 2. Create the buckets

```bash
# Scratch: public, delete 30 days after upload, CORS for STAC Browser
uvx b2 bucket create gaia-scratch allPublic \
  --lifecycle-rule '{"fileNamePrefix": "", "daysFromUploadingToHiding": 30, "daysFromHidingToDeleting": 1}' \
  --cors-rules "$(cat docs/cors-rules.json)"

# Long-term data: private, keep current versions forever, purge overwritten/deleted versions
uvx b2 bucket create gaia-data allPrivate \
  --lifecycle-rule '{"fileNamePrefix": "", "daysFromUploadingToHiding": null, "daysFromHidingToDeleting": 1}'

# Logs: private, delete 30 days after upload
uvx b2 bucket create gaia-logging allPrivate \
  --lifecycle-rule '{"fileNamePrefix": "", "daysFromUploadingToHiding": 30, "daysFromHidingToDeleting": 1}'
```

To change an existing bucket, use `uvx b2 bucket update <name> [allPublic|allPrivate] --lifecycle-rule ... --cors-rules ...`. **Each update replaces the whole lifecycle/CORS config**, so pass every rule you want to keep. Inspect with `uvx b2 bucket get <name>`.

### 3. Create scoped keys for day-to-day use

Don't use the admin key in `gaia`. Make a key restricted to one bucket:

```bash
uvx b2 key create --bucket gaia-scratch gaia-scratch-rw \
  listBuckets,listFiles,readFiles,writeFiles,deleteFiles
```

### 4. Configure an AWS profile

`gaia` resolves the endpoint from the AWS profile (see `src/gaia_cli/io.py`), so add a profile:

```ini
# ~/.aws/config
[profile b2]
region = us-west-004
endpoint_url = https://s3.us-west-004.backblazeb2.com

# ~/.aws/credentials
[b2]
aws_access_key_id = <keyID>
aws_secret_access_key = <applicationKey>
```

Then:

```bash
AWS_PROFILE=b2 aws s3 ls s3://gaia-scratch/
AWS_PROFILE=b2 pixi run gaia stage prism ... -o s3://gaia-scratch/skagit-prism
```

## Enable public read access

Set the bucket type to `allPublic` at creation (above), or change it later:

```bash
uvx b2 bucket update gaia-scratch allPublic
```

S3 API equivalent: `aws s3api put-bucket-acl --bucket gaia-scratch --acl public-read`. B2 only supports the canned `private` / `public-read` ACLs at the **bucket** level. Objects inherit the bucket ACL, and per-object ACLs return 403.

The first time you do this, you'll be prompted for email verification + a payment (see above). This is a one-time step; afterwards the account stays on the free tier unless you exceed it. Without a payment, `b2 bucket create ... allPublic` fails with `Account has no payment history ... (no_payment_history)`. In that case, create the bucket `allPrivate` and see [Serving a private bucket through Cloudflare](#serving-a-private-bucket-through-cloudflare).

Public objects are readable without credentials at any of:

```
https://s3.us-west-004.backblazeb2.com/gaia-scratch/<key>      # S3 path-style (use this for STAC hrefs)
https://gaia-scratch.s3.us-west-004.backblazeb2.com/<key>      # S3 virtual-host style
https://f004.backblazeb2.com/file/gaia-scratch/<key>           # B2 native download URL
```

## Adding a lifecycle policy

B2 lifecycle rules have three fields, applied per `fileNamePrefix`:
- `daysFromUploadingToHiding`: hide the current version N days after upload (the equivalent of S3 `Expiration`).
- `daysFromHidingToDeleting`: permanently delete hidden/overwritten versions N days later (S3 `NoncurrentVersionExpiration`).
- `daysFromStartingToCancelingUnfinishedLargeFiles`: clean up abandoned multipart uploads (S3 `AbortIncompleteMultipartUpload`).

So "delete after 30 days" is `daysFromUploadingToHiding: 30` + `daysFromHidingToDeleting: 1`, which means an object is actually gone after about 31 days. Lifecycle runs once a day, so timing is approximate.

The S3 API (`aws s3api put-bucket-lifecycle-configuration`) also works, with these restrictions:
- Each `Expiration.Days` rule must be paired with an `ExpiredObjectDeleteMarker: true` rule on the same prefix.
- No transitions, date-based rules, tag filters, or overlapping prefixes.
- Don't mix S3 and native/web-console edits on the same bucket, because that regenerates rule IDs.

## Setting CORS policy

Saved as [cors-rules.json](cors-rules.json) and passed with `--cors-rules` (native API). Native CORS rules use B2's own schema:

```json
[
  {
    "corsRuleName": "publicReadAnyOrigin",
    "allowedOrigins": ["*"],
    "allowedOperations": ["s3_get", "s3_head", "b2_download_file_by_name", "b2_download_file_by_id"],
    "allowedHeaders": ["*"],
    "exposeHeaders": ["ETag", "Content-Length", "Content-Range", "Content-Type", "Accept-Ranges", "Last-Modified"],
    "maxAgeSeconds": 3600
  }
]
```

Or via the S3 API with the AWS-style schema:

```
{
  "CORSRules": [
    {
      "ID": "public-read-any-origin",
      "AllowedOrigins": ["*"],
      "AllowedMethods": ["GET", "HEAD"],
      "AllowedHeaders": ["*"],
      "ExposeHeaders": ["ETag", "Content-Length", "Content-Range", "Content-Type", "Accept-Ranges", "Last-Modified"],
      "MaxAgeSeconds": 3600
    }
  ]
}
```

```bash
AWS_PROFILE=b2 aws s3api put-bucket-cors --bucket gaia-scratch --cors-configuration file://cors.json
```

These are two views of the same setting, so pick one. `MaxAgeSeconds` must be 0–86400, and B2 rejects empty origins.

## Bucket access logging

**Requires a paid account** (see table above). Once available, it's S3-API only. `b2-admin` is a second AWS profile set up like `b2`, but with the admin key:

```bash
cat > logging.json <<'EOF'
{"LoggingEnabled": {"TargetBucket": "gaia-logging", "TargetPrefix": "gaia-scratch/"}}
EOF
AWS_PROFILE=b2-admin aws s3api put-bucket-logging --bucket gaia-scratch --bucket-logging-status file://logging.json
AWS_PROFILE=b2-admin aws s3api get-bucket-logging --bucket gaia-scratch
```

- The target bucket must be in the same account and region, and must not have Object Lock.
- Never log a bucket into itself (that creates an infinite loop).
- The key needs `writeBucketLogging`/`readBucketLogging`.
- Logs use the AWS server-access-log format, are best-effort (may have gaps or duplicates), and count toward storage.
- To disable, send `{}` as the logging status.

## Deleting buckets

B2 won't delete a bucket that still contains files, including hidden versions. So empty it first, then delete it:

```bash
uvx b2 ls --recursive --versions b2://gaia-scratch-test   # preview what will be removed
uvx b2 rm --recursive --versions b2://gaia-scratch-test   # delete every version (irreversible)
uvx b2 bucket delete gaia-scratch-test
```

## Serving a private bucket through Cloudflare

This is a common alternative to a public bucket. The bucket stays private, and a Cloudflare Worker holds a read-only B2 key. The Worker signs each incoming request (AWS SigV4) and proxies it to B2's S3 endpoint. Clients see a plain public HTTPS URL. Backblaze maintains a ready-made Worker for this: [backblaze-b2-samples/cloudflare-b2](https://github.com/backblaze-b2-samples/cloudflare-b2).

Why bother:
- **No payment on B2 is needed**, because the bucket is never public.
- **Egress from B2 to Cloudflare is free** (Bandwidth Alliance), so downloads don't count against the 3x-storage egress allowance. Per third-party write-ups, this applies when users download through your own domain on Cloudflare; this wasn't confirmed in Backblaze's docs.
- **Access control and rate limiting live in the Worker.** You can turn off object listing, restrict to a prefix, etc.

Cloudflare free plan limits: 100,000 requests/day (1,000/min burst), 10 ms CPU per request (time waiting on B2 doesn't count), and 50 subrequests per request ([Workers limits](https://developers.cloudflare.com/workers/platform/limits/)). Range reads of a Zarr store or COG are one request each, so a busy STAC Browser/xarray session can use a lot of these, but it's fine for scratch-scale sharing.

### Setup

Requires Node.js ≥ 22 (`npx` is used for `wrangler`, Cloudflare's CLI).

1. **Create a read-only key scoped to the bucket.** Drop `listFiles` if you don't want listing.

   ```bash
   uvx b2 key create --bucket gaia-scratch gaia-scratch-cloudflare readFiles,listFiles
   ```

2. **Get the Worker code and log in to Cloudflare:**

   ```bash
   git clone https://github.com/backblaze-b2-samples/cloudflare-b2.git
   cd cloudflare-b2
   npm install
   npx wrangler login
   ```

3. **Configure `wrangler.toml` `[vars]`:**

   ```toml
   B2_APPLICATION_KEY_ID = "<keyID from step 1>"
   B2_ENDPOINT = "s3.us-west-004.backblazeb2.com"
   BUCKET_NAME = "gaia-scratch"     # or "$path" to serve several buckets as /<bucket>/<key>
   ALLOW_LIST_BUCKET = "false"      # "true" lets clients list objects (e.g. for obstore/fsspec ls)
   RCLONE_DOWNLOAD = "false"
   ```

   With the default `*.workers.dev` URL, `BUCKET_NAME` must be a fixed name or `$path`. The `$host` option needs a custom domain.

4. **Store the secret key (not in `wrangler.toml`) and deploy:**

   ```bash
   echo "<applicationKey from step 1>" | npx wrangler secret put B2_APPLICATION_KEY
   npx wrangler deploy
   ```

   The deploy prints a URL like `https://cloudflare-b2.<your-subdomain>.workers.dev`. Then `https://cloudflare-b2.<sub>.workers.dev/<key>` serves `s3://gaia-scratch/<key>`.

5. **Make responses cacheable.** Cloudflare won't cache responses to requests carrying an `Authorization` header, which the Worker adds. Set a default `Cache-Control` on the bucket:

   ```bash
   uvx b2 bucket update gaia-scratch allPrivate --bucket-info '{"Cache-Control": "public, max-age=3600"}'
   ```

   Keep `max-age` short for scratch data that gets overwritten.

### Notes and caveats

- **Use a custom domain for real use.** A `workers.dev` URL is fine for testing. Edge caching and (per third-party write-ups) the free-egress arrangement apply to traffic on a domain you've added to Cloudflare, e.g. a cheap domain via Cloudflare Registrar. Attach it under the Worker's *Settings → Domains & Routes*.
- **CORS:** the Worker forwards the browser's `Origin` header to B2, so the bucket's CORS rules ([cors-rules.json](cors-rules.json)) still apply. Test STAC Browser against the Worker URL. If preflight (`OPTIONS`) requests fail, add CORS handling in the Worker itself.
- **Large range reads:** for files over ~2 GB, Cloudflare may return the whole file instead of the range. The sample Worker detects this and retries.
- **STAC hrefs** should point at the Worker URL, not `s3.us-west-004.backblazeb2.com`, since the latter now requires credentials.
- **Writing** still goes directly to B2 with the `b2` AWS profile (see "Configure an AWS profile"). The Worker is read-only.

## References

- [Pricing](https://www.backblaze.com/cloud-storage/pricing)
- [Create and manage buckets](https://www.backblaze.com/docs/cloud-storage-create-and-manage-buckets) (public bucket requirements)
- [`b2 bucket create`](https://b2-command-line-tool.readthedocs.io/en/master/subcommands/bucket_create.html)
- [S3-compatible API intro](https://www.backblaze.com/apidocs/introduction-to-the-s3-compatible-api)
- [S3 Put Lifecycle Configuration](https://www.backblaze.com/apidocs/s3-put-lifecycle-configuration)
- [Enable CORS with the S3 API](https://www.backblaze.com/docs/cloud-storage-enable-cors-with-the-s3-compatible-api)
- [Bucket access logs](https://www.backblaze.com/docs/cloud-storage-bucket-access-logs)
- [Deliver private B2 content through Cloudflare Workers](https://www.backblaze.com/b2/docs/cloudflare_workers.html)
- [cloudflare-b2 sample Worker](https://github.com/backblaze-b2-samples/cloudflare-b2)
- [Cloudflare Workers limits](https://developers.cloudflare.com/workers/platform/limits/)
