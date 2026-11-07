# Earthmover marketplace datasets: access and usage terms

Notes from 2026-11-06, on how `gaia stage era5` reads Earthmover's ERA5 and whether it
could read it without an account. Open question for Earthmover at the end.

## How gaia reads it today

`era5.py` opens `earthmover-public/era5` through the Arraylake client:
`arraylake.Client().get_repo("earthmover-public/era5").readonly_session("main")`. The
client asks Earthmover's API for the repo, which needs a login: `ARRAYLAKE_TOKEN` (a
free Earthmover API key) or a cached `arraylake auth login`. The API returns the
storage location and credentials. See [era5.md](era5.md).

HRRR is different: `hrrr.py` opens dynamical.org's Icechunk repo straight from its
public S3 bucket with `icechunk.s3_storage(anonymous=True)`, with no account.

## Anonymous access works

The repo lives in a plain S3 bucket, which allows anonymous reads:

| | |
|---|---|
| Bucket | `earthmover-icechunk-era5` |
| Prefix | `icechunkV2` |
| Region | `us-east-1` |

With `ARRAYLAKE_TOKEN` unset, this opens the same snapshot as the Arraylake route
(`VKGPX45BPEQJMVZZXZEG` on 2026-11-06) and read 48 hours of `tp` at Peters Dome in 3 s:

```python
import icechunk
import xarray as xr

storage = icechunk.s3_storage(
    bucket="earthmover-icechunk-era5", prefix="icechunkV2",
    region="us-east-1", anonymous=True,
    endpoint_url="https://s3.us-east-1.amazonaws.com",  # only if an AWS profile sets another
)
session = icechunk.Repository.open(storage).readonly_session("main")
ds = xr.open_zarr(session.store, group="single/temporal", zarr_format=3, chunks=None)
ds.tp.sel(latitude=63.25, longitude=208.75, valid_time="2023-09-13")  # lon 0–360, m
```

## What Earthmover's docs say

Nothing documents anonymous S3 access, and nothing forbids it explicitly. Every
official route assumes a login:

- [ERA5 blog post](https://www.earthmover.io/blog/low-latency-era5-icechunk-marketplace):
  the comparison table gives the access model as "Arraylake", key limitation "Login
  required"; the quickstart has `client.login()  # requires free Earthmover account`.
- [Marketplace docs for data users](https://docs.earthmover.io/marketplace/data-users):
  "Subscriptions are not anonymous. When you subscribe to a dataset, the data provider
  can see your organization name and contact email." Marketplace repos need
  `arraylake >= 0.26.1`.
- [Marketplace FAQ](https://docs.earthmover.io/marketplace/faq): "Anyone with an
  Arraylake account can subscribe to free listings instantly and start querying data
  right away."

[Terms of service](https://docs.earthmover.io/terms-of-service):

- Part II §3.2: "Unless Listing Terms expressly state otherwise, Customer may not
  modify, create derivative works of, reverse engineer, resell, sublicense, or transfer
  any Product, or circumvent any technical controls applied to it." Whether an open
  bucket is a "technical control" is unclear.
- Part I §1.3(b)(v) forbids attempts "to gain unauthorized access to the Earthmover
  Platform". Reading the bucket directly isn't addressed.
- Part II §2.2: "Where a Data Provider does not specify Listing Terms for a Product,
  Part III provides that the Product is made available under the Creative Commons
  Attribution 4.0 license."
- Part II §2.4 refers to Marketplace Policies, not read.

Free tier (blog and [pricing](https://docs.earthmover.io/pricing)):

- "There are no artificial quotas or limits on consumption."
- "Free and open for science, backed by a commercial SLA for industry."
- "The only difference between the free and paid versions is the latency; the free
  version is updated every three months, rather than daily."
- The Community tier includes "Unlimited access to free datasets in the Arraylake
  Marketplace".

Not confirmed: the ERA5 listing page
(`app.earthmover.io/marketplace/6a19bcfe9aa6e97720a2fad2`) needs a login, so its
license field and attribution text are unknown.

`earthmover-public/era5-surface-aws`, in Earthmover's
[sample data docs](https://docs.earthmover.io/sample-data/era5), is a different
dataset: 18 surface variables, 1975–2024, sourced from NCAR. gaia doesn't use it.
Checked 2026-11-06 via `gaia/era5-surface-aws-subscription`: it is ERA5 on the same
0.25° grid (721 × 1440), not ERA5-Land, with no lapse-rate adjustment. Its 22 variables
(`t2`, `d2`, `sp`, `cp`/`lsp`, ...) have no `tp` and no surface geopotential.

## Decision

Keep the Arraylake route for now. Anonymous reads work, but they aren't the
documented or intended route, and the terms are ambiguous.

`gaia stage` warns up front when `ARRAYLAKE_TOKEN` isn't set (`TOKENS` in `era5.py`,
PR #58).

## Ask Earthmover

support@earthmover.io, in writing:

1. Are anonymous reads of the `earthmover-public/era5` bucket
   (`s3://earthmover-icechunk-era5/icechunkV2`) intended and allowed, e.g. for an
   open-source tool that reads small AOIs?
2. What are the ERA5 listing's license and attribution text?

If anonymous reads are fine, switching `era5.py` to `icechunk.s3_storage(anonymous=True)`
as in `hrrr.py` drops `ARRAYLAKE_TOKEN`, the `arraylake` dependency, the CI secret and
the `AWS_ENDPOINT_URL_S3` workaround. The Icechunk snapshot id stays in provenance.
