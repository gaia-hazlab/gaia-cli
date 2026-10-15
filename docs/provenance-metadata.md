# Provenance metadata in staged Zarr stores

Every store written by `gaia stage …` or `gaia align` records where its data
came from and how it was made. The metadata is stored as Zarr group attributes
named after STAC fields: core common metadata where a field exists, otherwise
the processing, version and scientific extensions. Facts with no STAC field use
a `gaia:` prefix. CF attributes on variables (`units`, `cell_methods`,
`grid_mapping`, `bounds`, …) are unchanged.

The code is in [`provenance.py`](../src/gaia_cli/provenance.py).

## Where it goes

There are two kinds of attributes.

- **Run attributes** say who made the store, when, and how. The `io.save_*`
  functions add them to the root group when they write the store.
- **Source attributes** describe one upstream product. Each product's `load()`
  sets them on the Dataset it returns.

| Store | Root group | Product groups |
|---|---|---|
| `gaia stage <product>` | run + source attributes | n/a (single group) |
| `gaia stage all` | run attributes, plus `gaia:license-note` | source attributes on `rasters/hrrr`, `rasters/prism`, `stations/synoptic` |
| `gaia align` | new run attributes; the input's run moves to `gaia:inputs` | carried through; time extent updated after aggregation |

## Run attributes (root)

| Attribute | Meaning |
|---|---|
| `created` | Write time, RFC 3339 UTC (STAC common metadata). |
| `gaia:cli-version` | Commit SHA of the gaia-cli that wrote the store (see below). |
| `gaia:command` | The `gaia …` command line. Only present when the CLI wrote the store, not when `io.save_*` is called from Python. |
| `stac_extensions` | Schema URLs for the STAC extensions whose fields appear in the root group (omitted if none). |
| `gaia:license-note` | Present when the product groups have different licenses. The root then has no `license`. |
| `gaia:inputs` | `gaia align` only: the run attributes and `href` of each input store, oldest first. |
| `links` | `gaia align` only: a `derived_from` link to the input store. |

`gaia:cli-version` is resolved in this order:

1. `git rev-parse HEAD`, used only if the package is in a gaia-cli checkout.
   This covers the usual editable pixi install.
2. The `vcs_info.commit_id` that `pip install git+https://…` records.
3. Otherwise, the package version, or `unknown` if that can't be found either.

The commit is recorded whether or not the working tree has uncommitted changes.
Once releases are tagged, the package version (step 3) becomes the semver tag.

## Source attributes (per product)

Every product carries `title`, `description`, `license` (an SPDX id, or
`other`), `providers`, `links`, `processing:lineage`, and
`start_datetime`/`end_datetime`. The time range comes from `time_bounds` when
present, so it covers each value's whole interval.

### HRRR

These come from dynamical.org's STAC collection
(`https://stac.dynamical.org/noaa-hrrr-analysis/collection.json`), fetched on
every run.

- **From the collection:** the title, description, providers, and the `license`
  and `about` links.
- **From the Icechunk store itself:** `version` (from `dataset_version`) and
  `license`. The collection describes dynamical's latest version, which can be
  newer than the one in `hrrr.DEFAULT_PREFIX`.
- **`gaia:icechunk-snapshot`:** the snapshot ID that was read. This pins the
  exact upstream state, because the `main` branch keeps changing.
- **Links:** a `derived_from` link to the collection and a `via` link to the
  Icechunk repository.

### PRISM

These come from the gaia PRISM cache's own STAC collection
(`<cache>/collection.json`), which already names PRISM as producer and
licensor, the gaia hazlab as host, and links to the PRISM terms of use. PRISM
release numbers vary by day, so they stay per time step in the variable's
`prism_release` attribute.

### Synoptic

Synoptic publishes no STAC collection, so these are written in
`synoptic.PROVENANCE`.

- `license` is `other`. `gaia:license-note` asks users to check Synoptic's terms
  and the station networks' terms before sharing.
- `version` is the API version (`v2`).
- The `via` link is the API endpoint without query parameters, so the token is
  never recorded. A test guards this.

### NLCD

These are written in `nlcd.source_attrs(year)`.

- **Providers:** USGS is producer and licensor; USGS HyTEST hosts the mosaic.
- **`license`:** `other`, linking to the USGS copyright policy, because U.S.
  public domain has no SPDX id.
- **`version`:** `C1V1`.
- **Time range:** the calendar year.
- **`via` link:** the COG that was read.

## Fetch failures

`provenance.from_stac_collection` logs a warning if the upstream collection
can't be fetched, and then writes only the `derived_from` link. A metadata
outage never stops a staging run.

## Example (`gaia stage hrrr`, root group)

```json
{
  "title": "NOAA HRRR analysis",
  "license": "CC-BY-4.0",
  "providers": [
    {"name": "NOAA NWS NCEP", "roles": ["producer"], "url": "https://www.weather.gov/ncep/"},
    {"name": "dynamical.org", "roles": ["processor", "host"], "url": "https://dynamical.org/"}
  ],
  "version": "0.2.0",
  "links": [
    {"rel": "license", "href": "https://creativecommons.org/licenses/by/4.0/", "type": "text/html", "title": "CC-BY-4.0"},
    {"rel": "about", "href": "https://dynamical.org/catalog/noaa-hrrr-analysis/", "type": "text/html"},
    {"rel": "derived_from", "href": "https://stac.dynamical.org/noaa-hrrr-analysis/collection.json", "type": "application/json"},
    {"rel": "via", "href": "s3://dynamical-noaa-hrrr/noaa-hrrr-analysis/v0.2.0.icechunk/", "title": "Icechunk repository"}
  ],
  "gaia:icechunk-snapshot": "MFB4FGH0W7XY3DECSAA0",
  "processing:lineage": "precipitation_surface clipped to the AOI on its native grid",
  "start_datetime": "2025-11-30T23:00:00Z",
  "end_datetime": "2025-12-01T23:00:00Z",
  "created": "2026-10-15T21:36:22Z",
  "gaia:cli-version": "f407c4fce48da63fe26200c5c4cb2f90230bf59e",
  "gaia:command": "gaia stage hrrr -i …/SkagitBoundary.json -s 2025-12-01 -e 2025-12-01 -o skagit_hrrr.zarr",
  "stac_extensions": ["https://stac-extensions.github.io/processing/v1.2.0/schema.json", "…"]
}
```

## Not done (by choice)

- No CF `history` attribute.
- No AOI geometry or link in the attributes.
- No versions of dependencies such as xarray or zarr; only gaia-cli is recorded.
- No flag for uncommitted changes in the working tree.

Each store also gets a STAC Item next to it, built from these attributes; see
[stac-catalog.md](stac-catalog.md).
