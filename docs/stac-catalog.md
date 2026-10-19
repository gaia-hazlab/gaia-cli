# STAC catalog of staged stores

Every `gaia stage …` and `gaia align` run writes a
[STAC Item](https://github.com/radiantearth/stac-spec/blob/master/item-spec/item-spec.md)
next to its Zarr store. `gaia stage` writes both into its output folder, with the
resolved `stage.yaml` ([stage-file.md](stage-file.md#output-folder)). `gaia catalog rebuild` collects these Items into a
`catalog.json` that you can open in tools such as
[stac-browser](https://radiantearth.github.io/stac-browser/). This is the same
setup the PRISM cache uses.

```
s3://gaia/scratch/
  catalog.json            # written by `gaia catalog rebuild -p s3://gaia/scratch`
  skagit-hrrr/             # written by `gaia stage hrrr -o …/skagit-hrrr`
    skagit-hrrr.zarr/
    skagit-hrrr.json      # STAC Item
    stage.yaml
  skagit-daily/            # written by `gaia align -o …/skagit-daily`
    skagit-daily.zarr/
    skagit-daily.json
```

## The Item

The Item is built from the store's provenance attributes (see
[provenance-metadata.md](provenance-metadata.md)) by `provenance.stac_item`.

| Field | Source |
|---|---|
| `id` | Store name without `.zarr` |
| `bbox`, `geometry` | WGS84 extent of the data: raster cell edges and station points across every group. The AOI isn't used, because clipping keeps every touched cell. |
| `properties` | Root attributes (`created`, `gaia:cli-version`, `gaia:command`, …), plus `start_datetime`/`end_datetime` spanning every group |
| `properties.license` | The root's license, or the single license shared by every group. With several licenses it is omitted, and `gaia:license-note` explains why. |
| `properties.providers` | Providers from every group, one entry per name |
| `links` | Every group's `license`, `about`, `via` and `derived_from` links, once each |
| `assets.data` | `./<name>.zarr`, type `application/vnd+zarr`, with [xarray-assets](https://github.com/stac-extensions/xarray-assets) `xarray:open_kwargs` |
| `assets.config` | `gaia stage` only: `./stage.yaml`, type `application/yaml`, role `metadata`, the resolved stage file that reproduces the store |
| `stac_extensions` | Only extensions whose fields are present (their schemas reject an Item that declares an extension without its fields) |

## Rebuilding the catalog

```bash
gaia catalog rebuild -p s3://gaia/scratch    # needs AWS_PROFILE with write access
gaia catalog rebuild -p ./outputs            # local directories work too
```

The rebuild does four things:

1. Walks the prefix for `.json` files, without descending into `*.zarr`
   directories.
2. Keeps the gaia-cli Items, skipping any other JSON.
3. Drops Items whose store no longer has a `zarr.json`, for example scratch data
   that was deleted.
4. Writes `catalog.json`.

All hrefs are relative: the asset points at `./<name>.zarr` and the catalog
links to `./path/<name>.json`. The catalog therefore works however the prefix
is served (Kopah, B2, AWS, or a local directory), with no `s3://` to HTTPS
mapping. Browsing public buckets in stac-browser needs CORS on the bucket (see
[cors-rules.json](cors-rules.json)).

Staging never updates `catalog.json` itself. Each run writes only its own Item,
so concurrent runs can't overwrite each other's catalog links. Run
`gaia catalog rebuild` whenever the catalog should catch up, by hand or on a
schedule.
