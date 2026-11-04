# Viewing stores in a browser

Notes from opening a staged store in the Source Cooperative
[zarr-viewer](https://source-cooperative.github.io/zarr-viewer/) (tested 2026-11-02),
and options for opening it from stac-browser.

## What the zarr-viewer needs

The viewer reads Zarr in the browser with [zarrita](https://github.com/manzt/zarrita.js).
For example, the aligned Pioneer store (the `min_zoom` part is explained below):

<https://source-cooperative.github.io/zarr-viewer/?url=https%3A%2F%2Fs3.kopah.uw.edu%2Fgaia-scratch%2Fgaia-cli%2Fpioneer-10m%2Fpioneer-10m.zarr&min_zoom=8>

- **404 for missing keys.** zarrita's `zarr.open` works out the Zarr version by trying v2
  (`.zattrs`, `.zgroup`) before v3, and moves on only after a 404. Its `FetchStore`
  throws on any other status. Anonymous requests to `s3://gaia/scratch/public/` get 403
  for missing keys, because anonymous users may list only that prefix, not the whole
  bucket. So the viewer fails there with "403 Forbidden". The fully public
  `gaia-scratch` bucket returns 404 and works. Icechunk has the same problem; see
  nisar-soilmoisture-viewer's `Kopah+icechunk.md`.
- **CORS** allowing GET from any origin, with the `Range` request header allowed and
  `Content-Range` and `Content-Length` exposed. `gaia-scratch` on Kopah sends these.
- **GeoZarr grid attributes.** The default `scalar-grid` profile places a grid using the
  `proj:` and `spatial:` attributes on the root group or the variable's parent group.
  Without them it assumes the `y`/`x` coordinates are EPSG:4326 degrees.
  `io._geozarr_attrs` writes these on every raster group. CF `grid_mapping` and `_CRS`
  are kept for xarray and GDAL.
- **`min_zoom`.** The viewer works out one minimum zoom from the first variable, assuming
  a full-resolution store with no overviews. Its budget is about 8 MB and 16 chunks per
  256 px tile. It doesn't account for how small the store is, so a 10 m store gates at
  z12. A store whose first variable is in EPSG:4269 gates at z24: the viewer recognizes
  only 4326 as a degrees CRS, so it reads degrees as metres. Add `&min_zoom=8` to the
  URL. A store whose spatial plane is a single chunk skips the gate.
  - The proper fix is overviews, using the GeoZarr
    [multiscales](https://github.com/zarr-conventions/multiscales) convention, which the
    viewer's `multiscale-grid` profile renders. This is not implemented.

## Opening the viewer from stac-browser

stac-browser shows the store's footprint but can't render Zarr. There are two ways to
get from an Item to the viewer.

1. **A link on the Item. Works on the public stac-browser.** stac-browser lists every
   Item link and opens any link that isn't to a STAC object in a new tab
   ([metadata.md](https://github.com/radiantearth/stac-browser/blob/main/docs/metadata.md#customizing-links)):

   ```json
   {"rel": "related", "type": "text/html", "title": "Open in zarr-viewer",
    "href": "https://source-cooperative.github.io/zarr-viewer/?url=<encoded store URL>&min_zoom=8"}
   ```

   - The link appears in the link list, not on the Zarr asset.
   - Don't use `rel: "alternate"`: stac-browser reads it as the language extension.
   - The viewer needs the store's absolute HTTPS URL, but Item hrefs are relative
     ([stac-catalog.md](stac-catalog.md)). gaia would need the public base URL, for
     example a `--public-url` option on `gaia catalog rebuild` or on stage and align.
2. **A button on the Zarr asset. Needs a stac-browser change.** stac-browser calls these
   buttons [actions](https://github.com/radiantearth/stac-browser/blob/main/docs/actions.md).
   They are built in and enabled at build time in `assetActions.config.js`, and none
   handles Zarr.
   - A `ZarrViewer` action would be about 20 lines, like `src/actions/assets/GeoJsonIo.js`.
     It would show for `asset.type === "application/vnd+zarr"` and pass the asset's
     absolute href as `?url=`.
   - Upstream: a PR to radiantearth/stac-browser that enables it by default. Several
     third-party viewers are already enabled this way, and the public instance would
     then get the button.
   - Otherwise, deploy our own stac-browser with the action enabled.
   - A generic action can't choose `min_zoom` for each store. Without overviews,
     stores would load with the zoom gate.

Option 1 is the quickest. Try option 2 upstream before deploying our own stac-browser.
