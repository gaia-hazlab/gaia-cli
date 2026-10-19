# Soil properties

Module: `src/gaia_cli/soil.py` · Command: `gaia stage soil`

## Source

`--product solus100` (the default and, for now, only product) is
[SOLUS100](https://doi.org/10.15482/USDA.ADC/25033856), the USDA Soil Landscapes of
the United States maps ([Nauman et al., 2024](https://doi.org/10.1002/saj2.20769)).
They are public Cloud-Optimized GeoTIFFs in a Google Cloud Storage bucket
(`https://storage.googleapis.com/solus100pub/`, [index](https://storage.googleapis.com/solus100pub/index.html)),
read anonymously over HTTPS with rioxarray. No key, no new dependency. License
CC-BY-4.0.

Every layer shares one grid: EPSG:5070 (CONUS Albers), 100 m. Coverage is CONUS only.
The extent check skips AOIs outside it.

## Layers

SOLUS publishes about 500 layers: about 20 properties at 7 depths, each with a
prediction and three uncertainty layers. `soil` stages only the seven the
[landlab workflow](https://github.com/gaia-hazlab/landlab_data_prep/blob/master/docs/landlab_workflow_contract.md)
reads: the predictions (`_p`) at 0 cm, plus depth to bedrock, which has no depth.
Variables keep SOLUS's property names:

| variable | SOLUS file | stored as | units |
|---|---|---|---|
| `anylithicdpt` | `anylithicdpt_cm_p.tif` | uint8 | cm |
| `cec7` | `cec7_0_cm_p.tif` | uint16, ×10 | meq/100g |
| `claytotal` | `claytotal_0_cm_p.tif` | uint8 | percent (mass) |
| `dbovendry` | `dbovendry_0_cm_p.tif` | uint16, ×100 | g cm-3 |
| `ph1to1h2o` | `ph1to1h2o_0_cm_p.tif` | uint16, ×100 | pH |
| `sandtotal` | `sandtotal_0_cm_p.tif` | uint8 | percent (mass) |
| `silttotal` | `silttotal_0_cm_p.tif` | uint8 | percent (mass) |

Each variable has `long_name`, `units` and `depth` (`0 cm`, or `all` for
`anylithicdpt`) attributes.

## Scaling

The COGs store integers multiplied by a scalar, for example pH 5.64 as 564. Only the
layer table, [`Final_Layer_Table_20231215.csv`](https://storage.googleapis.com/solus100pub/Final_Layer_Table_20231215.csv),
records the scalar. The files' own scale is 1, so reading a COG directly gives values
10 or 100 times too large. `load()` divides by the scalar. In the store, each
variable is packed back to its source integers with CF `scale_factor` and
`_FillValue`, so it is lossless and as compact as the source. xarray decodes it on
read to float64 in physical units, the same values `load()` returns. Nodata (255 or
65535) and cells outside the AOI are NaN.

The layer table's file names for depth to bedrock include `_2D_`
(`anylithicdpt_cm_2D_p.tif`), but the bucket's don't. `soil.py` lists the URLs
directly.

## Notes

- Like every raster, soil stays on its native grid, and `gaia align` puts it on a
  common grid. Its default regrid method is `bilinear` (`REGRID_METHOD` in
  `soil.py`), as in the landlab contract; `--method` overrides it. The contract
  also closes a node wherever any one layer is missing. That rule belongs in
  `gaia convert landlab` (roadmap step 3.1), as do the Landlab field names
  (`clay__total`, …) and the conversion of `soil__thickness` to metres.
- Not staged: other depths, other properties, and the prediction intervals (`_l`,
  `_h`, `_rpi`). A `--variables` option can add them if a workflow needs them. The
  [solus-stac](https://github.com/gaia-hazlab/solus-stac) prototype catalogs all
  of them. It isn't hosted, and seven fixed URLs don't need a catalog.
- SOLUS100 is one release built from soil observations of many dates, so the STAC
  Item has no time range, and its `datetime` is when the data was staged, as for
  `dem/usgs10m`.
