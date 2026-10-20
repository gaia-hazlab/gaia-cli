# Soil properties

Module: `src/gaia_cli/datasets/soil.py` · Command: `gaia stage soil`

## Source

Two products, each its own group (`soil/<product>`), so a stage file can list both:

- `--product solus100` (the default): the seven SOLUS100 layers. See below.
- `--product statsgo`: the STATSGO K-factor. See [STATSGO K-factor](#statsgo-k-factor).

`--product solus100` is
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

## STATSGO K-factor

`--product statsgo` stages one variable, `kffact`, the soil erodibility factor the
[landlab contract](https://github.com/gaia-hazlab/landlab_data_prep/blob/master/docs/landlab_workflow_contract.md)
requires. It is the same input the USGS postfire debris-flow likelihood model uses: the
M1 model's soil term is KF-factor × rainfall
([Staley et al., 2016](https://pubs.usgs.gov/of/2016/1106/ofr20161106.pdf);
[Staley et al., 2017](https://doi.org/10.1016/j.geomorph.2016.10.019)), and USGS's
[pfdf](https://ghsc.code-pages.usgs.gov/lhp/pfdf/api/data/usgs/statsgo.html) reads it
from this dataset.

### Source

[King (2025)](https://doi.org/10.5066/P13WAPYV), a USGS data release: the `KFFACT`
field of [Schwarz and Alexander (1995)](https://doi.org/10.5066/P94JAULO) rasterized as
one CONUS Cloud-Optimized GeoTIFF. It's read anonymously over HTTPS from the release's
public S3 copy:

`https://prod-is-usgs-sb-prod-publish.s3.amazonaws.com/6750c172d34ed8d3858534d8/statsgo-KFFACT.tif`

The ScienceBase file link (`www.sciencebase.gov/catalog/file/get/6750c172d34ed8d3858534d8?name=statsgo-KFFACT.tif`,
as in landlab_data_prep's `base.example.yaml`) returns 404 since ScienceBase moved its
files to `sciencebase.usgs.gov`. The S3 URL is the item's `publishedS3Uri`. U.S. public
domain.

| | |
|---|---|
| grid | 30 m, float32, NaN nodata, nearest-neighbour overviews |
| CRS | NAD27 / CONUS Albers, EPSG:5069 (see [Datum](#datum)) |
| source scale | 1:250,000 polygons (STATSGO, October 1994), about 78,500 in CONUS |
| values | 0.02 to 0.62 (FGDC metadata range), and the code -0.1 |

### What KFFACT is

Schwarz and Alexander built one value per map unit: the `kffact` of each STATSGO
layer, averaged over the layers weighted by thickness, then over the map unit's
components weighted by their percentage. So it is a whole-profile average, not the
surface-horizon K that USLE/RUSLE use, and it is Kf (the fine-earth fraction, rock
fragments excluded), not Kw. A map unit is a 1:250,000 polygon, so an AOI usually
has only a few distinct values: five for the Pioneer Fire AOI.

**Missing values.** `-0.1` is the source's missing-value code: a map unit with no K
factor, which was set for every map unit whose components carry a water designation.
King (2025) describes it as marking large water bodies, but in the source it means
missing in general. `load()` sets it to NaN, as the landlab contract does (`-9999`).
On the Pioneer AOI aligned to NLCD, 87.5% of NLCD open-water cells are NaN in `kffact`
and 0.7% of the other cells are.

**Units.** `units` is `inches per hour`, copied from King (2025), which defines KFFACT
as "the saturated hydraulic conductivity of the fine soil (< 2mm) fraction". That
contradicts the source, which defines it as the "actual k factor used in water erosion
component of universal soil loss equation", and the NRCS definition of `kffact`, in
t·ac·h/(100·ac·ft·tonf·in) (× 0.1317 for t·ha·h/(ha·MJ·mm)). The values (0.02 to
0.62) are typical K factors.

That unit follows from USLE, A = R·K·LS·C·P: K is soil loss A per unit of rainfall
erosivity R. A is in short tons per acre (t/ac), and R is in hundreds of foot-tonf-inch
per acre-hour, the rainfall energy (foot-tons-force) times the peak 30-minute
intensity (inches per hour). Dividing gives (t/ac) / (100·ft·tonf·in/(ac·h)) =
t·ac·h/(100·ac·ft·tonf·in): t short tons, ac acres, h hours, ft·tonf foot-tons-force,
in inches. The acres cancel, but by convention the unit keeps both, to show where it
comes from.

A TODO in `soil.py` marks the attribute to check with the authors.

### Datum

The release really is NAD27, not a mislabelled NAD83. The STATSGO polygons were
compiled in 1994 to 1995 on NAD27, and the 1995 shapefiles' `.prj` says so. USGS's
`statsgo-cog` package rasterized them in the shapefiles' own CRS, and the 2025 FGDC
metadata records NAD27 / Clarke 1866.

NAD27 to NAD83 moves points by tens of metres in CONUS; at Pioneer, PROJ's NADCON
grid shifts them 35 m west and 70 m north. `gaia align` reprojects with that
transformation: rasterio and pyproj both use PROJ, and conda-forge's PROJ sets
`PROJ_NETWORK=ON`, so the grid is fetched on demand.

How much it matters: matching the `-0.1` water polygons to NLCD open water, the
offset left after alignment is smaller with the declared NAD27 than if the file were
read as NAD83 at Lake Tahoe (30 m vs 67 m) and Lake Chelan (42 vs 85 m), similar at
the Pioneer AOI (134 vs 153 m), and slightly larger at Flathead Lake (108 vs 90 m).
All of these are within the ~125 m horizontal accuracy of a 1:250,000 map, so the
datum is right but makes little practical difference. The polygon boundaries are
coarser than either.

### Alternatives

Checked in October 2026; none is a drop-in for the landlab and M1 input.

| source | K | resolution | access | notes |
|---|---|---|---|---|
| **King (2025), used here** | Kf, STATSGO 1995, whole profile | 30 m (1:250,000 polygons) | anonymous S3 COG | What M1 was calibrated on and what pfdf reads |
| [NRCS Soil Data Access](https://sdmdataaccess.sc.egov.usda.gov) (SSURGO, STATSGO2) | `kwfact`, `kffact` per horizon | survey polygons (SSURGO ~1:24,000) | anonymous SQL REST | Current. Surface K needs the dominant component and top horizon per map unit, then a join to a map unit raster |
| [Planetary Computer gNATSGO](https://planetarycomputer.microsoft.com/dataset/gnatsgo-rasters) | none pre-joined | 10 m map unit keys, NAD83 | anonymous SAS token | The raster to join SDA's K to. Its `chorizon` table, where K lives, is almost empty (about 1,000 rows), so SDA is needed |
| [ISSR-800](https://casoilresource.lawr.ucdavis.edu/soil-properties/) (UC Davis) | Kw, 0–25 cm | 800 m | WCS, CC-BY-4.0 | Coarse cross-check |
| [ESDAC global soil erodibility](https://esdac.jrc.ec.europa.eu/content/global-soil-erodibility) | Several K estimates | 1 km, global | registration, no redistribution | Not usable for anonymous staging |
| SOLUS100, POLARIS | none | | | K could be derived from texture and organic carbon (nomograph), as a new product |

A surface K from SSURGO (gNATSGO map unit keys plus Soil Data Access) is the natural
upgrade if a workflow wants current, finer K and accepts that it isn't what M1 was
calibrated on. It would also cover Alaska, which King (2025) doesn't. It stays
deferred until a workflow asks for it.

