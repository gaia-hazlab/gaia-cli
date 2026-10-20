# Burn severity

Module: `src/gaia_cli/datasets/burn_severity.py` · Command: `gaia stage burn-severity`

## Source

Soil burn severity (SBS) maps from USDA Forest Service Burned Area Emergency Response
(BAER) assessments: a satellite burn severity classification (BARC, from a dNBR image)
that the BAER team checks and corrects in the field. They are read from the national
[BAER Soil Burn Severity image service](https://imagery.geoplatform.gov/iipp/rest/services/Fire_Aviation/USFS_EDW_BAER_SoilBurnSeverityClassification/ImageServer)
on GeoPlatform, an ArcGIS ImageServer holding one raster per fire
([data.gov entry](https://catalog.data.gov/dataset/baer-soil-burn-severity-classification-image-service)).
Plain `requests`: no key, no new dependency. U.S. public domain.

Coverage is the Forest Service's own assessments, **2017 on** (about 800 fires in
2026, 60–150 a year), within CONUS. Fires assessed only by Interior agencies, and
older fires, are missing.

## How a fire is chosen

`burn-severity` is static and driven by the AOI, like the other static datasets:

1. The service's catalog is queried for the fire rasters whose footprints intersect
   the AOI.
2. If there is exactly one, it is used. If there are several (the Pioneer reference
   AOI meets five fires from 2017 to 2025), `fire` and `year` choose one. `fire` is
   part of the BAER name, ignoring case and with spaces as underscores
   (`bolt creek` matches `bolt_creek_sbs`). An ambiguous or empty choice is an error
   that lists the fires in the AOI.
3. The raster is exported with `exportImage` on its own grid (its CRS, origin and
   cell size from the raster's `/info`), nearest neighbour and locked to that one
   raster, so nothing is resampled. The export for Bolt Creek and for Pioneer is
   identical to the GeoTIFF in the USGS download zip within the fire's classes.

```bash
gaia stage burn-severity -i pioneer-aoi.geojson -o /tmp/pioneer-sbs --fire pioneer
```

```yaml
datasets:
  burn-severity: {fire: bolt creek}     # the Bolt Creek AOI also meets Irving Peak (2022)
```

Each fire keeps its own grid: usually 20 m (Sentinel-2) or 30 m (Landsat), in a UTM or
regional Albers CRS. `gaia align` puts it on a common grid.

## Values

The variable `burn_severity` has four classes, marked categorical with CF flags
(`flag_values`, `flag_meanings`, and `flag_colors` from the BAER legend), so `gaia align`
regrids it with `mode`:

| value | class |
|---|---|
| 1 | unburned to very low |
| 2 | low |
| 3 | moderate |
| 4 | high |

Everything else is outside the mapped area: the service writes 15 there, and some
source files also use 0. Those cells, and cells outside the AOI, are NaN. The store
keeps the source bytes: uint8 with `_FillValue` 15.

The STAC Item's time range is the fire's year, the only date the service records. The
pre- and post-fire image dates are in the metadata text of the USGS download zips,
not in the service.

## Notes

- **Not used: the USGS zip downloads**
  (`https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/baer/<name>.zip`),
  which landlab reads. The server refuses directory listings, so a file can only be
  found by guessing its name, from the fire name and sometimes the incident id
  (`pioneer_wa4818312053120240608_sbs.zip`). Older names follow no pattern
  (`Tucker_SBS`, `Final_SBS_Brooklyn`). In a sample of the service's catalog, only
  half the names had a zip.
- The landlab workflow sets every value except 2–4, missing data included, to 1
  (unburned). That rule belongs in `gaia convert landlab` (roadmap step 3.1).
- dNBR is not staged; nothing reads it.
- MTBS is a different product (Landsat burn severity of the vegetation, 1984 on, not
  field-checked soil severity), with
  [similar image services](https://imagery.geoplatform.gov/iipp/rest/services/Fire_Aviation/USFS_EDW_MTBS_CONUS/ImageServer).
  It could fill in for older fires if a workflow needs it.
- The service refuses (403) a spatial reference given as WKT in the query string, so
  requests use the raster's WKID.
