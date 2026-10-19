# Digital elevation models

Module: `src/gaia_cli/dem.py` · Command: `gaia stage dem`

## Sources

Both products are public VRT mosaics of Cloud-Optimized GeoTIFFs, read anonymously
over HTTPS with rioxarray (GDAL). They need no API key and no new dependency.

| `--product` | source | grid | heights | coverage |
|---|---|---|---|---|
| `cop30` (default) | [Copernicus DEM GLO-30](https://doi.org/10.5069/G9028PQB), hosted by OpenTopography: `https://opentopography.s3.sdsc.edu/raster/COP30/COP30_hh.vrt` | 1 arc-second (~30 m), EPSG:4326 | EGM2008 (EPSG:3855) | global |
| `usgs10m` | [USGS 3DEP 1/3 arc-second](https://doi.org/10.5069/G98K778D), hosted by USGS (The National Map): `https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/USGS_Seamless_DEM_13.vrt` (`s3://prd-tnm/...`) | 1/3 arc-second (~10 m), EPSG:4269 (NAD83) | NAVD88 (EPSG:5703) | CONUS, Hawaii, territories, parts of Alaska |

The variable is `elevation` in metres (`standard_name: surface_altitude`). Cop30 is a
surface model (it includes trees and buildings). 3DEP is bare earth. Cop30 sets no
nodata value. 3DEP's `-999999` is read as NaN. An AOI outside 3DEP's coverage gives
an all-NaN array.

```bash
gaia stage dem -i $AOI -o /tmp/bolt-creek-dem                   # cop30
gaia stage dem -i $AOI -o /tmp/bolt-creek-dem10 --product usgs10m
```

In a stage file, `dem: {product: usgs10m}`. One stage file holds one `dem` entry
(group `rasters/dem`). To compare the two products, stage each into its own store.

## Vertical datum

Heights are relative to each product's geoid. They are not converted. The datum's EPSG
code is the variable's `vertical_datum` attribute. EGM2008 and NAVD88 differ by about a
metre over CONUS, so subtracting one product from the other includes that offset.
`gaia align` regrids horizontally only, and a TODO in `align.py` marks where a vertical
datum conversion would go.

## Notes

- Like every raster, the DEM stays on its native grid. `gaia align` puts it on a
  common grid.
- `processing:lineage`, `sci:doi`, `sci:citation` and a `cite-as` link record the
  source. The STAC Item's time range for Cop30 is the acquisition period
  (2011–2015). 3DEP blends source data of many dates, so it has none.
