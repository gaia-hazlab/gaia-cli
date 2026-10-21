#!/usr/bin/env bash
# Landlab inputs from an aligned store, with GDAL only: one GeoTIFF and one ESRI
# ASCII grid per layer, float32 with -9999 for nodata, following the landlab
# contract (https://github.com/gaia-hazlab/landlab_data_prep/blob/master/docs/landlab_workflow_contract.md).
#
#   scripts/landlab-inputs.sh STORE OUTDIR
#
# STORE can be remote, e.g. /vsis3/gaia/scratch/public/pioneer-10m/pioneer-10m.zarr with
# AWS_S3_ENDPOINT=s3.kopah.uw.edu AWS_NO_SIGN_REQUEST=YES AWS_VIRTUAL_HOSTING=FALSE.
#
# STORE must already be on the contract grid, e.g.
#   gaia align -i /tmp/pioneer/pioneer.zarr -o /tmp/pioneer-10m --crs utm --res 10 --aoi $AOI \
#     -m dem/usgs10m=cubic -m soil/solus100=bilinear \
#     -m rasters/nlcd=nearest -m rasters/burn-severity=nearest
# so this only converts formats and applies the contract's nodata rules. Layers
# the model reads directly get their Landlab field names; the soil layers keep
# SOLUS's names, since landlab_data_prep derives the soil fields from them.
# Needs GDAL >= 3.12 (the `gdal` command).
set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 STORE OUTDIR" >&2
  exit 1
fi
store=$1
out=$2
mkdir -p "$out"
# Statistics would otherwise write .aux.xml files into the store
export GDAL_PAM_ENABLED=NO

# GDAL's Zarr driver looks for <dim>/zarr.json in every directory above the store
# and warns for each one an S3 server refuses (dev/gdal-zarr/upstream-issue.md),
# as it does for missing external overviews (zarr.json.<name>.ovr).
# Drop those lines only, so real errors still show.
gdal() {
  command gdal "$@" 2> >(grep -vE '^(Warning 1: HTTP response code on .*/zarr\.json(\.[^/]*\.ovr)?: (0|403)|Warning 6: Unsupported option: https?|ERROR 5: Missing url parameter)$' >&2)
}

# output name, then the array in the store
layers=(
  "topographic__elevation dem/usgs10m/elevation"
  "landcover rasters/nlcd/landcover"
  "kffact soil/statsgo/kffact"
  "anylithicdpt soil/solus100/anylithicdpt"
  "cec7 soil/solus100/cec7"
  "claytotal soil/solus100/claytotal"
  "dbovendry soil/solus100/dbovendry"
  "ph1to1h2o soil/solus100/ph1to1h2o"
  "sandtotal soil/solus100/sandtotal"
  "silttotal soil/solus100/silttotal"
)

for layer in "${layers[@]}"; do
  read -r name array <<<"$layer"
  echo "$array -> $out/$name.{tif,asc}"
  # unscale applies the scale_factor of packed integers (soil), else only casts;
  # calc then writes every nodata cell (NaN or the integer fill) as -9999
  gdal raster pipeline -q read "ZARR:\"$store\":/$array" \
    ! unscale --ot Float32 \
    ! write --overwrite "$out/$name.vrt"
  gdal raster calc -q --overwrite -i "A=$out/$name.vrt" --calc A \
    --ot Float32 --nodata -9999 --propagate-nodata "$out/$name.tif"
  gdal raster convert -q --overwrite "$out/$name.tif" "$out/$name.asc"
  rm "$out/$name.vrt"
done

# Burn severity has no nodata: classes 2-4 are kept, and everything else
# (unburned to very low, unmapped, outside the AOI) becomes 1. materialize
# fixes the classes before edit removes the nodata reclassify reads.
echo "rasters/burn-severity/burn_severity -> $out/burn__severity.{tif,asc}"
gdal raster pipeline -q read "ZARR:\"$store\":/rasters/burn-severity/burn_severity" \
  ! reclassify --ot Float32 -m "[2,4]=PASS_THROUGH; NO_DATA=1; DEFAULT=1" \
  ! materialize \
  ! edit --nodata none \
  ! write --overwrite "$out/burn__severity.tif"
gdal raster convert -q --overwrite "$out/burn__severity.tif" "$out/burn__severity.asc"
