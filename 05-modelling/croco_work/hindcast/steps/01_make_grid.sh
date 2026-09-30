#!/bin/bash
# =============================================================================
# Step 01 -- the model grid
#   hindcast/steps/01_make_grid.sh NAME
#
#   grid.ini      <- gtools/config/make_grid_config.py  (box + 1/RES_INV from domain.cfg;
#                    etopo2 bathymetry + GSHHS coastline from DATASETS_CROCOTOOLS)
#   croco_grd.nc  <- code/croco_pytools/prepro/make_grid.py (mask, smoothing rx0<=0.2)
# Output: hindcast/configs/NAME/grid.ini, hindcast/scratch/NAME/CROCO_FILES/croco_grd.nc
# The real size (xi_rho x eta_rho) is read from the file by the later steps --
# param.h gets LLm0 = xi_rho-2, MMm0 = eta_rho-2 automatically.
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"

for f in "${CROCO_DATA_ROOT}/DATASETS_CROCOTOOLS/Topo/etopo2.nc" \
         "${CROCO_DATA_ROOT}/DATASETS_CROCOTOOLS/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp"; do
    [ -f "$f" ] || die "missing $f -- install DATASETS_CROCOTOOLS (install/README.md, section 3)"
done
ls "${CROCO_PYTOOLS_DIR}"/prepro/Modules/toolsf*.so >/dev/null 2>&1 \
    || die "croco_pytools Fortran tools not compiled -- install/README.md (install/04)"

say "grid.ini for ${CONFIG_NAME}: lon ${LON_MIN}..${LON_MAX}, lat ${LAT_MIN}..${LAT_MAX}, 1/${RES_INV} deg"
cd "${GTOOLS_DIR}/config"
python3 make_grid_config.py "${CONFIG_NAME}" ${LON_MIN} ${LON_MAX} ${LAT_MIN} ${LAT_MAX} ${RES} ${RES}

say "make_grid.py (topography interpolation + smoothing iterations) ..."
rm -f "${CF}/croco_grd.nc"
cd "${CROCO_PYTOOLS_DIR}/prepro"
python3 make_grid.py "${CONFIG_DIR}/grid.ini" 2>&1 | tail -6
[ -f "${CF}/croco_grd.nc" ] || die "croco_grd.nc was not written"

say "CHECK -- the real dimensions (path must say hindcast/scratch):"
ls -l "${CF}/croco_grd.nc"
ncdump -h "${CF}/croco_grd.nc" | grep -E "xi_rho =|eta_rho ="
say "next: hindcast/steps/02_check_grid.sh ${CONFIG_NAME}"
