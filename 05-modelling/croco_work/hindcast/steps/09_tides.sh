#!/bin/bash
# =============================================================================
# Step 09 -- OPTIONAL: tidal forcing from TPXO
#   hindcast/steps/09_tides.sh NAME [RUN_DATE]      (default RUN_DATE = TEST_START)
#
# Builds CROCO_FILES/croco_frc.nc for a run starting at RUN_DATE:
#   * parameters: hindcast/configs/NAME/crocotools_param_tides.py (written by step 02;
#     TPXO7, reader key 'tpxo7_croco', 10 waves M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm,
#     currents + potential + nodal corrections)
#   * make_tides runs in its OWN gen dir (tide_gen/CROCO_FILES): its inputdata is a
#     TPXO key that would clash with the GLORYS 'mercator' crocotools_param.py
#   * phases are referenced to Yorig; the nodal corrections are evaluated at RUN_DATE,
#     so a tide file belongs to ONE run -- the driver rebuilds it for every cycle
# Checks: 8 variables x 10 waves, M2 amplitude, NO fill value on ocean cells.
# Tides are then switched on at compile time (07_compile.sh NAME plain_tides) and
# the output is retimed to hourly history / daily averages (08 and the driver).
# Not needed for a tide-free hindcast (TIDES=0).
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
RUN_DATE=${2:-${TEST_START}}
[ -f "${TPXO_DIR}/TPXO7.nc" ] || die "${TPXO_DIR}/TPXO7.nc missing (DATASETS_CROCOTOOLS/TPXO7)"
[ -f "${CONFIG_DIR}/crocotools_param_tides.py" ] || die "run steps/02_check_grid.sh ${CONFIG_NAME} first"

TGEN=${HCAST}/tide_gen/CROCO_FILES
mkdir -p "${TGEN}"
cp "${CF}/croco_grd.nc" "${TGEN}/"
cp "${CONFIG_DIR}/crocotools_param_tides.py" "${TGEN}/crocotools_param.py"
rm -f "${TGEN}/croco_frc.nc"

say "make_tides (TPXO7, run date ${RUN_DATE}, Yorig ${YORIG}) ..."
cd "${GTOOLS_DIR}"
python ggosss26.py make_tides --input_dir "${TPXO_DIR}" --output_dir "${TGEN}" \
    --run_date "${RUN_DATE} 00:00:00" --Yorig "${YORIG}" --fname_out croco_frc.nc > "${TGEN}/make_tides.log" 2>&1 \
    || { tail -15 "${TGEN}/make_tides.log"; die "make_tides failed -- ${TGEN}/make_tides.log"; }
grep -E "Processing \*" "${TGEN}/make_tides.log" | tr -s ' \n' ' '; echo

say "CHECK tide file:"
python3 - <<PYEOF
import xarray as xr, numpy as np, sys
d = xr.open_dataset("${TGEN}/croco_frc.nc", decode_times=False)
g = xr.open_dataset("${TGEN}/croco_grd.nc")
m = g.mask_rho.values.astype(bool)
print("vars:", list(d.data_vars)); print("dims:", dict(d.sizes))
m2 = np.where(m, d.tide_Eamp.isel(tide_period=0).values, np.nan)
print("M2 amplitude on ocean: mean %.3f m  max %.3f m" % (np.nanmean(m2), np.nanmax(m2)))
bad = sum(int(((~np.isfinite(d[v].values)) | (np.abs(d[v].values) > 1e30))[:, m].sum()) for v in d.data_vars)
print("fill/NaN values on ocean cells (all variables, all waves): %d  (must be 0)" % bad)
sys.exit(1 if (bad or d.sizes["tide_period"] != 10) else 0)
PYEOF
cp "${TGEN}/croco_frc.nc" "${CF}/croco_frc.nc"
say "croco_frc.nc -> ${CF}; next: steps/07_compile.sh ${CONFIG_NAME} plain_tides (or plain_tides_rivers)"
