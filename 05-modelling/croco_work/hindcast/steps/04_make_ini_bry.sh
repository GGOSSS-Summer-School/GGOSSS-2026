#!/bin/bash
# =============================================================================
# Step 04 -- initial + boundary conditions from the ocean source (GLORYS or Mercator) for the proof run
#            The cycling driver builds its own per cycle.
#   hindcast/steps/04_make_ini_bry.sh NAME
#
#   ini  at TEST_START (day 2 of the period: GLORYS records on both sides)
#   bry  TEST_START-1 day .. TEST_END+1 day (CROCO needs a record bracketing
#        every step, else "ERROR in get_bry: cannot read variable 'bry_time'")
#   Only the open boundaries (obc_dict from step 02) get boundary data.
# Output: CROCO_FILES/croco_ini_<TAG>_Y..M..D...nc, croco_bry_<TAG>_Y.._to_Y...nc (TAG = GLORYS or MERCATOR)
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
need_obc
[ -f "${CONFIG_DIR}/crocotools_param.py" ] || die "run hindcast/steps/02_check_grid.sh ${CONFIG_NAME} first"
cp "${CONFIG_DIR}/crocotools_param.py" "${CF}/"
cd "${GTOOLS_DIR}"

say "initial condition (${OCEAN_SOURCE}, ${TEST_START}) ..."
python ggosss26.py make_ini_hindcast --input_dir "${GLORYS_DIR}" --output_dir "${CF}" \
    --date "${TEST_START}" --Yorig "${YORIG}"

BRY_START=$(date -u -d "${TEST_START} - 1 day" +%Y-%m-%d)
BRY_END=$(date -u -d "${TEST_END} + 1 day" +%Y-%m-%d)
say "boundaries (${OCEAN_SOURCE}, ${BRY_START} -> ${BRY_END}) ..."
python ggosss26.py make_bry_hindcast --input_dir "${GLORYS_DIR}" --output_dir "${CF}" \
    --start_date "${BRY_START}" --end_date "${BRY_END}" --Yorig "${YORIG}"

INI=${CF}/croco_ini_${OCEAN_TAG}_$(date -u -d "${TEST_START}" +Y%YM%mD%d).nc
BRY=${CF}/croco_bry_${OCEAN_TAG}_$(date -u -d "${BRY_START}" +Y%YM%mD%d)_to_$(date -u -d "${BRY_END}" +Y%YM%mD%d).nc
say "CHECK -- s_rho = ${N_LEVELS}, time units since ${YORIG}-01-01:"
ncdump -h "${INI}" | grep -E "s_rho = |scrum_time:units"
ncdump -h "${BRY}" | grep -E "bry_time = |bry_time:units"
say "CHECK -- boundary variables only on the open edges (obc: S=${OBC_SOUTH} W=${OBC_WEST} E=${OBC_EAST} N=${OBC_NORTH}):"
for e in south west east north; do printf "  %-5s %d variables\n" $e "$(ncdump -h "${BRY}" | grep -c "_${e}(")"; done
say "next: hindcast/steps/05_patch_source_files.sh ${CONFIG_NAME}"
