#!/bin/bash
# =============================================================================
# Step 06 -- croco.in, the run-time file (incl. the optional tide/river entries)
#   hindcast/steps/06_patch_croco_in.sh NAME
#
# Writes the recipe croco.in (hindcast/configs/NAME/croco.in) for the proof run,
# rewriting the line(s) BELOW each section header (sed '/^header:/{n; s/.*/.../}'):
#   title          <CONFIG_CPP> HINDCAST
#   time_stepping  NTIMES = TEST_DAYS*86400/DT, DT, NDTFAST, 1
#   S-coord        THETA_S THETA_B SIGMA_HC (= sigma_params = param.h N)
#   initial        NRREC 1 + CROCO_FILES/croco_ini_<TAG>_<TEST_START>.nc
#   boundary       CROCO_FILES/croco_bry_<TAG>_<TEST_START-1>_to_<TEST_END+1>.nc
#   forcing        CROCO_FILES/croco_frc.nc (the tide file; unused without TIDES)
#   sponge         SPONGE (default 0. 0.)
#   online         ERA5 form: byear bmonth 24 byearend bmonthend + for_croco/ path
#                  (bmonthend = the LAST month of the run)
#   psource_ncfile rivers, only if CROCO_FILES/croco_runoff.nc exists (step 10)
# The run/driver retime history/averages per build (6 h; 1 h/24 h with tides).
# start_date/end_date are left alone: USE_CALENDAR is off (harmless warning).
# Finally stages cppdefs.h, param.h, croco.in, jobcomp into the workbench HCAST.
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
F="${CONFIG_DIR}/croco.in"
[ -f "${CONFIG_DIR}/cppdefs.h" ] || die "run hindcast/steps/05_patch_source_files.sh ${CONFIG_NAME} first"
cp "${CROCO_MODEL_DIR}/OCEAN/croco.in" "$F"          # always start from the template

NTIMES=$(( TEST_DAYS * NTIMES_PER_DAY ))
BRY_START=$(date -u -d "${TEST_START} - 1 day" +%Y-%m-%d)
BRY_END=$(date -u -d "${TEST_END} + 1 day" +%Y-%m-%d)
INI=croco_ini_${OCEAN_TAG}_$(date -u -d "${TEST_START}" +Y%YM%mD%d).nc
BRY=croco_bry_${OCEAN_TAG}_$(date -u -d "${BRY_START}" +Y%YM%mD%d)_to_$(date -u -d "${BRY_END}" +Y%YM%mD%d).nc
BY=$(date -u -d "${TEST_START}" +%Y); BM=$(date -u -d "${TEST_START}" +%-m)
BYE=$(date -u -d "${TEST_END}" +%Y);  BME=$(date -u -d "${TEST_END}" +%-m)
read SP1 SP2 <<< "${SPONGE}"

sed -i "/^title:/{n; s/.*/        ${CONFIG_CPP} HINDCAST/}" "$F"
sed -i "/^time_stepping:/{n; s/.*/                ${NTIMES}      ${DT}      ${NDTFAST}      1/}" "$F"
SC=$(awk -v a="${THETA_S}" -v b="${THETA_B}" -v c="${SIGMA_HC}" 'BEGIN{printf "           %.1fd0     %.1fd0      %.1fd0", a, b, c}')
sed -i "/^S-coord:/{n; s/.*/${SC}/}" "$F"
sed -i "/^initial:/{n; s/.*/          1/}" "$F"
sed -i "/^initial:/{n; n; s|.*|    CROCO_FILES/${INI}|}" "$F"
sed -i "/^boundary:/{n; s|.*|    CROCO_FILES/${BRY}|}" "$F"
sed -i '/^forcing:/{n; s|.*|    CROCO_FILES/croco_frc.nc|}' "$F"
sed -i "/^sponge:/{n; s/.*/                    ${SP1}                ${SP2}/}" "$F"
sed -i "/^online:/{n; s/.*/           ${BY}   ${BM}      24            ${BYE}     ${BME}/}" "$F"
sed -i "/^online:/{n; n; s|.*|    ${ERA5_DIR}/for_croco/|}" "$F"
if [ -f "${CF}/croco_runoff.nc" ] && [ -f "${CF}/for_croco_in.txt" ]; then
    write_psource_block "$F"
fi

say "CHECK croco.in:"
sed -n '2p' "$F"
grep -A1 -E "^(time_stepping|S-coord|boundary|forcing|sponge):" "$F" | grep -v "^--"
grep -A2 -E "^(initial|online):" "$F" | grep -v "^--"
grep -n "XXX" "$F" && die "croco.in still has XXX" || echo "  no XXX left"
[ -f "${CF}/${INI}" ] && [ -f "${CF}/${BRY}" ] && echo "  ini/bry files exist" \
    || echo "  NOTE: ${INI} / ${BRY} not built yet (step 04)"

say "staging recipe -> workbench ${HCAST}"
cp "${CONFIG_DIR}"/{cppdefs.h,param.h,croco.in,jobcomp} "${HCAST}/"
say "next: hindcast/steps/07_compile.sh ${CONFIG_NAME} [build]"
