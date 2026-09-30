#!/bin/bash
# =============================================================================
# Step 08 -- the proof run (any build)
#   hindcast/steps/08_run_test.sh NAME [BUILD]
#
# Runs ./croco_<BUILD> for TEST_DAYS days from TEST_START with the recipe croco.in.
# Output retimed per build: history/averages every 6 h; with tides hourly history
# (resolves M2) and daily averages (removes it). Watch in the log: GET_INITIAL,
# GET_BRY, ONLINE_BULK -- Read file, (GET_PSOURCE with rivers), a bounded kinetic
# energy, and the first time[DAYS] = days since YORIG-01-01. Must end MAIN: DONE.
# Outputs are copied to hindcast/scratch/NAME/test_<BUILD>/ (runs never overwrite).
# Checks printed: time span, surface T/S; with tides the sea-level period (~12.4 h).
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
BUILD=${2:-${DEFAULT_BUILD}}
build_flags "${BUILD}"
cd "${HCAST}"
[ -x "./croco_${BUILD}" ] || die "./croco_${BUILD} missing -- run steps/07_compile.sh ${CONFIG_NAME} ${BUILD}"
[ "${B_TIDES}" = 1 ]  && { [ -f CROCO_FILES/croco_frc.nc ]    || die "CROCO_FILES/croco_frc.nc missing -- steps/09_tides.sh ${CONFIG_NAME}"; }
[ "${B_RIVERS}" = 1 ] && { [ -f CROCO_FILES/croco_runoff.nc ] || die "CROCO_FILES/croco_runoff.nc missing -- steps/10_rivers_glofas.sh ${CONFIG_NAME}"; }

cp "${CONFIG_DIR}/croco.in" croco.in
if [ "${B_RIVERS}" = 1 ]; then
    grep -A1 "^psource_ncfile:" croco.in | grep -q "CROCO_FILES/croco_runoff.nc" || write_psource_block croco.in
fi
if [ "${B_TIDES}" = 1 ]; then NWRT=$(( 3600 / DT )); NAVG=$(( 86400 / DT )); else NWRT=$(( 21600 / DT )); NAVG=$(( 21600 / DT )); fi
sed -i "/^history:/{n; s/.*/            T      ${NWRT}     0/}" croco.in
sed -i "/^averages:/{n; s/.*/            1      ${NAVG}     0/}" croco.in

say "build ${BUILD}, ${TEST_START} -> ${TEST_END}; expected first time[DAYS] = $(( ( $(date -u -d "${TEST_START}" +%s) - $(date -u -d "${YORIG}-01-01" +%s) ) / 86400 ))"
rm -f CROCO_FILES/croco_his.nc CROCO_FILES/croco_avg.nc CROCO_FILES/croco_rst.nc
"./croco_${BUILD}" croco.in > "run_${BUILD}.log" 2>&1 || true
tail -15 "run_${BUILD}.log"
if ! grep -q "MAIN: DONE" "run_${BUILD}.log"; then
    grep -nE "ERROR|BLOW|Abnormal|NaN" "run_${BUILD}.log" | head
    die "run did not reach MAIN: DONE -- see ${HCAST}/run_${BUILD}.log"
fi
mkdir -p "test_${BUILD}"
cp CROCO_FILES/croco_his.nc CROCO_FILES/croco_avg.nc "run_${BUILD}.log" "test_${BUILD}/"
say "OK: MAIN: DONE -> ${HCAST}/test_${BUILD}/"
grep -m2 "GET_PSOURCE" "run_${BUILD}.log" || true

python3 - <<PYEOF
import xarray as xr, numpy as np
d = xr.open_dataset("${HCAST}/test_${BUILD}/croco_his.nc", decode_times=False)
m = d.mask_rho.values > 0
t = d.scrum_time.values / 86400.
print("records %d, time %.3f .. %.3f days since ${YORIG}-01-01" % (len(t), t[0], t[-1]))
s = d.salt.isel(time=-1, s_rho=-1).values[m]; T = d.temp.isel(time=-1, s_rho=-1).values[m]
print("surface T %.2f..%.2f degC   surface S %.2f..%.2f (mean %.2f)" % (T.min(), T.max(), s.min(), s.max(), s.mean()))
if ${B_TIDES}:
    h = d.h.values; lon, lat = d.lon_rho.values, d.lat_rho.values
    shelf = m & (h > 30) & (h < 200)
    idx = np.argwhere(shelf)
    for j, i in (idx[np.linspace(0, len(idx) - 1, 3).astype(int)] if len(idx) else []):
        z = d.zeta.isel(eta_rho=j, xi_rho=i).values; zz = z - z.mean()
        up = np.where((zz[:-1] < 0) & (zz[1:] >= 0))[0]
        per = np.mean(np.diff(t[up])) * 24 if len(up) > 2 else float('nan')
        print("zeta at %.2fE %.2fN (h=%.0f m): range %.2f m, mean period %.1f h (M2 = 12.4 h)"
              % (lon[j, i], lat[j, i], h[j, i], z.max() - z.min(), per))
PYEOF
