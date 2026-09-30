#!/bin/bash
# =============================================================================
# Step 10 -- OPTIONAL: river forcing, discharge from GloFAS
#   hindcast/steps/10_rivers_glofas.sh NAME
#
# Instead of a climatology (e.g. Dai & Trenberth), a hindcast uses the REAL daily
# discharge of its period from GloFAS. Built ONCE for the whole period
# (GLORYS_MONTH_START..END), real dates, Yorig, no cycle:
#   1 download_rivers_hindcast  GloFAS v4 consolidated daily discharge from EWDS
#                               (request schema of 29 Jul 2026: year/month/day,
#                               timespan=time_mean, average_river_discharge_in_the_last_24_hours)
#   2 make_rivers_hindcast      river mouths/names from the Dai file (positions only), the
#                               outlet = max mean GloFAS discharge within RIVER_RADIUS,
#                               rivers below RIVER_QMIN dropped, no double counting
#                               -> CROCO_FILES/GLOFAS_RIVERS/<river>.txt + river_list.txt
#   3 make_river_run.py         -> CROCO_FILES/croco_runoff.nc + for_croco_in.txt
#   4 croco.in                  the psource_ncfile block (T=RIVER_TEMP, S=RIVER_SALT)
# Rivers are then switched on at compile time (07_compile.sh NAME plain_rivers or
# plain_tides_rivers: PSOURCE + PSOURCE_NCFILE, PSOURCE_NCFILE_TS stays undef).
# Prerequisite (once): ~/.ewdsapirc (EWDS url + ECMWF token) + GloFAS licence accepted.
# Not needed for a river-free hindcast (RIVERS=0).
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
RIVDIR=${CF}/GLOFAS_RIVERS
cd "${GTOOLS_DIR}"

say "[1] GloFAS ${GLORYS_MONTH_START} .. ${GLORYS_MONTH_END} ..."
python ggosss26.py download_rivers_hindcast --grid "${CF}/croco_grd.nc" \
    --month_start "${GLORYS_MONTH_START}" --month_end "${GLORYS_MONTH_END}" --outputDir "${GLOFAS_DIR}"

say "[2] river mouths + daily discharge ..."
rm -rf "${RIVDIR}"
EXTRA=(); [ -n "${RIVER_EXTRA}" ] && EXTRA=(--extra_list "${RIVER_EXTRA}")
python ggosss26.py make_rivers_hindcast --grid "${CF}/croco_grd.nc" --glofas_dir "${GLOFAS_DIR}" \
    --outputDir "${RIVDIR}" --month_start "${GLORYS_MONTH_START}" --month_end "${GLORYS_MONTH_END}" \
    --qmin "${RIVER_QMIN}" --radius "${RIVER_RADIUS}" --margin "${RIVER_MARGIN}" "${EXTRA[@]}"

say "[3] croco_runoff.nc (make_river_run.py, real dates: --rivers_cyl 0) ..."
YS=${GLORYS_MONTH_START%-*}; MS=$((10#${GLORYS_MONTH_START#*-}))
YE=${GLORYS_MONTH_END%-*};   ME=$((10#${GLORYS_MONTH_END#*-}))
rm -f "${CF}/croco_runoff.nc" "${CF}/for_croco_in.txt"
cd "${GTOOLS_DIR}/croco_pytools/prepro"          # the builder finds ./Modules from here
python make_river_run.py --input_dir "${RIVDIR}/" --croco_dir "${CF}/" \
    --Yorig "${YORIG}" --Ystart "${YS}" --Mstart "${MS}" --Yend "${YE}" --Mend "${ME}" \
    --rivers_cyl 0 --file_format FULL > "${RIVDIR}/make_river_run.log" 2>&1 \
    || { tail -15 "${RIVDIR}/make_river_run.log"; die "make_river_run.py failed"; }
echo "  $(grep -c "positionned in sea" "${RIVDIR}/make_river_run.log" || true) river(s) positioned on wet cells"

say "CHECK croco_runoff.nc:"
python3 - <<PYEOF
import xarray as xr, numpy as np
d = xr.open_dataset("${CF}/croco_runoff.nc", decode_times=False)
t = d.qbar_time.values; o = np.datetime64("${YORIG}-01-01")
f = lambda x: str(o + np.timedelta64(int(round(x * 86400)), 's'))[:16]
print("n_qbar %d, qbar_time %s -> %s (%d daily records, units '%s', cycle_length %s)"
      % (d.sizes["n_qbar"], f(t.min()), f(t.max()), len(t), d.qbar_time.attrs.get("units"),
         d.qbar_time.attrs.get("cycle_length", "none")))
for k, n in enumerate(d.runoff_name.values):
    q = d.Qbar[k].values; print("  %-16s Q mean %7.0f  min %7.0f  max %7.0f m3/s" % (str(n).strip(), q.mean(), q.min(), q.max()))
PYEOF

say "[4] croco.in psource_ncfile block (recipe):"
[ -f "${CONFIG_DIR}/croco.in" ] || die "no recipe croco.in -- run steps 05-06 first"
write_psource_block "${CONFIG_DIR}/croco.in"
cp "${CONFIG_DIR}/croco.in" "${HCAST}/croco.in"
grep -A$(( $(ncdump -h "${CF}/croco_runoff.nc" | awk '/n_qbar =/{print $3}') + 2 )) "^psource_ncfile:" "${CONFIG_DIR}/croco.in"
say "next: steps/07_compile.sh ${CONFIG_NAME} plain_rivers (or plain_tides_rivers)"
