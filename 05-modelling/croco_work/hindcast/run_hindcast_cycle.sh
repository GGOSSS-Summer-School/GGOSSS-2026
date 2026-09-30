#!/bin/bash
# =============================================================================
# run_hindcast_cycle.sh -- the operational cycling hindcast
#
#   hindcast/run_hindcast_cycle.sh NAME [--tides|--no-tides] [--rivers|--no-rivers]
#                                       [--start YYYY-MM-DD] [--ncycles N]
#
# Everything comes from hindcast/configs/NAME/domain.cfg (via lib/common.sh):
# period HC_START..HC_END -> NCYCLES cycles; boxes; Yorig; DT; the physics.
# Each cycle T (stride HCAST_DAYS, so the hindcast windows tile the period):
#   [0] --tides : croco_frc.nc built for this cycle at the spin-up start (TPXO)
#   [1] spin-up ini (GLORYS at T-SPINUP) + bry (T-SPINUP-1 .. T+1)
#   [2] spin-up run SPINUP_DAYS -> croco_rst.nc
#   [3] hindcast bry (T-1 .. T+HCAST_DAYS+1); IC = the spin-up restart (NRREC 1)
#   [4] hindcast run HCAST_DAYS -> croco_his.nc, croco_avg.nc
# Every cycle restarts from GLORYS through its own spin-up (anchored to the
# reanalysis, no drift). Missing GLORYS/ERA5 months are downloaded on the fly.
# The online (ERA5) block spans months/years automatically (e.g. 2018 05 24 2018 06).
#
# Physics (compile-time => one binary per combination, built by build_config.sh
# or steps/07_compile.sh NAME <build>); default = TIDES/RIVERS in domain.cfg:
#   build plain | plain_tides | plain_rivers | plain_tides_rivers -> croco_<build>
#   --rivers stages the pre-built CROCO_FILES/croco_runoff.nc (step 10) in both phases
#   --tides  retimes the output to hourly history / daily averages
# Output: hindcast/model-runs/NAME/<YYYYMMDD>/ (plain) or <YYYYMMDD>_<build>/
#   spinup/ hcast/ gen_spinup/ gen_hcast/ [gen_tides/]
# Run it in the background:  nohup hindcast/run_hindcast_cycle.sh NAME > NAME.log 2>&1 &
# Re-running skips cycles whose hcast already ended with MAIN: DONE.
# =============================================================================
set -e -o pipefail
HC_HOME="$(cd "$(dirname "$0")" && pwd)"
[ -n "$1" ] && [[ "$1" != --* ]] || { sed -n '2,30p' "$0"; exit 1; }
NAME=$1; shift
source "${HC_HOME}/lib/common.sh" "${NAME}"

USE_TIDES=${TIDES}; USE_RIVERS=${RIVERS}; START_DATE=${HC_START}; NCYC=${NCYCLES}
while [ $# -gt 0 ]; do
    case "$1" in
        --tides)     USE_TIDES=1; shift ;;
        --no-tides)  USE_TIDES=0; shift ;;
        --rivers)    USE_RIVERS=1; shift ;;
        --no-rivers) USE_RIVERS=0; shift ;;
        --start)     START_DATE=$2; shift 2 ;;
        --ncycles)   NCYC=$2; shift 2 ;;
        -h|--help)   sed -n '2,30p' "$0"; exit 0 ;;
        *) die "unknown option $1" ;;
    esac
done
BUILD="plain"; [ "${USE_TIDES}" = 1 ] && BUILD="${BUILD}_tides"; [ "${USE_RIVERS}" = 1 ] && BUILD="${BUILD}_rivers"

# ---- guards ---------------------------------------------------------------------
RUN_ROOT="${HCAST}"
CROCO_BIN="${RUN_ROOT}/croco_${BUILD}"
[ -x "${CROCO_BIN}" ] || die "binary not found: ${CROCO_BIN} -- build it: hindcast/steps/07_compile.sh ${NAME} ${BUILD}"
[ -f "${CF}/croco_grd.nc" ] || die "grid not found -- hindcast/build_config.sh ${NAME}"
[ -f "${CONFIG_DIR}/croco.in" ] && [ -f "${CONFIG_DIR}/crocotools_param.py" ] || die "recipe incomplete -- run steps 02, 05, 06"
if [ "${USE_TIDES}" = 1 ]; then
    [ -f "${CONFIG_DIR}/crocotools_param_tides.py" ] || die "crocotools_param_tides.py missing -- step 02"
    [ -f "${TPXO_DIR}/TPXO7.nc" ] || die "TPXO not found: ${TPXO_DIR}/TPXO7.nc"
fi
if [ "${USE_RIVERS}" = 1 ]; then
    [ -f "${CF}/croco_runoff.nc" ] || die "CROCO_FILES/croco_runoff.nc missing -- hindcast/steps/10_rivers_glofas.sh ${NAME}"
fi

ERA5_CROCO_DIR="${ERA5_DIR}/for_croco"
NRST=$(( 86400 / DT )); NWRT=$(( 21600 / DT )); NAVG=$(( 21600 / DT ))
[ "${USE_TIDES}" = 1 ] && { NWRT=$(( 3600 / DT )); NAVG=$(( 86400 / DT )); }
add_days () { date -u -d "$1 + $2 days" +%Y-%m-%d; }

# ---- patch_croco_in <file> <sdate> <edate> <days> ------------------------------
patch_croco_in () {
    local f="$1" sdate="$2" edate="$3" days="$4"
    local ntimes=$(( days * NTIMES_PER_DAY ))
    local by=$(date -u -d "${sdate}" +%Y) bm=$(date -u -d "${sdate}" +%m)
    local bye=$(date -u -d "${edate}" +%Y) bme=$(date -u -d "${edate}" +%m)
    sed -i "/^start_date:/{n; s/.*/${sdate} 00:00:00/}" "$f"
    sed -i "/^end_date:/{n; s/.*/${edate} 00:00:00/}"   "$f"
    sed -i "/^time_stepping:/{n; s/.*/                ${ntimes}      ${DT}      ${NDTFAST}      1/}" "$f"
    sed -i "/^restart:/{n; s/.*/                   ${NRST}    -1/}" "$f"
    sed -i "/^history:/{n; s/.*/            T      ${NWRT}     0/}" "$f"
    sed -i "/^averages:/{n; s/.*/            1      ${NAVG}     0/}" "$f"
    sed -i "/^initial:/{n; s/.*/          1/}" "$f"
    sed -i "/^initial:/{n; n; s|.*|    CROCO_FILES/croco_ini.nc|}" "$f"
    sed -i "/^boundary:/{n; s|.*|    CROCO_FILES/croco_bry.nc|}" "$f"
    sed -i '/^forcing:/{n; s|.*|    CROCO_FILES/croco_frc.nc|}' "$f"
    sed -i "/^online:/{n; s/.*/           ${by}   ${bm}      24            ${bye}     ${bme}/}" "$f"
    sed -i "/^online:/{n; n; s|.*|    ${ERA5_CROCO_DIR}/|}" "$f"
    if [ "${USE_RIVERS}" = 1 ]; then
        grep -A1 "^psource_ncfile:" "$f" | grep -q "CROCO_FILES/croco_runoff.nc" || write_psource_block "$f"
    fi
    return 0
}

# ---- stage <run dir>: binary, recipe, grid, optional physics ---------------------
stage () {
    local d="$1"
    cp "${CONFIG_DIR}"/{cppdefs.h,param.h,jobcomp,croco.in} "$d/"
    cp "${CROCO_BIN}" "$d/croco"
    cp "${CF}/croco_grd.nc" "$d/CROCO_FILES/croco_grd.nc"
    [ "${USE_TIDES}"  = 1 ] && cp "${TIDE_FRC}" "$d/CROCO_FILES/croco_frc.nc"
    [ "${USE_RIVERS}" = 1 ] && cp "${CF}/croco_runoff.nc" "$d/CROCO_FILES/croco_runoff.nc"
    return 0
}

# ---- ensure_data <date0> <date1>: GLORYS + ERA5 for every month spanned ---------
ensure_data () {
    local d="$(date -u -d "$1" +%Y-%m-01)" last="$(date -u -d "$2" +%Y-%m-01)"
    cd "${GTOOLS_DIR}"
    while [ "$(date -u -d "$d" +%s)" -le "$(date -u -d "$last" +%s)" ]; do
        local ym=$(date -u -d "$d" +%Y-%m); local y=${ym%-*} m=${ym#*-}
        if [ ! -s "${GLORYS_DIR}/${y}_${m}.nc" ] && [ -n "${PARENT}" ]; then
            die "nested in ${PARENT}: no parent data for ${ym} in ${GLORYS_DIR} -- run hindcast/steps/03_download_data.sh ${CONFIG_NAME} (after the ${PARENT} hindcast)"
        fi
        if [ ! -s "${ERA5_CROCO_DIR}/T2M_Y${y}M${m}.nc" ] && [ -n "${PARENT}" ]; then
            die "nested in ${PARENT}: its ERA5 has no ${ym} (${ERA5_CROCO_DIR})"
        fi
        if [ ! -s "${GLORYS_DIR}/${y}_${m}.nc" ]; then
            say "downloading GLORYS ${ym} ..."
            python ggosss26.py download_ocean_hindcast --domain="${EXTENTS}" --month_start "${ym}" \
                --month_end "${ym}" --product_id "${GLORYS_PRODUCT}" --source "${OCEAN_SOURCE}" --outputDir "${GLORYS_DIR}"
        fi
        if [ ! -s "${ERA5_CROCO_DIR}/T2M_Y${y}M${m}.nc" ]; then
            say "downloading + converting ERA5 ${ym} ..."
            python ggosss26.py download_atmosphere_hindcast --domain="${ERA5_BOX}" --month_start "${ym}" \
                --month_end "${ym}" --outputDir "${ERA5_DIR}" --Yorig "${YORIG}"
        fi
        d=$(add_days "$d" 32); d=$(date -u -d "$d" +%Y-%m-01)
    done
}

show_config
echo " driver: build ${BUILD} (binary ${CROCO_BIN##*/}), ${NCYC} cycle(s) from ${START_DATE}"

T="${START_DATE}"
for (( c=1; c<=NCYC; c++ )); do
    TAG=$(date -u -d "$T" +%Y%m%d)
    SPIN_START=$(add_days "$T" -${SPINUP_DAYS}); SPIN_END="$T"
    HC_S="$T"; HC_E=$(add_days "$T" ${HCAST_DAYS})
    if [ "${BUILD}" = plain ]; then CYCLE_ROOT="${OUTPUT_ROOT}/${TAG}"; else CYCLE_ROOT="${OUTPUT_ROOT}/${TAG}_${BUILD}"; fi
    SPIN_DIR="${CYCLE_ROOT}/spinup"; HC_DIR="${CYCLE_ROOT}/hcast"
    SPIN_GEN="${CYCLE_ROOT}/gen_spinup/CROCO_FILES"; HC_GEN="${CYCLE_ROOT}/gen_hcast/CROCO_FILES"

    echo ""
    echo "############################################################"
    echo "# CYCLE ${c}/${NCYC}  T=${TAG}  (${BUILD})"
    echo "#   spin-up : ${SPIN_START} -> ${SPIN_END}"
    echo "#   hindcast: ${HC_S} -> ${HC_E}"
    echo "############################################################"
    if grep -qs "MAIN: DONE" "${HC_DIR}/croco_hcast.out"; then
        echo "  already done (${HC_DIR}/croco_hcast.out) -- skipping"; T="${HC_E}"; continue
    fi
    mkdir -p "${SPIN_DIR}/CROCO_FILES" "${HC_DIR}/CROCO_FILES" "${SPIN_GEN}" "${HC_GEN}"
    ensure_data "$(add_days "${SPIN_START}" -2)" "$(add_days "${HC_E}" 2)"

    if [ "${USE_TIDES}" = 1 ]; then
        say "[0/4] tide file (TPXO) at ${SPIN_START} ..."
        TGEN="${CYCLE_ROOT}/gen_tides/CROCO_FILES"; mkdir -p "${TGEN}"
        cp "${CF}/croco_grd.nc" "${TGEN}/"; cp "${CONFIG_DIR}/crocotools_param_tides.py" "${TGEN}/crocotools_param.py"
        cd "${GTOOLS_DIR}"
        python ggosss26.py make_tides --input_dir "${TPXO_DIR}" --output_dir "${TGEN}" \
            --run_date "${SPIN_START} 00:00:00" --Yorig "${YORIG}" --fname_out croco_frc.nc > "${TGEN}/make_tides.log" 2>&1 \
            || die "make_tides failed -- ${TGEN}/make_tides.log"
        TIDE_FRC="${TGEN}/croco_frc.nc"
    fi

    say "[1/4] spin-up ini + bry (${OCEAN_TAG}) ..."
    cp "${CF}/croco_grd.nc" "${CONFIG_DIR}/crocotools_param.py" "${SPIN_GEN}/"
    cd "${GTOOLS_DIR}"
    python ggosss26.py make_ini_hindcast --input_dir "${GLORYS_DIR}" --output_dir "${SPIN_GEN}" \
        --date "${SPIN_START}" --Yorig "${YORIG}" > "${SPIN_GEN}/make_ini.log" 2>&1 || die "make_ini_hindcast failed -- ${SPIN_GEN}/make_ini.log"
    python ggosss26.py make_bry_hindcast --input_dir "${GLORYS_DIR}" --output_dir "${SPIN_GEN}" \
        --start_date "$(add_days "${SPIN_START}" -1)" --end_date "$(add_days "${SPIN_END}" 1)" --Yorig "${YORIG}" \
        > "${SPIN_GEN}/make_bry.log" 2>&1 || die "make_bry_hindcast failed -- ${SPIN_GEN}/make_bry.log"
    SPIN_INI=$(ls -t "${SPIN_GEN}"/croco_ini_*.nc 2>/dev/null | head -1)
    SPIN_BRY=$(ls -t "${SPIN_GEN}"/croco_bry_*.nc 2>/dev/null | head -1)
    [ -n "${SPIN_INI}" ] && [ -n "${SPIN_BRY}" ] || die "spin-up ini/bry not produced"
    stage "${SPIN_DIR}"
    cp -L "${SPIN_INI}" "${SPIN_DIR}/CROCO_FILES/croco_ini.nc"
    cp -L "${SPIN_BRY}" "${SPIN_DIR}/CROCO_FILES/croco_bry.nc"
    patch_croco_in "${SPIN_DIR}/croco.in" "${SPIN_START}" "${SPIN_END}" "${SPINUP_DAYS}"

    say "[2/4] spin-up run (${SPINUP_DAYS} d) ..."
    ( cd "${SPIN_DIR}" && ./croco croco.in > croco_spinup.out 2>&1 ) || true
    grep -q "MAIN: DONE" "${SPIN_DIR}/croco_spinup.out" || { tail -20 "${SPIN_DIR}/croco_spinup.out"; die "spin-up did not reach MAIN: DONE -- ${SPIN_DIR}/croco_spinup.out"; }
    SPIN_RST="${SPIN_DIR}/CROCO_FILES/croco_rst.nc"
    [ -e "${SPIN_RST}" ] || die "spin-up produced no restart"

    say "[3/4] hindcast bry (${OCEAN_TAG}); IC = spin-up restart ..."
    cp "${CF}/croco_grd.nc" "${CONFIG_DIR}/crocotools_param.py" "${HC_GEN}/"
    cd "${GTOOLS_DIR}"
    python ggosss26.py make_bry_hindcast --input_dir "${GLORYS_DIR}" --output_dir "${HC_GEN}" \
        --start_date "$(add_days "${HC_S}" -1)" --end_date "$(add_days "${HC_E}" 1)" --Yorig "${YORIG}" \
        > "${HC_GEN}/make_bry.log" 2>&1 || die "make_bry_hindcast failed -- ${HC_GEN}/make_bry.log"
    HC_BRY=$(ls -t "${HC_GEN}"/croco_bry_*.nc 2>/dev/null | head -1)
    [ -n "${HC_BRY}" ] || die "hindcast bry not produced"
    stage "${HC_DIR}"
    cp "${SPIN_RST}" "${HC_DIR}/CROCO_FILES/croco_ini.nc"
    cp -L "${HC_BRY}" "${HC_DIR}/CROCO_FILES/croco_bry.nc"
    patch_croco_in "${HC_DIR}/croco.in" "${HC_S}" "${HC_E}" "${HCAST_DAYS}"

    say "[4/4] hindcast run (${HCAST_DAYS} d) ..."
    ( cd "${HC_DIR}" && ./croco croco.in > croco_hcast.out 2>&1 ) || true
    grep -q "MAIN: DONE" "${HC_DIR}/croco_hcast.out" || { tail -20 "${HC_DIR}/croco_hcast.out"; die "hindcast did not reach MAIN: DONE -- ${HC_DIR}/croco_hcast.out"; }
    echo "  cycle ${TAG} done -> ${HC_DIR}/CROCO_FILES/croco_his.nc"
    T="${HC_E}"
done

echo ""
echo "============================================================"
echo " HINDCAST DONE  ${CONFIG_NAME}  build ${BUILD}  (${NCYC} cycles from ${START_DATE})"
echo "   outputs: ${OUTPUT_ROOT}/<YYYYMMDD>$( [ "${BUILD}" = plain ] || echo "_${BUILD}")/hcast/CROCO_FILES/croco_his.nc"
echo "   check:   python3 ${CROCO_ROOT}/notebooks/verify_config.py ${CONFIG_NAME} ${BUILD}"
echo "============================================================"
