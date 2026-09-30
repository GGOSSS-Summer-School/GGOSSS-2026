#!/bin/bash
# =============================================================================
# Step 03 -- ocean + atmosphere data for the WHOLE period
#   hindcast/steps/03_download_data.sh NAME [--ocean-only | --atmos-only]
#
#   GLORYS (CMEMS)  daily reanalysis P1D-m, one file YYYY_MM.nc per month, box =
#                   grid + 1.5 deg, months GLORYS_MONTH_START..END (period plus the
#                   spin-up and the +/-1 day boundary pad -- derived, see lib/common.sh)
#   ERA5 (CDS)      10 variables, hourly, grid box (+2 deg added by the downloader),
#                   converted to CROCO online files for_croco/<VAR>_YyyyyMmm.nc
# Prerequisites (once): `copernicusmarine login`; ~/.cdsapirc + ERA5 licences.
# Re-running is safe: existing months / files are skipped. It can take hours
# (CDS queue, bandwidth); the cycling driver would also fetch missing months.
# Tunables (environment): CMEMS_STALL_MIN=10, CMEMS_MAX_RETRIES=5 (GLORYS pieces), CMEMS_SERVICE= (empty =
#   arco-time-series), ERA5_ARCO_PARALLEL=3, ERA5_N_PARALLEL=2 (CDS path), ERA5_LAG_DAYS=6.
#   ERA5 (default ERA5_SOURCE=arco): era5_arco.py reads the CDS ARCO Zarr store directly, no
#   queue. ERA5_SOURCE=cds: era5_fast.py (CDS API, single levels). Old path: --legacy.
# =============================================================================
#
# NESTED child (PARENT=<parent> in domain.cfg): nothing is downloaded. Instead
#   ocean  the parent's daily-mean outputs (spinup/ and hcast/ croco_avg.nc of every
#          parent cycle overlapping the child period) are converted to GLORYS/Mercator
#          format -- sigma levels -> the 50 standard z-levels, currents rotated to
#          east/north, names thetao/so/uo/vo/zos -- and written as monthly YYYY_MM.nc
#          in downloaded_data/PARENT_<parent> (gtools/nesting.py: parent_to_monthly).
#   ERA5   the parent's for_croco files are used as they are (checked here).
#   The parent hindcast must have been run first (hindcast/run_hindcast_cycle.sh <parent>).
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
WHAT=${2:-all}
cd "${GTOOLS_DIR}"

if [ -n "${PARENT}" ]; then
    # parent build folder suffix and cycle length, from the parent's domain.cfg
    read P_TIDES P_RIVERS P_HCAST <<< "$(bash -c "TIDES=0; RIVERS=0; HCAST_DAYS=5; source '${CROCO_CONFIGS_ROOT}/${PARENT}/domain.cfg'; echo \$TIDES \$RIVERS \$HCAST_DAYS")"
    P_BUILD=plain; [ "${P_TIDES}" = 1 ] && P_BUILD=${P_BUILD}_tides; [ "${P_RIVERS}" = 1 ] && P_BUILD=${P_BUILD}_rivers
    SUFFIX=""; [ "${P_BUILD}" != plain ] && SUFFIX="_${P_BUILD}"
    if [ -n "${PARENT_RUN}" ]; then
        RUNS=("${PARENT_RUN}")
    else   # every parent cycle whose window [T-spinup, T+HCAST_DAYS] overlaps the child's
        RUNS=()
        c0=$(date -u -d "${SPIN_FIRST} - 1 day" +%s); c1=$(date -u -d "${HC_END} + 1 day" +%s)
        for d in "${PARENT_OUTPUT_ROOT}"/[0-9]*"${SUFFIX}"; do
            tag=$(basename "$d"); tag=${tag%%_*}; [ "${#tag}" = 8 ] || continue
            [ -n "${SUFFIX}" ] || [ "$(basename "$d")" = "${tag}" ] || continue
            t0=$(date -u -d "${tag} - 3 days" +%s); t1=$(date -u -d "${tag} + ${P_HCAST} days" +%s)
            [ "${t1}" -ge "${c0}" ] && [ "${t0}" -le "${c1}" ] && RUNS+=("$d")
        done
    fi
    [ "${#RUNS[@]}" -gt 0 ] || die "no ${PARENT} run (${P_BUILD}) covers ${SPIN_FIRST}..${HC_END} in ${PARENT_OUTPUT_ROOT} -- run hindcast/run_hindcast_cycle.sh ${PARENT} first"
    FILES=()
    for d in "${RUNS[@]}"; do
        grep -qs "MAIN: DONE" "$d"/hcast/croco_hcast.out || die "parent run not finished: $d/hcast/croco_hcast.out has no MAIN: DONE"
        for ph in spinup hcast; do
            [ -f "$d/$ph/CROCO_FILES/croco_avg.nc" ] && FILES+=("$d/$ph/CROCO_FILES/croco_avg.nc")
        done
    done
    say "nesting in ${PARENT} (${P_BUILD}): converting ${#FILES[@]} daily-mean file(s) to GLORYS format"
    printf '    %s\n' "${FILES[@]}"
    rm -f "${GLORYS_DIR}"/20*.nc
    python -c "
import sys; sys.path.insert(0, '.')
import nesting
nesting.parent_to_monthly(sys.argv[2:], sys.argv[1], Yorig=${YORIG})
" "${GLORYS_DIR}" "${FILES[@]}"
    say "ERA5: the parent's files in ${ERA5_DIR}/for_croco"
    for ym in $(seq 0 12); do
        m=$(date -u -d "${ERA5_MONTH_START}-01 + ${ym} months" +%Y-%m)
        [[ "$m" > "${ERA5_MONTH_END}" ]] && break
        y=${m%-*}; mm=${m#*-}
        [ -s "${ERA5_DIR}/for_croco/T2M_Y${y}M${mm}.nc" ] || die "parent ERA5 missing for ${m} -- hindcast/steps/03_download_data.sh ${PARENT} --atmos-only"
        echo "  ${m}: $(ls "${ERA5_DIR}"/for_croco/*_Y${y}M${mm}.nc | wc -l) files"
    done
    say "next: hindcast/steps/04_make_ini_bry.sh ${CONFIG_NAME}"
    exit 0
fi

if [ "${WHAT}" != "--atmos-only" ]; then
    say "GLORYS ${GLORYS_MONTH_START} .. ${GLORYS_MONTH_END}, box ${EXTENTS}"
    python ggosss26.py download_ocean_hindcast --domain="${EXTENTS}" \
        --month_start "${GLORYS_MONTH_START}" --month_end "${GLORYS_MONTH_END}" \
        --product_id "${GLORYS_PRODUCT}" --source "${OCEAN_SOURCE}" --outputDir "${GLORYS_DIR}"
    say "CHECK -- each month: time = 30/31 (daily), depth = 50"
    for f in "${GLORYS_DIR}"/*.nc; do
        echo "  $(basename "$f"): $(ncdump -h "$f" | grep -E '^\s*(time|depth) = ' | tr -s ' \t\n' ' ')"
    done
fi

if [ "${WHAT}" != "--ocean-only" ]; then
    say "ERA5 ${ERA5_MONTH_START} .. ${ERA5_MONTH_END}, box ${ERA5_BOX} (+2 deg)"
    python ggosss26.py download_atmosphere_hindcast --domain="${ERA5_BOX}" \
        --month_start "${ERA5_MONTH_START}" --month_end "${ERA5_MONTH_END}" \
        --outputDir "${ERA5_DIR}" --Yorig "${YORIG}"
    NM=$(( ( $(date -u -d "${ERA5_MONTH_END}-01" +%Y)*12 + 10#$(date -u -d "${ERA5_MONTH_END}-01" +%m) ) \
         - ( $(date -u -d "${ERA5_MONTH_START}-01" +%Y)*12 + 10#$(date -u -d "${ERA5_MONTH_START}-01" +%m) ) + 1 ))
    say "CHECK -- expect $(( NM * 10 )) files (10 per month: T2M Q TP SSR STRD U10M V10M msl SST LSM):"
    ls "${ERA5_DIR}"/for_croco/*.nc | wc -l
fi
say "next: hindcast/steps/04_make_ini_bry.sh ${CONFIG_NAME}"
