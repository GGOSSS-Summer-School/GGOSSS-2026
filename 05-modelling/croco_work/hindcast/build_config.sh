#!/bin/bash
# =============================================================================
# build_config.sh -- build (and prove) a hindcast configuration in one command.
#
#   hindcast/build_config.sh NAME [--from N] [--only N] [--no-download] [--no-test]
#
#   Runs the steps in hindcast/steps/ in order, with everything derived from
#   hindcast/configs/NAME/domain.cfg (create it with new_config.sh):
#     01 grid          02 boundaries/params/card   03 GLORYS + ERA5 downloads
#     04 ini + bry     05 cppdefs.h/param.h/jobcomp 06 croco.in
#     07 compile plain 08 7-day proof run (plain)
#   and, only if switched on in domain.cfg:
#     TIDES=1  -> 09 tide file
#     RIVERS=1 -> 10 GloFAS rivers (runoff file + croco.in block)
#     then 07 + 08 again for the final build (plain_tides / plain_rivers / plain_tides_rivers)
#
#   --from N       start at step N (e.g. --from 4 after the downloads finished)
#   --only N       run only step N
#   --no-download  skip step 03 (data already there, or let the driver fetch it)
#   --no-test      skip the proof runs (08)
# Each step is logged to hindcast/scratch/NAME/build_logs/stepNN.log and the build
# stops at the first failing step (fix, then --from that step).
# After this: hindcast/run_hindcast_cycle.sh NAME   (the full cycled hindcast)
# =============================================================================
set -e -o pipefail
HC_HOME="$(cd "$(dirname "$0")" && pwd)"
[ -n "$1" ] && [[ "$1" != --* ]] || { sed -n '2,27p' "$0"; exit 1; }
NAME=$1; shift
FROM=1; ONLY=""; NODL=0; NOTEST=0
while [ $# -gt 0 ]; do
    case "$1" in
        --from) FROM=$2; shift 2 ;;
        --only) ONLY=$2; shift 2 ;;
        --no-download) NODL=1; shift ;;
        --no-test) NOTEST=1; shift ;;
        *) echo "!! unknown option $1"; exit 1 ;;
    esac
done
source "${HC_HOME}/lib/common.sh" "${NAME}"
show_config
LOGS="${HCAST}/build_logs"; mkdir -p "${LOGS}"
S="${HC_HOME}/steps"

run () {   # run <step number> <script> [args...]
    local n=$1; shift
    if [ -n "${ONLY}" ]; then [ "$n" -eq "${ONLY}" ] || return 0
    else [ "$n" -ge "${FROM}" ] || return 0; fi
    local tag="step$(printf %02d "$n")$( [ -n "$3" ] && echo "_$3" )"
    echo ""; echo "################  ${tag}: $(basename "$1") $2 $3  ################"
    if ! bash "$@" 2>&1 | tee "${LOGS}/${tag}.log"; then
        echo "!! ${tag} failed -- see ${LOGS}/${tag}.log; fix it, then: $0 ${NAME} --from ${n}"
        exit 1
    fi
}

run 1 "$S/01_make_grid.sh"        "${NAME}"
run 2 "$S/02_check_grid.sh"       "${NAME}"
[ "${NODL}" = 1 ] || run 3 "$S/03_download_data.sh" "${NAME}"
run 4 "$S/04_make_ini_bry.sh"     "${NAME}"
run 5 "$S/05_patch_source_files.sh" "${NAME}"
run 6 "$S/06_patch_croco_in.sh"   "${NAME}"
run 7 "$S/07_compile.sh"          "${NAME}" plain
[ "${NOTEST}" = 1 ] || run 8 "$S/08_run_test.sh" "${NAME}" plain
[ "${TIDES}" = 1 ]  && run 9  "$S/09_tides.sh"        "${NAME}"
[ "${RIVERS}" = 1 ] && run 10 "$S/10_rivers_glofas.sh" "${NAME}"
if [ "${DEFAULT_BUILD}" != plain ]; then
    run 7 "$S/07_compile.sh" "${NAME}" "${DEFAULT_BUILD}"
    [ "${NOTEST}" = 1 ] || run 8 "$S/08_run_test.sh" "${NAME}" "${DEFAULT_BUILD}"
fi

[ -n "${ONLY}" ] && { echo ">>> step ${ONLY} done (log ${LOGS})"; exit 0; }
echo ""
echo "============================================================"
echo " ${CONFIG_NAME} built. Binaries: $(ls "${HCAST}"/croco_plain* 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"
echo " region card: ${CONFIG_DIR}/CARD.md   logs: ${LOGS}"
echo " next: ${HC_HOME}/run_hindcast_cycle.sh ${CONFIG_NAME}      (build ${DEFAULT_BUILD})"
echo "============================================================"
