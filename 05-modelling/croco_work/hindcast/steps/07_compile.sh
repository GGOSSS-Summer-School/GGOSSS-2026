#!/bin/bash
# =============================================================================
# Step 07 -- compile a CROCO binary (plain, or with the optional tides/rivers)
#   hindcast/steps/07_compile.sh NAME [BUILD]
#     BUILD = plain | plain_tides | plain_rivers | plain_tides_rivers
#             (default: what TIDES/RIVERS in domain.cfg ask for)
#
# Tides and rivers are COMPILE-TIME options, so each combination is its own
# binary, kept as hindcast/scratch/NAME/croco_<BUILD> (the name the driver selects).
# The recipe files (step 05) are staged into the workbench, then the TIDES,
# PSOURCE and PSOURCE_NCFILE switches of the <CONFIG_CPP> block are set for BUILD
# (PSOURCE_NCFILE_TS stays undef: the runoff file carries flow only).
# The compile runs WITHOUT conda automatically (conda's NetCDF would break the
# link, "CURL_OPENSSL" errors); nf-config must resolve to opt_seq.
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
BUILD=${2:-${DEFAULT_BUILD}}
build_flags "${BUILD}"

for f in cppdefs.h param.h jobcomp; do
    [ -f "${CONFIG_DIR}/$f" ] || die "${CONFIG_DIR}/$f missing -- run step 05"
done
cd "${HCAST}"
cp "${CONFIG_DIR}"/{cppdefs.h,param.h,jobcomp} .
[ -f croco.in ] || cp "${CONFIG_DIR}/croco.in" . 2>/dev/null || true

L0=$(grep -n "^# define ${CONFIG_CPP}\$" cppdefs.h | head -1 | cut -d: -f1)
[ -n "$L0" ] || die "'# define ${CONFIG_CPP}' not in cppdefs.h -- re-run step 05"
LT=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +TIDES$")
LP=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +PSOURCE$")
LN=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +PSOURCE_NCFILE$")
LTS=$(cpp_line cppdefs.h "$L0" "^#  (define|undef) +PSOURCE_NCFILE_TS$")
onoff () { [ "$2" = 1 ] && echo "# define $1" || echo "# undef  $1"; }
sed -i "${LT}s/.*/$(onoff TIDES ${B_TIDES})/;${LP}s/.*/$(onoff PSOURCE ${B_RIVERS})/;${LN}s/.*/$(onoff PSOURCE_NCFILE ${B_RIVERS})/" cppdefs.h
sed -i "${LTS}s/.*/#  undef PSOURCE_NCFILE_TS/" cppdefs.h
SW="$(sed -n "${LT}p" cppdefs.h | tr -s " ") | $(sed -n "${LP}p" cppdefs.h | tr -s " ") | $(sed -n "${LN}p" cppdefs.h | tr -s " ")"

# ---- OpenMP (NTHREADS in domain.cfg, or 3rd argument): shared-memory tiles
#      SPLITTING_X x SPLITTING_ETA = NTHREADS (param.h: NPP, NSUB_X, NSUB_E);
#      jobcomp adds -fopenmp by itself when OPENMP is defined.
NT=${3:-${NTHREADS}}
LOMP=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +OPENMP$")
[ -n "$LOMP" ] || die "anchor OPENMP not found in cppdefs.h"
if [ "${NT}" -gt 1 ]; then
    read SX SE < <(awk -v n="${NT}" 'BEGIN{b=1; for(i=1;i*i<=n;i++) if(n%i==0) b=i; print b, n/b}')
    sed -i "${LOMP}s/.*/# define OPENMP\n#  define SPLITTING_X ${SX}\n#  define SPLITTING_ETA ${SE}/" cppdefs.h
    OMPTXT="OpenMP ${NT} threads (${SX} x ${SE} tiles)"
else
    sed -i "${LOMP}s/.*/# undef  OPENMP/" cppdefs.h
    OMPTXT="sequential"
fi
say "build ${BUILD}: ${SW} | ${OMPTXT}"

rm -rf croco Compile        # full rebuild: flags differ between builds (OpenMP, keys)
without_conda bash -c '
    source "${CROCO_ROOT}/env.sh" >/dev/null
    NFC=$(which nf-config || true)
    [ "${NFC}" = "${CROCO_PREFIX}/bin/nf-config" ] || { echo "!! nf-config is \"${NFC}\", expected ${CROCO_PREFIX}/bin/nf-config"; exit 1; }
    echo ">>> nf-config ${NFC} (prefix $(nf-config --prefix)); compiling ..."
    ./jobcomp > "compile_'"${BUILD}"'.log" 2>&1 || true
    tail -4 "compile_'"${BUILD}"'.log"'
[ -x croco ] || die "no croco binary -- read ${HCAST}/compile_${BUILD}.log"
cp croco "croco_${BUILD}"
ls -lh "${HCAST}/croco_${BUILD}"
