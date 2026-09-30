# =============================================================================
# hindcast/lib/common.sh -- loaded by EVERY hindcast step, the builder and the
# cycling driver:
#
#     source "${HC_LIB}/common.sh" <CONFIG_NAME>
#
# It (1) sets the session (env.sh + hindcast/track.sh + the Python env),
# (2) reads the ONE hand-written file of a configuration,
#         hindcast/configs/<CONFIG_NAME>/domain.cfg
# and (3) DERIVES everything else, so no number is typed twice:
#
#   derived                     from
#   --------------------------  ------------------------------------------------
#   CONFIG_CPP (GOG_12)         CONFIG_NAME, upper-cased (cppdefs.h/param.h key)
#   CONFIG_DIR, HCAST, CF       hindcast/configs/<NAME>, hindcast/scratch/<NAME>
#   RES (deg)                   1/RES_INV, full precision
#   EXTENTS (GLORYS box)        grid box + GLORYS_PAD (1.5 deg)
#   ERA5_BOX                    grid box (the ERA5 downloader adds its own 2 deg)
#   DT, NDTFAST                 resolution + latitude (300 s at 1/12 deg), unless DT set
#   GLORYS_MONTH_START/END      period - spin-up - 2 d  ...  period end + 2 d
#   ERA5_MONTH_START/END        period - spin-up        ...  period end
#   NCYCLES                     period length / HCAST_DAYS (rounded up)
#   TEST_START/TEST_END         proof run: day 2 of the period, TEST_DAYS long
#   OBC_SOUTH/WEST/EAST/NORTH   the land mask (obc.cfg, written by step 02), or OBC override
#   DEFAULT_BUILD               plain[_tides][_rivers] from TIDES / RIVERS
# =============================================================================

_cfg_name="${1:-${CONFIG_NAME:-}}"
if [ -z "${_cfg_name}" ]; then
    echo "!! give the configuration name, e.g.:  $0 GoG_12"
    echo "   (create one with: hindcast/new_config.sh <NAME> <lon_min> <lon_max> <lat_min> <lat_max>)"
    return 1 2>/dev/null || exit 1
fi

die () { echo "!! $*" >&2; exit 1; }
say () { echo ">>> $*"; }

# ---- 1. session: shared paths (env.sh) + hindcast track + Python env ---------
HC_LIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HC_HOME="$(dirname "${HC_LIB}")"                          # .../hindcast
export CROCO_ROOT="${CROCO_ROOT:-$(dirname "${HC_HOME}")}"
source "${CROCO_ROOT}/env.sh" >/dev/null
source "${HC_HOME}/track.sh" >/dev/null
export HC_LIB HC_HOME

# the Python tools need the conda env (xarray, netCDF4, copernicusmarine, cdsapi)
CONDA_ENV_NAME="${CONDA_ENV_NAME:-ggosss26}"
if [ "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV_NAME}" ] && [ -z "${NO_CONDA:-}" ]; then
    for _c in "${CONDA_EXE%/bin/conda}" "${HOME}/miniconda3" "${HOME}/anaconda3" "${HOME}/mambaforge"; do
        if [ -n "${_c}" ] && [ -f "${_c}/etc/profile.d/conda.sh" ]; then
            source "${_c}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_NAME}" && break
        fi
    done
fi

# ---- 2. defaults, then the configuration file --------------------------------
CONFIG_NAME="${_cfg_name}"
RES_INV=12                  # resolution = 1/RES_INV degree
THETA_S=7; THETA_B=2; N_LEVELS=50; SIGMA_HC=200     # vertical grid (sigma_params)
SPINUP_DAYS=2; HCAST_DAYS=5; TEST_DAYS=7
YORIG=1993                  # GLORYS/ERA5 reanalysis time origin
GLORYS_PAD=1.5              # deg added around the grid for the GLORYS download
GLORYS_PRODUCT=cmems_mod_glo_phy_my_0.083deg_P1D-m
OCEAN_SOURCE=glorys        # glorys (daily reanalysis) or mercator (global analysis, for recent periods)
NTHREADS=1                 # >1 = OpenMP build + run on that many threads
GLOFAS_PRODUCT_TYPE=auto   # auto / consolidated / intermediate; GLOFAS_VERSION=version_4_0
GLOFAS_VERSION=version_4_0
NDTFAST=60; DT=""           # DT empty = automatic
OBC=""                      # empty = from the land mask; or "south west east north" e.g. "1 1 0 0"
SPONGE="0. 0."              # X_SPONGE[m] V_SPONGE[m2/s]; "50000. 400." if energy piles up at an open edge
TIDES=0; RIVERS=0           # optional physics, see domain.cfg
TPXO_DIR=""                 # empty = ${CROCO_DATA_ROOT}/DATASETS_CROCOTOOLS/TPXO7
RIVER_TEMP=27.0; RIVER_SALT=0.0; RIVER_QMIN=100; RIVER_RADIUS=0.25; RIVER_MARGIN=0.0; RIVER_EXTRA=""
LON_MIN=""; LON_MAX=""; LAT_MIN=""; LAT_MAX=""; HC_START=""; HC_END=""
PARENT=""                   # offline nesting: name of the parent configuration (empty = not nested)
PARENT_RUN=""               # parent cycle folder to nest in (empty = the one starting on HC_START)

export CONFIG_DIR="${CROCO_CONFIGS_ROOT}/${CONFIG_NAME}"
CFG_FILE="${CONFIG_DIR}/domain.cfg"
[ -f "${CFG_FILE}" ] || die "${CFG_FILE} not found -- create it: ${HC_HOME}/new_config.sh ${CONFIG_NAME} <lon_min> <lon_max> <lat_min> <lat_max>"
source "${CFG_FILE}"
[ "${CONFIG_NAME}" = "${_cfg_name}" ] || die "CONFIG_NAME in ${CFG_FILE} (${CONFIG_NAME}) differs from the folder name (${_cfg_name})"

for v in LON_MIN LON_MAX LAT_MIN LAT_MAX HC_START HC_END; do
    [ -n "${!v}" ] || die "${v} is not set in ${CFG_FILE}"
done
awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v c="${LAT_MIN}" -v d="${LAT_MAX}" \
    'BEGIN{exit !(a<b && c<d && a>=-180 && b<=180 && c>=-80 && d<=90)}' \
    || die "bad box in ${CFG_FILE}: lon ${LON_MIN}..${LON_MAX}, lat ${LAT_MIN}..${LAT_MAX} (need min<max, -180..180)"
date -u -d "${HC_START}" >/dev/null 2>&1 && date -u -d "${HC_END}" >/dev/null 2>&1 \
    || die "HC_START/HC_END must be YYYY-MM-DD in ${CFG_FILE}"
[ "$(date -u -d "${HC_END}" +%s)" -gt "$(date -u -d "${HC_START}" +%s)" ] || die "HC_END must be after HC_START"

# ---- 3. derived ---------------------------------------------------------------
export CONFIG_NAME
export CONFIG_CPP=$(echo "${CONFIG_NAME}" | tr 'a-z' 'A-Z' | sed 's/[^A-Z0-9_]/_/g; s/^\([0-9]\)/C\1/')
export HCAST="${CROCO_RUNS_ROOT}/${CONFIG_NAME}"
export CF="${HCAST}/CROCO_FILES"
export OUTPUT_ROOT="${CROCO_ROOT}/hindcast/model-runs/${CONFIG_NAME}"
export LON_MIN LON_MAX LAT_MIN LAT_MAX RES_INV THETA_S THETA_B N_LEVELS SIGMA_HC
export SPINUP_DAYS HCAST_DAYS TEST_DAYS YORIG GLORYS_PRODUCT NDTFAST SPONGE
export TIDES RIVERS RIVER_TEMP RIVER_SALT RIVER_QMIN RIVER_RADIUS RIVER_MARGIN RIVER_EXTRA
export HC_START HC_END OBC GLORYS_PAD OCEAN_SOURCE NTHREADS GLOFAS_PRODUCT_TYPE GLOFAS_VERSION
export RES=$(echo "1/${RES_INV}" | bc -l)
export EXTENTS=$(awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v c="${LAT_MIN}" -v d="${LAT_MAX}" -v p="${GLORYS_PAD}" \
    'BEGIN{printf "%.2f,%.2f,%.2f,%.2f", a-p, b+p, c-p, d+p}')
export ERA5_BOX="${LON_MIN},${LON_MAX},${LAT_MIN},${LAT_MAX}"
export TPXO_DIR="${TPXO_DIR:-${CROCO_DATA_ROOT}/DATASETS_CROCOTOOLS/TPXO7}"
export GLOFAS_DIR="${HCAST}/downloaded_data/GLOFAS"
export ERA5_DIR="${HCAST}/downloaded_data/ERA5"
# ---- offline nesting: PARENT set -> the ocean comes from the parent's own output --
#   ocean  the parent's daily averages, converted to GLORYS/Mercator format by step 03
#          (gtools/nesting.py), in downloaded_data/PARENT_<parent> as YYYY_MM.nc files:
#          the child's ini/bry are then built by the usual steps from those files
#   ERA5   the parent's ERA5 files are reused (CROCO interpolates them online onto
#          the finer grid; the parent's ERA5 box covers the child, which lies inside)
if [ -n "${PARENT}" ]; then
    [ -f "${CROCO_CONFIGS_ROOT}/${PARENT}/domain.cfg" ] || die "PARENT=${PARENT}: no ${CROCO_CONFIGS_ROOT}/${PARENT}/domain.cfg"
    # the child box must lie inside the parent box (checked on the domain.cfg boxes)
    read P_LON_MIN P_LON_MAX P_LAT_MIN P_LAT_MAX <<< "$(bash -c "source '${CROCO_CONFIGS_ROOT}/${PARENT}/domain.cfg'; echo \$LON_MIN \$LON_MAX \$LAT_MIN \$LAT_MAX")"
    awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v c="${LAT_MIN}" -v d="${LAT_MAX}" \
        -v A="${P_LON_MIN}" -v B="${P_LON_MAX}" -v C="${P_LAT_MIN}" -v D="${P_LAT_MAX}" \
        'BEGIN{exit !(a>=A && b<=B && c>=C && d<=D)}' \
        || die "child box ${LON_MIN}..${LON_MAX} / ${LAT_MIN}..${LAT_MAX} is not inside the parent ${PARENT} box ${P_LON_MIN}..${P_LON_MAX} / ${P_LAT_MIN}..${P_LAT_MAX}"
    # ... and at least NEST_MARGIN deg inside on EVERY side (land edges too): the ini/bry reader
    # (ibc_class.py) asks the source data to cover the child grid + a 0.2 deg buffer, else
    # "ERROR: The data does not cover the entire grid" at step 04
    NEST_MARGIN=${NEST_MARGIN:-0.3}
    _bad=$(awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v c="${LAT_MIN}" -v d="${LAT_MAX}" \
        -v A="${P_LON_MIN}" -v B="${P_LON_MAX}" -v C="${P_LAT_MIN}" -v D="${P_LAT_MAX}" -v m="${NEST_MARGIN}" 'BEGIN{
        if (a-A<m) printf " LON_MIN<=%g", A+m; if (B-b<m) printf " LON_MAX<=%g", B-m
        if (c-C<m) printf " LAT_MIN>=%g", C+m; if (D-d<m) printf " LAT_MAX<=%g", D-m }' | sed 's/LON_MIN<=/LON_MIN>=/')
    [ -z "${_bad}" ] || die "child box too close to the parent ${PARENT} edge (need ${NEST_MARGIN} deg on every side, land edges too:
   the ini/bry reader wants the child grid + 0.2 deg covered by the parent data). Set in ${CFG_FILE}:${_bad}
   then rebuild the child grid: build_config.sh ${CONFIG_NAME} --only 1, --only 2, then --from 4"
    export PARENT PARENT_RUN P_LON_MIN P_LON_MAX P_LAT_MIN P_LAT_MAX
    export PARENT_OUTPUT_ROOT="${CROCO_ROOT}/hindcast/model-runs/${PARENT}"
    OCEAN_SOURCE=parent
    export ERA5_DIR="${CROCO_RUNS_ROOT}/${PARENT}/downloaded_data/ERA5"
fi
case "${OCEAN_SOURCE}" in
    glorys)   export OCEAN_TAG=GLORYS ;;
    mercator) export OCEAN_TAG=MERCATOR ;;
    parent)   export OCEAN_TAG=NEST ;;
    *) die "OCEAN_SOURCE must be glorys or mercator (domain.cfg)" ;;
esac
export OCEAN_SOURCE
# ocean data folder (variable kept as GLORYS_DIR for all steps):
# downloaded_data/GLORYS, /MERCATOR, or /PARENT_<parent> for a nested child
if [ "${OCEAN_SOURCE}" = parent ]; then
    export GLORYS_DIR="${HCAST}/downloaded_data/PARENT_${PARENT}"
else
    export GLORYS_DIR="${HCAST}/downloaded_data/${OCEAN_TAG}"
fi
# OpenMP run environment (the binary is built with OPENMP by step 07 when NTHREADS > 1)
if [ "${NTHREADS}" -gt 1 ] 2>/dev/null; then
    export OMP_NUM_THREADS=${NTHREADS} OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
    ulimit -s unlimited 2>/dev/null || true
fi

# time step: 300 s at 1/12 deg up to 35 deg latitude (a well-tested 1/12 deg value),
# scaled with the cell size (1/RES_INV, and cos(lat) beyond 35 deg), then rounded
# DOWN to a divisor of 3600 so that hourly/daily output and NTIMES are exact.
if [ -z "${DT}" ]; then
    DT=$(awk -v r="${RES_INV}" -v a="${LAT_MIN}" -v b="${LAT_MAX}" 'BEGIN{
        pi=atan2(0,-1); lat=(a<0?-a:a); if ((b<0?-b:b)>lat) lat=(b<0?-b:b)
        f=cos(lat*pi/180)/cos(35*pi/180); if (f>1) f=1
        raw=300*(12/r)*f
        n=split("600 480 450 400 360 300 240 225 200 180 150 144 120 100 90 80 75 72 60 50 48 45 40 36 30 25 24 20 18 16 15 12 10",L," ")
        for (i=1;i<=n;i++) if (L[i]+0<=raw) {print L[i]; exit}
        print 10}')
fi
export DT
export NTIMES_PER_DAY=$(( 86400 / DT ))

_d () { date -u -d "$1" +%Y-%m-%d; }
export SPIN_FIRST=$(_d "${HC_START} - ${SPINUP_DAYS} days")
export GLORYS_MONTH_START=$(date -u -d "${SPIN_FIRST} - 2 days" +%Y-%m)
export GLORYS_MONTH_END=$(date -u -d "${HC_END} + 2 days" +%Y-%m)
export ERA5_MONTH_START=$(date -u -d "${SPIN_FIRST}" +%Y-%m)
export ERA5_MONTH_END=$(date -u -d "${HC_END}" +%Y-%m)
_ndays=$(( ( $(date -u -d "${HC_END}" +%s) - $(date -u -d "${HC_START}" +%s) ) / 86400 ))
export NCYCLES=$(( (_ndays + HCAST_DAYS - 1) / HCAST_DAYS ))
export TEST_START=$(_d "${HC_START} + 1 day")                  # not the 1st: records on both sides
export TEST_END=$(_d "${TEST_START} + ${TEST_DAYS} days")

# open boundaries: explicit override, else the mask verdict written by step 02
if [ -n "${OBC}" ]; then
    read OBC_SOUTH OBC_WEST OBC_EAST OBC_NORTH <<< "${OBC}"
elif [ -f "${CONFIG_DIR}/obc.cfg" ]; then
    source "${CONFIG_DIR}/obc.cfg"
fi
export OBC_SOUTH OBC_WEST OBC_EAST OBC_NORTH

_b="plain"; [ "${TIDES}" = 1 ] && _b="${_b}_tides"; [ "${RIVERS}" = 1 ] && _b="${_b}_rivers"
export DEFAULT_BUILD="${_b}"          # what TIDES/RIVERS in domain.cfg ask for

mkdir -p "${CONFIG_DIR}" "${CF}" "${GLORYS_DIR}" "${ERA5_DIR}/for_croco"

# ---- helpers ------------------------------------------------------------------
# build name -> the two switches
build_flags () {   # usage: build_flags plain_tides  -> sets B_TIDES B_RIVERS
    case "$1" in
        plain)              B_TIDES=0; B_RIVERS=0 ;;
        plain_tides)        B_TIDES=1; B_RIVERS=0 ;;
        plain_rivers)       B_TIDES=0; B_RIVERS=1 ;;
        plain_tides_rivers) B_TIDES=1; B_RIVERS=1 ;;
        *) die "unknown build '$1' (plain | plain_tides | plain_rivers | plain_tides_rivers)" ;;
    esac
}
# first line >= line of the config name matching an extended regex (cppdefs.h)
cpp_line () { awk -v s="$2" -v p="$3" 'NR>=s && $0 ~ p {print NR; exit}' "$1"; }
need_obc () {
    [ -n "${OBC_SOUTH}" ] || die "open boundaries unknown -- run steps/02_check_grid.sh ${CONFIG_NAME} (or set OBC in domain.cfg)"
}
show_config () {
    cat <<EOF
============================================================
 configuration ${CONFIG_NAME}   (cpp key ${CONFIG_CPP})
   box         lon ${LON_MIN} .. ${LON_MAX}   lat ${LAT_MIN} .. ${LAT_MAX}
   resolution  1/${RES_INV} deg, ${N_LEVELS} levels (theta_s ${THETA_S}, theta_b ${THETA_B}, hc ${SIGMA_HC})
   time step   dt ${DT} s, NDTFAST ${NDTFAST}  (${NTIMES_PER_DAY} steps/day)
   period      ${HC_START} -> ${HC_END}  = ${NCYCLES} cycles of ${SPINUP_DAYS} d spin-up + ${HCAST_DAYS} d
   test run    ${TEST_START} -> ${TEST_END}
   ocean       $( [ -n "${PARENT}" ] && echo "NESTED in ${PARENT} (its daily means, converted by step 03)" || echo "${OCEAN_SOURCE} (${OCEAN_TAG}), box ${EXTENTS}" ), months ${GLORYS_MONTH_START} .. ${GLORYS_MONTH_END}
   ERA5        box ${ERA5_BOX} (+2 deg), months ${ERA5_MONTH_START} .. ${ERA5_MONTH_END}
   boundaries  south=${OBC_SOUTH:-?} west=${OBC_WEST:-?} east=${OBC_EAST:-?} north=${OBC_NORTH:-?}  $( [ -n "${OBC}" ] && echo "(override)" || echo "(from mask)")
   physics     tides=${TIDES} rivers=${RIVERS}  -> build ${DEFAULT_BUILD}, threads ${NTHREADS}
   folders     recipe ${CONFIG_DIR}
               work   ${HCAST}
               runs   ${OUTPUT_ROOT}
============================================================
EOF
}

# write the psource_ncfile block of a croco.in from CROCO_FILES/for_croco_in.txt
# (make_river_run.py output), replacing the WHOLE old block (a stale block would keep old river indices):
# short path CROCO_FILES/croco_runoff.nc, Nsrc, one row per river with RIVER_TEMP/SALT
write_psource_block () {   # usage: write_psource_block <croco.in>
    [ -f "${CF}/for_croco_in.txt" ] || die "no ${CF}/for_croco_in.txt -- run steps/10_rivers_glofas.sh ${CONFIG_NAME}"
    CROCO_IN="$1" python3 - <<'PYEOF'
import os
E = os.environ
lines = open(os.path.join(E["CF"], "for_croco_in.txt")).read().split("\n")
k = [i for i, l in enumerate(lines) if l.startswith("psource_ncfile:")][0]
n = int(lines[k + 2].split()[0])
rows = []
for l in lines[k + 3:k + 3 + n]:
    ii, jj, dsrc, qdir = l.split()[:4]
    rows.append("                      %5s %5s   %s     %3s    T T   %.1f  %.1f"
                % (ii, jj, dsrc, qdir, float(E["RIVER_TEMP"]), float(E["RIVER_SALT"])))
block = ["psource_ncfile:   Nsrc  Isrc  Jsrc  Dsrc qbardir  Lsrc  Tsrc   runoff file name",
         "    CROCO_FILES/croco_runoff.nc", "                 %d" % n] + rows + [""]
src = open(E["CROCO_IN"]).read().split("\n")
a = [i for i, l in enumerate(src) if l.startswith("psource_ncfile:")][0]
b = a + 1
while b < len(src) and src[b].strip() != "":
    b += 1
src[a:b + 1] = block
open(E["CROCO_IN"], "w").write("\n".join(src))
print("psource_ncfile block: %d river(s), T=%s S=%s" % (n, E["RIVER_TEMP"], E["RIVER_SALT"]))
PYEOF
}

# run a command with conda removed from the environment (compile with opt_seq NetCDF)
without_conda () {
    env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV -u CONDA_SHLVL -u CONDA_PROMPT_MODIFIER \
        PATH="$(echo "${PATH}" | tr ':' '\n' | grep -v -E "conda3?/(envs|bin|condabin)|/(anaconda3|mambaforge)/" | paste -sd:)" \
        NO_CONDA=1 "$@"
}
