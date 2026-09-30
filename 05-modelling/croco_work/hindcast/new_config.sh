#!/bin/bash
# =============================================================================
# new_config.sh -- create a new hindcast configuration (its domain.cfg).
#
#   hindcast/new_config.sh NAME LON_MIN LON_MAX LAT_MIN LAT_MAX [options]
#
#   options (all optional):
#     --res N            resolution 1/N degree            (default 12  -> 1/12 deg, ~9 km)
#     --start YYYY-MM-DD first hindcast day               (default 2018-06-01)
#     --end   YYYY-MM-DD end of the period, 00:00 of that day (default 2018-08-30;
#                        2026-09-01 covers all of June-August for --start 2026-06-01)
#     --tides            switch tides on  (TPXO tidal forcing)
#     --rivers           switch rivers on (GloFAS river discharge)
#     --ocean glorys|mercator  ocean source (default glorys = daily reanalysis; mercator =
#                        Mercator global analysis, for periods GLORYS does not cover yet)
#     --threads N        OpenMP threads (default 1 = sequential)
#     --hcast-days N     hindcast days per cycle (default 5; use a divisor of the period)
#     --parent NAME      OFFLINE NESTING: make this a child of configuration NAME. The
#                        ocean (ini/bry) then comes from the parent's own hindcast output
#                        instead of GLORYS, and the parent's ERA5 is reused. The box must
#                        lie inside the parent's; start/end/hcast-days/tides/rivers/threads
#                        default to the parent's; 1-day spin-up; a boundary sponge is set.
#     --force            overwrite an existing domain.cfg
#
# examples
#   hindcast/new_config.sh GoG_12  -10 12.5 -6 8
#   hindcast/new_config.sh IGOG_25  4  12.5 -6 5.5 --res 25 --tides --rivers
#   hindcast/new_config.sh Canary_12 -22 -15.5 14 24 --start 2025-12-02 --end 2026-01-09
#   hindcast/new_config.sh IGOG_36   4  12.5 -6 5.5 --res 36 --parent GoG_12   (nested child)
#
# This writes hindcast/configs/NAME/domain.cfg -- the ONLY file you edit by
# hand. Everything else (download boxes, time step, data months, cycles,
# open boundaries, parameter files, cppdefs.h/param.h/croco.in edits) is
# derived from it by hindcast/lib/common.sh and the steps. Next:
#   hindcast/build_config.sh NAME
# =============================================================================
set -e
usage () { sed -n '2,36p' "$0"; exit "${1:-1}"; }
[ $# -ge 5 ] || usage
NAME=$1; LON_MIN=$2; LON_MAX=$3; LAT_MIN=$4; LAT_MAX=$5; shift 5
RES_INV=12; START=2018-06-01; END=2018-08-30; TIDES=0; RIVERS=0; FORCE=0; OCEAN=glorys; THREADS=1; HCD=5
PARENT=""; SPINUP=2; SPONGE_LINE='#SPONGE="50000. 400."   # X_SPONGE [m], V_SPONGE [m2/s]; default "0. 0." (off)'
# a nested child inherits the parent's period and physics unless given explicitly
for ((i=1; i<=$#; i++)); do [ "${!i}" = --parent ] && { j=$((i+1)); PARENT=${!j}; }; done
if [ -n "${PARENT}" ]; then
    PCFG="$(cd "$(dirname "$0")" && pwd)/configs/${PARENT}/domain.cfg"
    [ -f "${PCFG}" ] || { echo "!! --parent ${PARENT}: ${PCFG} not found"; exit 1; }
    eval "$(bash -c "TIDES=0; RIVERS=0; NTHREADS=1; HCAST_DAYS=5; source '${PCFG}'; echo START=\$HC_START END=\$HC_END TIDES=\$TIDES RIVERS=\$RIVERS THREADS=\$NTHREADS HCD=\$HCAST_DAYS; echo PLON=\$LON_MIN,\$LON_MAX,\$LAT_MIN,\$LAT_MAX")"
    # the parent's daily means are stamped at noon: its last one (parent end - 12 h) must
    # bracket the child's last time step, so the child ends one day before the parent
    END=$(date -u -d "${END} - 1 day" +%Y-%m-%d)
    _n=$(( ( $(date -u -d "${END}" +%s) - $(date -u -d "${START}" +%s) ) / 86400 ))
    [ "${_n}" -lt "${HCD}" ] && HCD=${_n}
fi
while [ $# -gt 0 ]; do
    case "$1" in
        --res)    RES_INV=$2; shift 2 ;;
        --start)  START=$2; shift 2 ;;
        --end)    END=$2; shift 2 ;;
        --tides)  TIDES=1; shift ;;
        --rivers) RIVERS=1; shift ;;
        --ocean)  OCEAN=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --hcast-days) HCD=$2; shift 2 ;;
        --parent) shift 2 ;;          # read above
        --force)  FORCE=1; shift ;;
        -h|--help) usage 0 ;;
        *) echo "!! unknown option $1"; usage ;;
    esac
done

[[ "${OCEAN}" == glorys || "${OCEAN}" == mercator ]] || { echo "!! --ocean must be glorys or mercator"; exit 1; }
[[ "${NAME}" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || { echo "!! NAME must start with a letter and contain only letters, digits, _ (it becomes a cpp key)"; exit 1; }
for v in "${LON_MIN}" "${LON_MAX}" "${LAT_MIN}" "${LAT_MAX}" "${RES_INV}"; do
    [[ "$v" =~ ^-?[0-9]+(\.[0-9]+)?$ ]] || { echo "!! '$v' is not a number"; exit 1; }
done
awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v c="${LAT_MIN}" -v d="${LAT_MAX}" \
    'BEGIN{exit !(a<b && c<d && a>=-180 && b<=180 && c>=-80 && d<=90)}' \
    || { echo "!! need LON_MIN < LON_MAX (-180..180) and LAT_MIN < LAT_MAX"; exit 1; }
date -u -d "${START}" >/dev/null && date -u -d "${END}" >/dev/null
if [ -n "${PARENT}" ]; then
    IFS=, read A B C D <<< "${PLON}"
    awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v c="${LAT_MIN}" -v d="${LAT_MAX}" -v A="$A" -v B="$B" -v C="$C" -v D="$D" \
        'BEGIN{exit !(a>=A && b<=B && c>=C && d<=D)}' \
        || { echo "!! the child box must lie inside the parent ${PARENT} box (lon $A..$B, lat $C..$D)"; exit 1; }
    SPINUP=1
    # sponge ~7 child cells wide (e.g. 1/36 deg -> ~20 km): damps the sharper parent
    # structures arriving at the open edges so they do not reflect inward
    SPW=$(awk -v r="${RES_INV}" 'BEGIN{w=7*111000/r; printf "%d", int(w/5000+0.5)*5000}')
    SPONGE_LINE="SPONGE=\"${SPW}. 400.\"      # X_SPONGE [m], V_SPONGE [m2/s]: on for a nested child"
fi

HC_HOME="$(cd "$(dirname "$0")" && pwd)"
CROCO_ROOT="${CROCO_ROOT:-$(dirname "${HC_HOME}")}"
DIR="${CROCO_ROOT}/hindcast/configs/${NAME}"
if [ -f "${DIR}/domain.cfg" ] && [ "${FORCE}" != 1 ]; then
    echo "!! ${DIR}/domain.cfg exists -- edit it, or use --force to overwrite"; exit 1
fi
mkdir -p "${DIR}"

# rough size estimate (the real size comes from the grid, step 01)
NX=$(awk -v a="${LON_MIN}" -v b="${LON_MAX}" -v r="${RES_INV}" 'BEGIN{printf "%d", (b-a)*r+1.5+2}')
NY=$(awk -v a="${LAT_MIN}" -v b="${LAT_MAX}" -v r="${RES_INV}" 'BEGIN{printf "%d", (b-a)*r+1.5+2}')

cat > "${DIR}/domain.cfg" <<EOF
# =============================================================================
# hindcast/configs/${NAME}/domain.cfg
# The ONE hand-written file of this configuration (bash syntax: no spaces
# around '='). Created by new_config.sh on $(date -u +%Y-%m-%d); edit freely.
# Everything else is derived from it -- see hindcast/lib/common.sh and
# hindcast/README.md. After editing the box or resolution, rebuild from
# step 01 (hindcast/build_config.sh ${NAME}).
# =============================================================================

# ---- name ---------------------------------------------------------------------
# Folder name, run folders, and (upper-cased: ${NAME^^}) the cppdefs.h/param.h key.
CONFIG_NAME=${NAME}

# ---- horizontal domain (degrees; longitudes -180..180, east positive) ----------
# The grid box. Estimated size ~${NX} x ${NY} rho points (the exact size is read
# from croco_grd.nc after step 01). The GLORYS download box (+1.5 deg) and the
# ERA5 box are derived from these four numbers.
LON_MIN=${LON_MIN}
LON_MAX=${LON_MAX}
LAT_MIN=${LAT_MIN}
LAT_MAX=${LAT_MAX}

# ---- resolution -----------------------------------------------------------------
# Grid spacing = 1/RES_INV degree (12 -> 1/12 deg ~ 9.3 km at the equator).
# The time step DT is derived: 300 s at 1/12 deg (scaled by 12/RES_INV, and by
# cos(latitude) poleward of 35 deg). Step 02 prints the resulting CFL numbers.
RES_INV=${RES_INV}
#DT=300                 # uncomment to force the baroclinic time step [s]
#NDTFAST=60             # barotropic sub-steps per DT

# ---- vertical grid (must be the same everywhere; written for you into
#      crocotools_param.py, croco.in S-coord and param.h) --------------------------
THETA_S=7
THETA_B=2
N_LEVELS=50
SIGMA_HC=200

# ---- hindcast period ---------------------------------------------------------------
# Cycled as SPINUP_DAYS spin-up + HCAST_DAYS hindcast per cycle; the number of
# cycles and the GLORYS/ERA5 months to download are derived (one neighbour
# month is added automatically for the spin-up and the +/-1 day boundary pad).
HC_START=${START}
HC_END=${END}
SPINUP_DAYS=${SPINUP}
HCAST_DAYS=${HCD}
TEST_DAYS=7             # length of the single proof run (steps 06-08), from day 2

# ---- open boundaries -----------------------------------------------------------------
# Empty = decided automatically from the land mask by step 02 (an edge is open
# if more than half of it is ocean) and written to obc.cfg. To force it, give
# "south west east north" as 1 (open) / 0 (closed), e.g. OBC="1 1 0 0".
OBC=""
${SPONGE_LINE}

# ---- optional physics (compile-time: each combination is its own binary) -------------
# TIDES=1  : TPXO tidal forcing. Needs DATASETS_CROCOTOOLS/TPXO7.
# RIVERS=1 : river point sources with GloFAS daily discharge (the real discharge
#            of the period, not a climatology). Needs ~/.ewdsapirc.
# Both 0 = the plain hindcast (always built and tested first).
# The driver uses these by default; --tides/--no-tides/--rivers/--no-rivers override.
TIDES=${TIDES}
RIVERS=${RIVERS}

# ---- ocean source ----------------------------------------------------------------------
# glorys   : GLORYS12 daily reanalysis (default; ends a few months before the present)
# mercator : Mercator global analysis (anfc, 1/12 deg, same variables), for recent periods
OCEAN_SOURCE=${OCEAN}

# ---- offline nesting ------------------------------------------------------------------------
# PARENT=<configuration> makes this a CHILD: its ini/bry are built from the parent's own
# hindcast output (daily means converted to GLORYS format by step 03) instead of GLORYS, and
# the parent's ERA5 files are reused (OCEAN_SOURCE is then ignored). Run the parent first.
# PARENT_RUN=<folder> picks one parent cycle folder (default: those overlapping the period).
PARENT=${PARENT}
#PARENT_RUN=

# ---- parallel run --------------------------------------------------------------------------
# NTHREADS > 1 builds and runs CROCO with OpenMP on that many threads (<= cores: nproc)
NTHREADS=${THREADS}

# ---- GloFAS (rivers) -----------------------------------------------------------------------
# auto = consolidated if it covers the whole period, else intermediate (recent months)
GLOFAS_PRODUCT_TYPE=auto

# ---- river options (only used when RIVERS=1) ------------------------------------------
RIVER_TEMP=27.0         # temperature of river water [degC] (tropical; ~15-20 at mid-latitudes)
RIVER_SALT=0.0          # fresh water
RIVER_QMIN=100          # ignore rivers with mean GloFAS discharge below this [m3/s]
RIVER_RADIUS=0.25       # search radius around a river mouth for the GloFAS outlet cell [deg]
RIVER_MARGIN=0.0        # accept mouths this far outside the grid [deg]
RIVER_EXTRA=""          # optional file with extra mouths: lines "name lon lat"
EOF
echo ">>> wrote ${DIR}/domain.cfg"
echo "    box lon ${LON_MIN}..${LON_MAX} lat ${LAT_MIN}..${LAT_MAX}, 1/${RES_INV} deg (~${NX} x ${NY}), ${START} -> ${END}, tides=${TIDES} rivers=${RIVERS}, ocean=${OCEAN}, threads=${THREADS}"
[ -n "${PARENT}" ] && echo "    NESTED in ${PARENT} (1-day spin-up, sponge on). Run the ${PARENT} hindcast first."
echo "    next: ${HC_HOME}/build_config.sh ${NAME}"
