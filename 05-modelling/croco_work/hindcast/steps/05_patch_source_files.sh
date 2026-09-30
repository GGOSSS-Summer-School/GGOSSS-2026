#!/bin/bash
# =============================================================================
# Step 05 -- compile-time files: cppdefs.h, param.h, jobcomp
#   hindcast/steps/05_patch_source_files.sh NAME
#
# Copies CROCO's templates (code/croco/OCEAN) into the recipe folder
# hindcast/configs/NAME and makes the configuration edits, with every
# value taken from domain.cfg / obc.cfg / croco_grd.nc:
#   cppdefs.h  '# define BENGUELA_LR' -> '# define <CONFIG_CPP>'
#              OBC_EAST/WEST/NORTH/SOUTH of the block RIGHT BELOW that line = obc.cfg
#              ONLINE define, AROME undef, ERA_ECMWF define (ERA5 online forcing)
#              TIDES / PSOURCE / PSOURCE_NCFILE left undef here: step 07 sets them
#              per build (plain, plain_tides, plain_rivers, plain_tides_rivers)
#   param.h    '# elif defined <CONFIG_CPP>' + parameter (LLm0=xi_rho-2, MMm0=eta_rho-2,
#              N=N_LEVELS) inserted above '# else'
#   jobcomp    SOURCE1=${CROCO_MODEL_DIR}/OCEAN
# Only lines inside the REGIONAL block, after the config-name line, are touched --
# never the OBC blocks of the other example configurations further down.
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
need_obc
[ -f "${CF}/croco_grd.nc" ] || die "no grid -- run step 01"

cd "${CONFIG_DIR}"
cp "${CROCO_MODEL_DIR}/OCEAN/cppdefs.h" "${CROCO_MODEL_DIR}/OCEAN/param.h" \
   "${CROCO_MODEL_DIR}/OCEAN/croco.in"  "${CROCO_MODEL_DIR}/OCEAN/jobcomp" .

# ---------------- cppdefs.h ----------------
L0=$(grep -n "^# define BENGUELA_LR" cppdefs.h | head -1 | cut -d: -f1)
[ -n "$L0" ] || die "'# define BENGUELA_LR' not found in the cppdefs.h template"
LE=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +OBC_EAST");  LW=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +OBC_WEST")
LN=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +OBC_NORTH"); LS=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +OBC_SOUTH")
LON=$(cpp_line cppdefs.h "$L0" "^#  (define|undef) +ONLINE$")
LAR=$(cpp_line cppdefs.h "$L0" "^#   (define|undef) +AROME")
LER=$(cpp_line cppdefs.h "$L0" "^#   (define|undef) +ERA_ECMWF")
LT=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +TIDES$")
LP=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +PSOURCE$")
LPN=$(cpp_line cppdefs.h "$L0" "^# (define|undef) +PSOURCE_NCFILE$")
for v in LE LW LN LS LON LAR LER LT LP LPN; do [ -n "${!v}" ] || die "anchor $v not found in cppdefs.h"; done
[ $((LE - L0)) -lt 40 ] || die "OBC block is not right below the config name (line $LE)"

onoff () { [ "$2" = 1 ] && echo "# define $1" || echo "# undef  $1"; }
sed -i "${L0}s/.*/# define ${CONFIG_CPP}/" cppdefs.h
sed -i "${LE}s/.*/$(onoff OBC_EAST ${OBC_EAST})/;${LW}s/.*/$(onoff OBC_WEST ${OBC_WEST})/" cppdefs.h
sed -i "${LN}s/.*/$(onoff OBC_NORTH ${OBC_NORTH})/;${LS}s/.*/$(onoff OBC_SOUTH ${OBC_SOUTH})/" cppdefs.h
sed -i "${LON}s/.*/#  define ONLINE/;${LAR}s/.*/#   undef  AROME/;${LER}s/.*/#   define ERA_ECMWF/" cppdefs.h
sed -i "${LT}s/.*/# undef  TIDES/;${LP}s/.*/# undef  PSOURCE/;${LPN}s/.*/# undef  PSOURCE_NCFILE/" cppdefs.h

say "CHECK cppdefs.h -- config name, the OBC block below it, forcing switches:"
sed -n "${L0}p;${LT},${LS}p;${LON}p;${LAR}p;${LER}p" cppdefs.h

# ---------------- param.h ----------------
read XI ETA < <(python3 -c "
import xarray as xr; g = xr.open_dataset('${CF}/croco_grd.nc'); print(g.sizes['xi_rho'], g.sizes['eta_rho'])")
LLM0=$((XI - 2)); MMM0=$((ETA - 2))
LX=$(grep -n "YOUR REGIONAL CONFIG" param.h | head -1 | cut -d: -f1)
[ -n "$LX" ] || die "'YOUR REGIONAL CONFIG' not found in param.h"
sed -n "$((LX - 1))p" param.h | grep -q "^# else" || die "line above 'YOUR REGIONAL CONFIG' is not '# else'"
sed -i "$((LX - 1))i\\
# elif defined  ${CONFIG_CPP}\\
      parameter (LLm0=${LLM0},   MMm0=${MMM0},   N=${N_LEVELS})   ! ${CONFIG_NAME} ${XI}x${ETA}" param.h
say "CHECK param.h (expect LLm0=${LLM0}, MMm0=${MMM0}, N=${N_LEVELS}):"
cpp -P -DREGIONAL -D${CONFIG_CPP} param.h 2>/dev/null | grep "parameter (LLm0" | head -1 | sed 's/!.*//'

# ---------------- jobcomp ----------------
sed -i "s|^SOURCE1=.*|SOURCE1=${CROCO_MODEL_DIR}/OCEAN|" jobcomp
chmod +x jobcomp
say "CHECK jobcomp: $(grep '^SOURCE1=' jobcomp)"
ls "${CROCO_MODEL_DIR}/OCEAN/cppdefs.h" >/dev/null || die "SOURCE1 path does not exist"
say "next: hindcast/steps/06_patch_croco_in.sh ${CONFIG_NAME}"
