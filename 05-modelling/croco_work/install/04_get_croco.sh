#!/bin/bash
# CROCO install — 4: get the CROCO model + croco_pytools
# Versions pinned HERE only (download detail); install into clean folders:
#     code/croco  and  code/croco_pytools
set -e
: "${CROCO_ROOT:=${HOME}/croco_work}"
C="${CROCO_ROOT}/code"; mkdir -p "$C"; cd "$C"

CROCO_VERSION="2.1.3"
CROCO_PYTOOLS_VERSION="2.0.4"

echo ">>> get CROCO v${CROCO_VERSION} and croco_pytools v${CROCO_PYTOOLS_VERSION}"
echo ">>> install into (clean names, no version): ${C}/croco  and  ${C}/croco_pytools"

# --- CROCO model ---
if [ ! -d "${C}/croco" ]; then
    wget -nc "https://gitlab.inria.fr/croco-ocean/croco/-/archive/v${CROCO_VERSION}/croco-v${CROCO_VERSION}.tar.gz"
    tar -xzf "croco-v${CROCO_VERSION}.tar.gz"
    mv "croco-v${CROCO_VERSION}" croco
    rm -f "croco-v${CROCO_VERSION}.tar.gz"
    echo "    -> ${C}/croco"
else
    echo "    croco/ already present — skipping"
fi

# --- croco_pytools pre-processing toolbox ---
if [ ! -d "${C}/croco_pytools" ]; then
    wget -nc "https://gitlab.inria.fr/croco-ocean/croco_pytools/-/archive/v${CROCO_PYTOOLS_VERSION}/croco_pytools-v${CROCO_PYTOOLS_VERSION}.tar.gz"
    tar -xzf "croco_pytools-v${CROCO_PYTOOLS_VERSION}.tar.gz"
    mv "croco_pytools-v${CROCO_PYTOOLS_VERSION}" croco_pytools
    rm -f "croco_pytools-v${CROCO_PYTOOLS_VERSION}.tar.gz"
    echo "    -> ${C}/croco_pytools"
else
    echo "    croco_pytools/ already present — skipping"
fi

# --- compile the croco_pytools Fortran helpers (needs the ggosss26 conda env active) ---
echo ">>> compiling croco_pytools Fortran tools ..."
TOOLS="${C}/croco_pytools/prepro/Modules/tools_fort_routines"
if [ -d "$TOOLS" ]; then
    ( cd "$TOOLS" && make clean && make )
    ls "${C}"/croco_pytools/prepro/Modules/toolsf*.so >/dev/null 2>&1 \
        && echo ">>> croco_pytools tools compiled" \
        || echo "!! croco_pytools tools did NOT compile — activate the 'ggosss26' conda env and run 'make' in ${TOOLS}"
else
    echo "!! ${TOOLS} not found — check the croco_pytools download"
fi

echo ">>> CROCO + croco_pytools ready in ${C}/"
