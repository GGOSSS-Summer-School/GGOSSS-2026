# CROCO environment -- ~/croco_work
# Sourced by every script (you can also source it yourself in a terminal):
#   source ~/croco_work/env.sh
#
# This file sets only the SHARED paths + compilers. The per-track paths
# (configs / runs) come from hindcast/track.sh.

export CROCO_ROOT=${CROCO_ROOT:-${HOME}/croco_work}

# --- shared: model, tools, code, data, CLI ---
export CROCO_MODEL_DIR=${CROCO_ROOT}/code/croco
export CROCO_PYTOOLS_DIR=${CROCO_ROOT}/code/croco_pytools
export CROCO_DATA_ROOT=${CROCO_ROOT}/data
export GTOOLS_DIR=${CROCO_ROOT}/gtools          # the Python tools; the CLI is gtools/ggosss26.py

# --- compilers (for building CROCO and the croco_pytools Fortran helpers) ---
export CC=gcc
export FC=gfortran
export F90=gfortran
export F77=gfortran
# gfortran 10+ rejects argument mismatches that older versions only warned
# about; croco_pytools' legacy Fortran needs the older behaviour to build
export FFLAGS="-ffree-line-length-none -fallow-argument-mismatch -O3"

# --- NetCDF stack for compiling CROCO (sequential build in the repo) ---
export CROCO_PREFIX=${CROCO_ROOT}/opt_seq
export NETCDF=${CROCO_PREFIX}
export PATH=${CROCO_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CROCO_PREFIX}/lib

echo "CROCO environment set (root: ${CROCO_ROOT})"
