#!/bin/bash
# CROCO install — 1: HDF5 (sequential)
set -e
: "${CROCO_ROOT:=${HOME}/croco_work}"
# Pick a safe NJOBS unless the user set one: min(cores, RAM_GB/2), at least 1
# (each parallel compile job needs ~2 GB RAM).
if [ -z "${NJOBS:-}" ]; then
    CORES=$(nproc)
    RAM_GB=$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo)
    NJOBS=$(( RAM_GB / 2 ))
    [ "$NJOBS" -gt "$CORES" ] && NJOBS=$CORES
    [ "$NJOBS" -lt 1 ] && NJOBS=1
    echo ">>> auto NJOBS=$NJOBS  (cores=$CORES, RAM=${RAM_GB}GB, ~2GB/job)"
fi
D="${CROCO_ROOT}/install"; OPT="${CROCO_ROOT}/opt_seq"
cd "$D"; tar -xvf hdf5-1.14.6.tar.gz; cd hdf5-1.14.6; mkdir -p build; cd build
CC=gcc FC=gfortran ../configure --prefix="$OPT" --enable-fortran --with-zlib=/usr 2>&1 | tee configure.log
make -j "$NJOBS" all 2>&1 | tee make.log; make install 2>&1 | tee install.log
echo ">>> next: bash install/02_build_netcdf_c.sh"
