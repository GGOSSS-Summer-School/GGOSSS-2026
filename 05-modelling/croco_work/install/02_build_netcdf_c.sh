#!/bin/bash
# CROCO install — 2: NetCDF-C (needs: sudo apt install libcurl4-openssl-dev m4)
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
cd "$D"; tar -xvf netcdf-c-4.10.0.tar.gz; cd netcdf-c-4.10.0; mkdir -p build; cd build
CC=gcc FC=gfortran CPPFLAGS="-I$OPT/include" LDFLAGS="-L$OPT/lib" LIBS="-ldl" \
  ../configure --prefix="$OPT" --enable-hdf5 --disable-libxml2 --enable-curl 2>&1 | tee configure.log
make -j "$NJOBS" all 2>&1 | tee make.log; make install 2>&1 | tee install.log
echo ">>> next: bash install/03_build_netcdf_fortran.sh"
