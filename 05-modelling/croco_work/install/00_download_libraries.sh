#!/bin/bash
# CROCO install — 0: download library sources
set -e
: "${CROCO_ROOT:=${HOME}/croco_work}"
D="${CROCO_ROOT}/install"; mkdir -p "$D"; cd "$D"
wget -nc https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
wget -nc https://downloads.unidata.ucar.edu/netcdf-c/4.10.0/netcdf-c-4.10.0.tar.gz
wget -nc https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.2/netcdf-fortran-4.6.2.tar.gz
echo ">>> next: bash install/01_build_hdf5.sh"
