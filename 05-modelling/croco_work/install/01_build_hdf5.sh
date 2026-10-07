#!/bin/bash
# CROCO install — 1: HDF5 (sequential)
set -e

# Définition du répertoire de travail par défaut
: "${CROCO_ROOT:=${HOME}/croco_work}"

# Détermination automatique d'un nombre de jobs de compilation sûr (min(cores, RAM_GB/2))
if [ -z "${NJOBS:-}" ]; then
    CORES=$(nproc)
    RAM_GB=$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo)
    NJOBS=$(( RAM_GB / 2 ))
    [ "$NJOBS" -gt "$CORES" ] && NJOBS=$CORES
    [ "$NJOBS" -lt 1 ] && NJOBS=1
    echo ">>> auto NJOBS=$NJOBS  (cores=$CORES, RAM=${RAM_GB}GB, ~2GB/job)"
fi

D="${CROCO_ROOT}/install"
OPT="${CROCO_ROOT}/opt_seq"

cd "$D"

# 1. NETTOYAGE CRITIQUE : On supprime toute tentative précédente corrompue
echo ">>> Nettoyage des anciennes compilations..."
rm -rf hdf5-1.14.6

# 2. EXTRACTION
echo ">>> Extraction de HDF5..."
tar -xf hdf5-1.14.6.tar.gz
cd hdf5-1.14.6
mkdir -p build
cd build

# 3. CONFIGURATION SÉCURISÉE (Contournement strict de Conda pour GCC, GFortran et Zlib)
echo ">>> Configuration de HDF5 avec les outils système d'Ubuntu..."
CC=/usr/bin/gcc \
FC=/usr/bin/gfortran \
../configure \
    --prefix="$OPT" \
    --enable-fortran \
    --with-zlib=/usr \
    CPPFLAGS="-I/usr/include" \
    LDFLAGS="-L/usr/lib" \
    2>&1 | tee configure.log

# 4. COMPILATION
echo ">>> Compilation en cours (Jobs: $NJOBS)..."
make -j "$NJOBS" all 2>&1 | tee make.log

# 5. INSTALLATION
echo ">>> Installation dans $OPT..."
make install 2>&1 | tee install.log

echo ">>> HDF5 installé avec succès !"
echo ">>> next: bash install/02_build_netcdf_c.sh"


# zlib1g-dev libsz2-dev
