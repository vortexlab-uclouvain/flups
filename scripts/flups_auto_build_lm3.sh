#!/bin/bash
source $2

export FLUPS_MPICXX = mpic++
export FLUPS_MPICC = mpicc
export FLUPS_CXXFLAGS = -fopenmp -O3 -g -DNDEBUG -std=c++11 -mtune=skylake -DPROF
export FLUPS_CCFLAGS = -fopenmp -O3 -g -DNDEBUG -std=c99
export FLUPS_LDFLAGS = -fopenmp -lstdc++

MPICC=${FLUPS_MPICC} MPICXX=${FLUPS_MPICXX} \
CXXFLAGS=$(FLUPS_CXXFLAGS) CCFLAGS=$(FLUPS_CCFLAGS) LDFLAGS=$(FLUPS_LDFLAGS) \
HDF5_DIR=${HDF5_DIR} FFTW_DIR=${EBROOTFFTW} \
make -C $1 install -j 4

FILE1=libflups_a2a.a
if test -f "$1/lib/$FILE1"; then
    echo "$FILE1 exists."
    exit 0
else
    echo "$FILE1 does not exist."
    exit 1
fi