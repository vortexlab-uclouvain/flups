#!/bin/bash
cd $1 

source $1/scripts/flups_loadmodule_lm3.sh

## Download and compile h3lpr
export H3LPR_MPICXX=mpic++
export H3LPR_MPICC=mpicc
export H3LPR_CXXFLAGS="-O3 -g -ggdb -fopenmp -DCOLOR_PROF"
export H3LPR_LDFLAGS="-fopenmp -lstdc++ -lm"
export H3LPR_PREFIX="${2}"

cd h3lpr
CXX=${H3LPR_MPICXX} CC=${H3LPR_MPICC} \
    CXXFLAGS=${H3LPR_CXXFLAGS} LDFLAGS=${H3LPR_LDFLAGS} \
    PREFIX=${H3LPR_PREFIX} \
    make install -j 
cd .. 
cd ${H3LPR_PREFIX} && rm -r lib/libh3lpr.so && cd - 
rm -rf h3lpr

## Compile flups
export FLUPS_MPICXX=mpic++
export FLUPS_MPICC=mpicc
export FLUPS_CXXFLAGS="-fopenmp -O3 -g -DNDEBUG -std=c++11 -mtune=skylake -DPROF -DHAVE_WISDOM=\\"$1/flups-wisdom/vortexbot.wsdm\\""
export FLUPS_CCFLAGS="-fopenmp -O3 -g -DNDEBUG -std=c99"
export FLUPS_LDFLAGS="-fopenmp -lstdc++"

ARCH_FILE=make_arch/make.default \
    CC=${FLUPS_MPICC} CXX=${FLUPS_MPICXX} \
    CXXFLAGS=${FLUPS_CXXFLAGS} CCFLAGS=${FLUPS_CCFLAGS} LDFLAGS=${FLUPS_LDFLAGS} \
    HDF5_DIR=${HDF5_DIR} FFTW_DIR=${EBROOTFFTW} H3LPR_DIR=${H3LPR_PREFIX} \
    make -C $1 install -j 4

FILE1=libflups_a2a.a
if test -f "$1/lib/$FILE1"; then
    echo "$FILE1 exists."
    exit 0
else
    echo "$FILE1 does not exist."
    exit 1
fi