#!/bin/bash
source $2

# Export the variables needed for flups
export FLUPS_MPICXX=mpic++
export FLUPS_MPICC=mpicc
export FLUPS_CXXFLAGS="-fopenmp -O3 -g -DNDEBUG -std=c++11 -mtune=skylake -DPROF"
export FLUPS_CCFLAGS="-fopenmp -O3 -g -DNDEBUG -std=c99"
export FLUPS_LDFLAGS="-fopenmp -lstdc++"


# Clone and compile google source test 
git clone https://github.com/google/googletest.git $1/googletest/
cd $1/googletest/
cmake . -DCMAKE_INSTALL_PREFIX=$1/googletest-lib -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpicxx
make install -j
cd $1 
rm -rf /googletest;


# sed -i '/GTEST_/d' make_arch/make.lm3_gtest
# echo "GTEST_LIB := $1/googletest-lib/lib64/" >> make_arch/make.lm3_gtest
# echo "GTEST_INC := $1/googletest-lib/include/" >> make_arch/make.lm3_gtest

# Export the variables needed for the gtest library
export FLUPS_GTEST_INC=$1/googletest-lib/include/
export FLUPS_GTEST_LIB=$1/googletest-lib/lib64/


# Compile the tests 
cd $1/test/
MPICC=${FLUPS_MPICC} MPICXX=${FLUPS_MPICXX} \
CXXFLAGS=${FLUPS_CXXFLAGS} CCFLAGS=${FLUPS_CCFLAGS} LDFLAGS=${FLUPS_LDFLAGS} \
HDF5_DIR=${HDF5_DIR} FFTW_DIR=${EBROOTFFTW} \
GTEST_INC=${FLUPS_GTEST_INC} GTEST_LIB=${FLUPS_GTEST_LIB} \ 
make -j 

FILE1=flups_test_a2a
if test -f "$FILE1"; then
    echo "$FILE1 exists."
    exit 0
else
    echo "$FILE1 does not exist."
    exit 1
fi