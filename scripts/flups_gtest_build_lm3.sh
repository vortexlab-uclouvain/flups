#!/bin/bash
source $2

git clone https://github.com/google/googletest.git $1/googletest/
cd $1/googletest/
cmake . -DCMAKE_INSTALL_PREFIX=$1/googletest-lib -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpicxx
make install -j

cd $1 
rm -rf /googletest;
sed -i '/GTEST_/d' make_arch/make.lm3_gtest
echo "GTEST_LIB := $1/googletest-lib/lib64/" >> make_arch/make.lm3_gtest
echo "GTEST_INC := $1/googletest-lib/include/" >> make_arch/make.lm3_gtest


cd $1/test/
make -j4

FILE1=flups_test_a2a
if test -f "$FILE1"; then
    echo "$FILE1 exists."
    exit 0
else
    echo "$FILE1 does not exist."
    exit 1
fi