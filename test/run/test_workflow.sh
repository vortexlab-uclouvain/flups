#!/bin/sh
KERNELS='CHAT2 LGF2 HEJ2 HEJ4 HEJ6 HEJ8 HEJ10 HEJ0'
CENTERS='Node Cell'
VERSION='nb'

export LIBDIR=$1

for center in ${CENTERS}
do
    for kernel in ${KERNELS}
    do 
        export EXEC_FLUPS=flups_test_${VERSION}
        export REPORT=xml:${center}_${kernel}_report.xml
        export FILE_OUT=${center}_${kernel}
        export TESTS=${center}${kernel}/ConvergenceTest.AllBoundaryConditions
        echo "Submitting job with command:  sbatch --export=LIBPATH=${LIBDIR},REPORT=${REPORT},FILE_OUT=${FILE_OUT},TESTS=${TESTS},EXEC=${EXEC_FLUPS} test_run.sh" 
        sbatch --export=LIBPATH=${LIBDIR},REPORT=${REPORT},FILE_OUT=${FILE_OUT},TESTS=${TESTS},EXEC=${EXEC_FLUPS} test_run.sh
    done
done