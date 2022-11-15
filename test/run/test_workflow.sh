#!/bin/sh
KERNELS='CHAT2 LGF2 HEJ2 HEJ4 HEJ6 HEJ8 HEJ10 HEJ0'
CENTERS='Node Cell'
VERSIONS='nb a2a isr'

export LIBDIR=$1

for version in ${VERSIONS}
do
    export EXEC=flups_test_${version}
    export REPORT=xml:report_test_${version}.xml
    export FILE_OUT=std_out_${version}
    export TESTS=AllTest/ConvergenceTest.AllBoundaryConditions
    export LIBPATH=${LIBDIR}
    echo "Submitting job with command:  sbatch --export=LIBPATH=${LIBDIR},REPORT=${REPORT},FILE_OUT=${FILE_OUT},TESTS=${TESTS},EXEC=${EXEC_FLUPS} test_convergence.sh" 
    sbatch test_convergence.sh
done

