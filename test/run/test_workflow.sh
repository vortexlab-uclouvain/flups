#!/bin/sh
CENTERS='node cell'
VERSIONS='nb a2a isr'
KERNELS='chat2 lgf2 lgf4 lgf6 lgf8 hej2 hej4 hej6 hej8 hej10 hej0 mehr4l mehr6l mehr4f mehr6f'

export LIBDIR=$1

for version in ${VERSIONS}; do
    for center in ${CENTERS}; do
        for kernel in ${KERNELS}; do
            export EXEC=flups_test_${version}
            export REPORT=xml:report_test_${version}_${center}_${kernel}.xml
            export FILE_OUT=std_out_${version}_${center}_${kernel}
            export TESTS=AllTest/ConvergenceTest.AllBoundaryConditions/${center}_${kernel}
            export LIBPATH=${LIBDIR}
            echo "Submitting job with command: sbatch --export=LIBPATH=${LIBDIR},REPORT=${REPORT},FILE_OUT=${FILE_OUT},TESTS=${TESTS},EXEC=${EXEC} test_convergence.sh"
            sbatch test_convergence.sh
            echo "   "
        done
    done
done