#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --account=d2202-040-users
#SBATCH --hint=nomultithread
##SBATCH --exclusive

WISDOM_DIR=${HOME}/flups/fftw-wisdom/wisdom

mkdir -p ${WISDOM_DIR}
fftw-wisdom -v -c -o ${WISDOM_DIR}/vega.wsdm -t 10
