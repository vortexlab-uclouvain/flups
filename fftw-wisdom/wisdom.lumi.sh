#!/bin/bash -l
#SBATCH --account=project_465000098
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --hint=nomultithread

module load cray-mpich/8.1.12
module load cray-fftw/3.3.8.12

module list

WISDOM_DIR=${HOME}/flups/fftw-wisdom/wisdom

mkdir -p ${WISDOM_DIR}
fftw-wisdom -v -c -o ${WISDOM_DIR}/lumi.wsdm -t 10
