#!/bin/bash
# Submission script for Lemaitre3 

#SBATCH --job-name=Wisdom-generation
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=2625
#SBATCH --partition=batch 
#SBATCH --output=wisdom%j.out
#SBATCH --error=wisdom%j.err

# ---------------------------------------------------------------------------------------
module purge
module load releases/2021b
module load OpenMPI/4.1.2-GCC-11.2.0
module load FFTW/3.3.10-gompi-2021b
module load HDF5/1.12.1-gompi-2021b
module load CMake/3.21.1-GCCcore-11.2.0
# ---------------------------------------------------------------------------------------

WISDOM_DIR=${HOME}/dev-flups/flups/fftw-wisdom/wisdom

mkdir -p ${WISDOM_DIR}
fftw-wisdom -v -c -o ${WISDOM_DIR}/lm3.wsdm -t 10
