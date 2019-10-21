#!/bin/bash
# Submission script for Marenostrum
#SBATCH --job-name=scaling
#SBATCH --time=00:10:00
#
#SBATCH --output=flups_%j.out
#SBATCH --error=flups_%j.err
#SBATCH --qos=prace

export OMP_NUM_THREADS=${MY_NTHREADS}

echo "----------------- Load modules -----------------"
module purge
module load intel/2018.4
module load impi/2018.4
module load mkl/2018.4
module load hdf5/1.10.5
module load fftw/3.3.6
module list

echo "----------------- launching job -----------------"
echo "launch command: srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${SIZE_X} ${SIZE_Y} ${SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0"

srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${SIZE_X} ${SIZE_Y} ${SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0