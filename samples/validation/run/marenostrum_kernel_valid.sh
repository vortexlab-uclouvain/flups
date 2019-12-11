#!/bin/bash
# Submission script for Marenostrum
#SBATCH --job-name=scaling
#
#SBATCH --output=flups_%j.out
#SBATCH --error=flups_%j.err
#SBATCH --qos=prace
#SBATCH --exclude=s07r2b[01-24]
#--> one of these failed at FFTW plans alloc
#SBATCH --exclude=s05r1b[01-24]
#--> proc s05r1b16 gave invalid address or slot during writev

export OMP_NUM_THREADS=${MY_NTHREADS}

echo "----------------- Load modules -----------------"
module purge
# module load intel/2018.4
# module load impi/2018.4
# module load mkl/2018.4
module laod intel/2017.4
module load impi/2017.4

module load hdf5/1.10.5
module load fftw/3.3.6
module list

echo "----------------- launching job -----------------"
echo "launch command: srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${SIZE_X} ${SIZE_Y} ${SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0"

srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${SIZE_X} ${SIZE_Y} ${SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0

scontrol show job ${SLURM_JOB_ID}