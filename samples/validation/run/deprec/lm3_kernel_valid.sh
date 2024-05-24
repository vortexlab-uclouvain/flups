#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=NAPS-like
#SBATCH --time=01:00:00
#
##SBATCH --ntasks=${MY_NTASKS}
##SBATCH --cpus-per-task=${MY_NTHREADS}
#SBATCH --mem-per-cpu=2625
#SBATCH --partition=batch 

export OMP_NUM_THREADS=${MY_NTHREADS}

echo "----------------- Load modules -----------------"
module purge
module load intel/2018b
module load HDF5/1.10.2-intel-2018b
module list


#HOME_FLUPS=/home/ucl/tfl/tgillis/flups
EXEC_FLUPS=flups_validation
#
#SCRATCH=$GLOBALSCRATCH/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
#
#mkdir -p $SCRATCH
#mkdir -p $SCRATCH/data
#mkdir -p $SCRATCH/prof
#
#cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
#
#cd $SCRATCH
#srun --label valgrind --leak-check=full --track-origins=yes ./flups
srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZE} -nres 1 -ns 100 -k 0
