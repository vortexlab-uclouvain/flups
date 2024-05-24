#!/bin/bash
# Submission script for NIC5 
#SBATCH --job-name=vector-validation
#SBATCH --time=00:15:00
#
##SBATCH --ntasks=${MY_NTASKS}
##SBATCH --cpus-per-task=${MY_NTHREADS}
#SBATCH --mem-per-cpu=2625
#SBATCH --partition=batch

export OMP_NUM_THREADS=${MY_NTHREADS}

echo "----------------- Load modules -----------------"
module purge
module load FFTW/3.3.8-gompi-2020b HDF5/1.10.7-gompi-2020b
module list


#HOME_FLUPS=/home/ucl/tfl/tgillis/flups
EXEC_FLUPS=flups_validation_nb
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
mpirun ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZEX} ${MY_SIZEY} ${MY_SIZEZ} -nres ${MY_NRES} -ns ${MY_NSOLVE} -k ${MY_KERNEL} -c 0 -bc ${MY_BC}
