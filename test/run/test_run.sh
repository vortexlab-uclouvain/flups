#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=flups_auto_test
#SBATCH --time=6:00:00 # hh:mm:ss
#
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=4000 # megabytes 
#SBATCH --partition=batch
#
#SBATCH --comment=flups
#SBATCH --profile=all

#------------------------------------------------------------------------------
module purge

module load releases/2021b
module load OpenMPI/4.1.2-GCC-11.2.0
module load FFTW/3.3.10-gompi-2021b
module load HDF5/1.12.1-gompi-2021b
module load CMake/3.21.1-GCCcore-11.2.0
#------------------------------------------------------------------------------

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${LIBPATH}

echo "--------------------------------------------------------------------"
echo "      WELCOME ON LM3!! "
echo "    > running ${EXEC} "
echo "    > command = mpirun -n ${SLURM_NTASKS} ${EXEC}  --gtest_output="${REPORT}" --gtest_filter=${TESTS}/* > ${FILE_OUT}_$SLURM_JOB_ID"
echo "--------------------------------------------------------------------"

#------------------------------------------------------------------------------
# run the simulation
mpirun -n ${SLURM_NTASKS} ${EXEC} --gtest_output="${REPORT}" --gtest_filter=${TESTS}/* > ${FILE_OUT}_$SLURM_JOB_ID
