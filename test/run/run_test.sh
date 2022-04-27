#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=flups_auto_test
#SBATCH --time=5:50:00 # hh:mm:ss
#
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4000 # megabytes 
#SBATCH --partition=debug
#
#SBATCH --comment=flups
#SBATCH --profile=all

#------------------------------------------------------------------------------
module purge

module load OpenMPI/4.0.5-GCC-10.2.0
module load FFTW/3.3.8-gompi-2020b
module load HDF5/1.10.7-gompi-2020b
module load CMake/3.18.4-GCCcore-10.2.0
#------------------------------------------------------------------------------

echo "--------------------------------------------------------------------"
echo "      WELCOME ON LM3!! "
echo "    > running ${EXEC} "
echo "    > command =m mpirun -n ${SLURM_NTASKS} ${EXEC}  --gtest_output="${REPORT}" --gtest_filter=${TESTS}/* > ${TYPE}_$SLURM_JOB_ID"
echo "--------------------------------------------------------------------"

#------------------------------------------------------------------------------
# run the simulation
mpirun -n ${SLURM_NTASKS} ${EXEC}  --gtest_output="${REPORT}" --gtest_filter=${TESTS}/* > ${TYPE}_$SLURM_JOB_ID
