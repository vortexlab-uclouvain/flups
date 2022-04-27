#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=flups_auto_test
#SBATCH --time=2:30:00 # hh:mm:ss
#
#SBATCH --ntasks=64
#SBATCH --mem-per-cpu=4000 # megabytes 
#SBATCH --partition=debug
#
#SBATCH --comment=flups
#SBATCH --profile=all

BIN_DIR=/home/users/v/o/vortexbot/flups_test/flups_auto_test/ 
SCRATCH=.

#EXEC=flups_test_nb

#------------------------------------------------------------------------------
module purge

module load OpenMPI/4.0.5-GCC-10.2.0
module load FFTW/3.3.8-gompi-2020b
module load HDF5/1.10.7-gompi-2020b
module load CMake/3.18.4-GCCcore-10.2.0
#------------------------------------------------------------------------------

echo "--------------------------------------------------------------------"
echo "      WELCOME ON LM3!! "
echo "    > runing ${EXEC} "
echo "    >looking here: ${BIN_DIR}/${EXEC}"
echo "    > scratch dir = ${SCRATCH}"
echo "--------------------------------------------------------------------"

#Go to the scratch dir
cd ${SCRATCH}

#------------------------------------------------------------------------------
# run the simulation
mpirun -n ${SLURM_NTASKS} ${EXEC}  --gtest_output="${REPORT}" --gtest_filter=${TESTS}/* > ${TYPE}_$SLURM_JOB_ID
