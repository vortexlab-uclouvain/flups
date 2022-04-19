#!/bin/bash -l
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --ntasks-per-node=256
#SBATCH --time=00:30:00
#SBATCH --output=flups.%j.out
#SBATCH --error=flups.%j.err

#--------------------------------------------------------
module use /apps/USE/easybuild/staging/2021.1/modules/all
module load OpenMPI/4.1.3-GCC-10.3.0 
module load HDF5/1.12.1-gompi-2021a
module load FFTW/3.3.10-gompi-2021a 

cd ${SCRATCH_FLUPS}
cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

echo "----------------- launching job -----------------"
echo "srun ${EXEC_FLUPS} -np ${NPROC_X} ${NPROC_Z} ${NPROC_Z} -res ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0"
srun ${EXEC_FLUPS} -np ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -res ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0
cd -