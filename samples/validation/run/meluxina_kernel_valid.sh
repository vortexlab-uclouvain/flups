#!/bin/bash -l
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --ntasks-per-node=256
#SBATCH --time=00:30:00
#SBATCH --output=flups_OpenMPI5.0.0.%j.out
#SBATCH --error=flups_OpenMPI5.0.0.%j.err

#--------------------------------------------------------
# module use /apps/USE/easybuild/staging/2021.1/modules/all
# ## First load the librairies relying on the old version of OpenMPI
# module load HDF5/1.12.1-gompi-2021a
# module load FFTW/3.3.10-gompi-2021a 
# ## Then load the latest version of OpenMPI
# #module load OpenMPI/5.0.0-GCC-10.3.0
# module load OpenMPI/4.1.3-GCC-10.3.0
# #module load OpenMPI/4.1.1-GCC-10.3.0
source ${MODULES} ${OMPIVERSION}

cd ${SCRATCH_FLUPS}
cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

echo "----------------- launching job -----------------"
echo "srun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --L=${L_X},${L_Y},${L_Z} --nres=1 --ns=20 --kernel=0"
#srun ${EXEC_FLUPS} -np ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -res ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -L ${L_X} ${L_Y} ${L_Z} -bc 3 3 3 3 3 3 -nres 1 -ns 20 -k 0
srun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --nsolve=20 --kernel=0
cd -
