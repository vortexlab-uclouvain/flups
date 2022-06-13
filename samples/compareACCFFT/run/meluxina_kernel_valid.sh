#!/bin/bash -l
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --ntasks-per-node=128
#SBATCH --time=00:30:00
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

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
cp ${FLUPS_DIR}/samples/compareACCFFT/${EXEC_FLUPS} ${SCRATCH_FLUPS}

echo "----------------- launching job -----------------"
echo "OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --nproc=${NPROC_X},${NPROC_Y},${NPROC_Z} --nglob=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z}"
#OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --nproc=${NPROC_X},${NPROC_Y},${NPROC_Z} --nglob=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --profile
OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --nproc=${NPROC_X},${NPROC_Y},${NPROC_Z} --nglob=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --profile --warm=0
#OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --nproc=${NPROC_X},${NPROC_Y},${NPROC_Z} --nglob=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --warm=0
cd -
