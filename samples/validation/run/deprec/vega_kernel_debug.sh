#!/bin/bash -l
#SBATCH --account=d2202-040-users
#SBATCH --partition=cpu
#SBATCH --ntasks-per-node=128
#SBATCH --time=00:30:00
#SBATCH --hint=nomultithread
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#--------------------------------------------------------
source ${MODULES} ${OMPIVERSION}

cd ${SCRATCH_FLUPS}
cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

export UCX_TLS="sm,ud"
export OMPI_MCA_pml="ucx"
export OMPI_MCA_osc="ucx"

echo "----------------- launching job -----------------"

OMP_NUM_THREADS=1 mpirun -n ${SLURM_NTASKS} ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=${FLUPS_BC} --nres=1 --nsolve=20 --kernel=0 --center=1
OMP_NUM_THREADS=1 mpirun -n ${SLURM_NTASKS} ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=${FLUPS_BC} --nres=1 --nsolve=20 --kernel=0 --center=0
cd -
