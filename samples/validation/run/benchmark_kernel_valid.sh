#!/bin/bash -l
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ${SCRIPT_MODULE} ${MPI_VERSION}

for version in ${CODE_VERSION}
do
    export EXEC_FLUPS=flups_validation_${version}
    export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
    mkdir -p ${SCRATCH_FLUPS}/prof/

    cd ${SCRATCH_FLUPS}
    cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

    echo "----------------- launching job -----------------"
    if [[ ${LCOMMAND} == "srun" ]]; then
        echo "srun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --L=${L_X},${L_Y},${L_Z} --nres=1 --ns=20 --kernel=0"; 
        export OMPI_MCA_coll_hcoll_enable=0
        export OMPI_MCA_coll_hcoll_priority=0
        export OMPI_MCA_coll_base_verbose=10
        OMP_NUM_THREADS=1 srun  -l -u ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=3,3,3,3,3,3 --nres=1 --nsolve=100 --kernel=0
    elif [[ ${LCOMMAND} == "mpirun" ]]; then 
        echo "mpirun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --L=${L_X},${L_Y},${L_Z} --nres=1 --ns=20 --kernel=0"; 
        OMP_NUM_THREADS=1 mpirun ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=3,3,3,3,3,3 --nres=1 --nsolve=100 --kernel=0
    fi 
    echo "----------------- done           -----------------"

    cd -
done