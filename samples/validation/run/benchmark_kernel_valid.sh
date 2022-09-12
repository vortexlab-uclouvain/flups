#!/bin/bash -l
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ${SCRIPT_MODULE} ${MPI_VERSION}


echo "---------------- UCX Flags ---------------------"
echo "UCX_TLS == ${UCX_TLS}"
echo "OMPI_MCA_pml == ${OMPI_MCA_pml}"
echo "OMPI_MCA_osc == ${OMPI_MCA_osc}"

for version in ${CODE_VERSION}
do
    export EXEC_FLUPS=flups_validation_${version}
    export SCRATCH_FLUPS=${SCRATCH_DIR}/flups_${version}
    mkdir -p ${SCRATCH_FLUPS}/prof/
    mkdir -p ${SCRATCH_FLUPS}/data/
    cd ${SCRATCH_FLUPS}
    cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}
    
    for kernel in ${CODE_KERNEL}
    do
        for bcs in ${CODE_BCS}
        do
            echo "----------------- launching job -----------------"
            echo "Launching jobs with: "
            echo "OMP_NUM_THREADS=1 srun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=${bcs} --nres=${NRES} --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}"
            #OMP_NUM_THREADS=1 ${LAUNCH_COMMAND} ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=${bcs} --nres=${NRES} --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}
            OMP_NUM_THREADS=1 ${LAUNCH_COMMAND} ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=${bcs} --nres=${NRES} --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}
            echo "----------------- done           -----------------"  
        done
    done
    cd -
    
done
