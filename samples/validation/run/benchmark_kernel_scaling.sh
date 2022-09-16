#!/bin/bash -l
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ${SCRIPT_MODULE} ${MPI_VERSION}

for version in ${CODE_VERSION}
do
    for kernel in ${CODE_KERNEL}
    do 
        export BC_IDX=0 
        for bcs in ${CODE_BCS}
        do
            
            ## Set the name of the directory based on the boundary condition
            if [ $(($BC_IDX)) -eq 0 ]
            then 
                export BC_DIR='unbounded'
            fi 

            if [ $(($BC_IDX)) -eq 1 ]
            then 
                export BC_DIR='periodic'
            fi 

            export EXEC_FLUPS=flups_validation_${version}
            export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_${BC_DIR}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
            mkdir -p ${SCRATCH_FLUPS}/prof/

            cd ${SCRATCH_FLUPS}
            cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

            echo "---------------- UCX Flags ---------------------"
            echo "UCX_TLS == ${UCX_TLS}"
            echo "OMPI_MCA_pml == ${OMPI_MCA_pml}"
            echo "OMPI_MCA_osc == ${OMPI_MCA_osc}"

            echo "----------------- launching job -----------------"
            echo "Launching jobs with: "
            echo "OMP_NUM_THREADS=1 ${LAUNCH_COMMAND} ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=4,4,4,4,4,4 --nres=1 --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}"

            OMP_NUM_THREADS=1 ${LAUNCH_COMMAND} ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=${bcs} --nres=${NRES} --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}
            echo "----------------- done           -----------------"

            cd -

            BC_IDX=$(($BC_IDX + 1))
        done
    done 
done
