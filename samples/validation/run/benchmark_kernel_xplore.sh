#!/bin/bash -l
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ${SCRIPT_MODULE} ${MPI_VERSION}


for version in ${CODE_VERSION}
do
    for scratch_dir in ${SCRATCH_DIR_LIST}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/${scratch_dir}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
        mkdir -p ${SCRATCH_FLUPS}/prof/

        cd ${SCRATCH_FLUPS}
        cp ${SCRATCH_DIR}/${scratch_dir}/flups/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

        echo "---------------- UCX Flags ---------------------"
        echo "UCX_TLS == ${UCX_TLS}"
        echo "OMPI_MCA_pml == ${OMPI_MCA_pml}"
        echo "OMPI_MCA_osc == ${OMPI_MCA_osc}"

        echo "----------------- launching job -----------------"
        echo "Launching jobs with: "
        echo "OMP_NUM_THREADS=1 srun --mpi=pmix --cpu_bind=cores ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=4,4,4,4,4,4 --nres=1 --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}"

       OMP_NUM_THREADS=1 mpirun ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=4,4,4,4,4,4 --nres=${NRES} --nsolve=20 --kernel=${kernel} --center=${CODE_CENTER}
        # if [[ ${LCOMMAND} == "srun" ]]; then
        #     echo "srun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --L=${L_X},${L_Y},${L_Z} --nres=1 --ns=20 --kernel=0"; 
        #     OMP_NUM_THREADS=1 srun  -l -u ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=4,4,4,4,4,4 --nres=1 --nsolve=20 --kernel=0
        # elif [[ ${LCOMMAND} == "mpirun" ]]; then 
        #     echo "mpirun ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --L=${L_X},${L_Y},${L_Z} --nres=1 --ns=20 --kernel=0"; 
        #     OMP_NUM_THREADS=1 mpirun ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=4,4,4,4,4,4 --nres=1 --nsolve=20 --kernel=0
        # fi 
        echo "----------------- done           -----------------"

        cd -
    done 
done
