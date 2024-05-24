#!/bin/bash -l
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ${SCRIPT_MODULE} ${MPI_VERSION}


echo "---------------- UCX Flags ---------------------"
echo "UCX_TLS == ${UCX_TLS}"
echo "UCX_DC_MLX5_NUM_DCI == ${UCX_DC_MLX5_NUM_DCI}"
echo "OMPI_MCA_pml == ${OMPI_MCA_pml}"
echo "OMPI_MCA_osc == ${OMPI_MCA_osc}"
echo "CODE_VERSION == ${CODE_VERSION}"
echo "ORDER_DIFF   == ${ORDER_DIFF}"
echo "CODE_CENTER  == ${CODE_CENTER}"

echo "---------------- mpich info ---------------------"
which mpiexec
mpichversion

echo "--------------- ucx info -----------------------"
which ucx_info 
ucx_info -v 

for version in ${CODE_VERSION}
do
    for diff_order in ${ORDER_DIFF}
    do
        export EXEC_FLUPS=flups_tube_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/flups_${version}/order_${diff_order}
        mkdir -p ${SCRATCH_FLUPS}/prof/
        mkdir -p ${SCRATCH_FLUPS}/data/
        cd ${SCRATCH_FLUPS}
        cp ${FLUPS_DIR}/samples/solve_vtube/${EXEC_FLUPS} ${SCRATCH_FLUPS}

        for kernel in ${CODE_KERNEL}
        do
            echo "----------------- launching job -----------------"
            echo "Launching jobs with: "
            echo "OMP_NUM_THREADS=1 srun --mpi=pmix --cpu_bind=cores ${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --case=${FIELD_CASE} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --center=${CODE_CENTER} --kernel=${kernel} --order=${diff_order} ${FIELD_COMPACT}" 
            OMP_NUM_THREADS=1 ${LAUNCH_COMMAND} -np ${SLURM_NTASKS} ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --case=${FIELD_CASE} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --center=${CODE_CENTER} --kernel=${kernel} --order=${diff_order} ${FIELD_COMPACT} 
            echo "----------------- done           -----------------"  
        done
        cd -
    done
done
