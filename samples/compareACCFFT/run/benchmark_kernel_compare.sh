#!/bin/bash -l
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

source ${SCRIPT_MODULE} ${MPI_VERSION}

for version in ${CODE_VERSION}
do
    export EXEC_FLUPS=flups_vs_accfft_${version}
    export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}_NPCPU${npcpu}/
    mkdir -p ${SCRATCH_FLUPS}/prof/
    cd ${SCRATCH_FLUPS}
    cp ${FLUPS_DIR}/samples/compareACCFFT/${EXEC_FLUPS} ${SCRATCH_FLUPS}
    echo "---------------- UCX Flags ---------------------"
    echo "UCX_TLS == ${UCX_TLS}"
    echo "OMPI_MCA_pml == ${OMPI_MCA_pml}"
    echo "OMPI_MCA_osc == ${OMPI_MCA_osc}"
    echo "----------------- launching job -----------------"
    echo "Launching jobs with: "
    echo "OMP_NUM_THREADS=1 srun --mpi=pmix --cpu_bind=cores ./${EXEC_FLUPS} --nproc=${NPROC_X},${NPROC_Y},${NPROC_Z} --nglob=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --profile --warm=0"
    OMP_NUM_THREADS=1 mpirun ./${EXEC_FLUPS} --nproc=${NPROC_X},${NPROC_Y},${NPROC_Z} --nglob=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --profile --warm=5
    echo "----------------- done           -----------------"
    cd -
done
