#!/bin/bash -l
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=xplore_batch_send-${MPI_VERSION}-${TAG}

#-------------------------------------------------------------------------------
## Ceation of the scratch directory
export SCRATCH_DIR=${BASE_SCRATCHDIR}/${SUBMISSION_NAME}  
echo "scratch directory = ${SCRATCH_DIR}"

#-------------------------------------------------------------------------------
## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
export SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/${CLUSTER}_modules.sh

##-------------------------------------------------------------------------------------------------------------
## Go to the scratch directory and copy what's needed
export COMPILE_OPTIONS=(" " "-DMPI_DEFAULT_ORDER")
export COMPILE_PREFIXES=("priority_order" "default_order")
export SCRATCH_DIR_LIST='priority_order default_order'

#batch_length=(1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192)
#export SCRATCH_DIR_LIST='BATCH_SEND_1 BATCH_SEND_2 BATCH_SEND_4 BATCH_SEND_8 BATCH_SEND_16 BATCH_SEND_32 BATCH_SEND_64 BATCH_SEND_128 BATCH_SEND_256 BATCH_SEND_512 BATCH_SEND_1024 BATCH_SEND_2048 BATCH_SEND_4096 BATCH_SEND_8192'

batch_length=(1 2147483647)
export SCRATCH_DIR_LIST='BATCH_SEND_1 BATCH_SEND_2147483647'

for idx in "${!batch_length[@]}";
do
    echo "------ Copying what's needed ..."
    export H3LPR_DIR=${SCRATCH_DIR}/BATCH_SEND_${batch_length[$idx]}/h3lpr/
    export FLUPS_DIR=${SCRATCH_DIR}/BATCH_SEND_${batch_length[$idx]}/flups/

    mkdir -p ${H3LPR_DIR}
    mkdir -p ${FLUPS_DIR}


    cd ${SCRATCH_DIR}
    rsync -r ${HOME_FLUPS} ${FLUPS_DIR}
    rsync -r ${HOME_H3LPR} ${H3LPR_DIR}
    echo " ------ ... done ! "

    #-------------------------------------------------------------------------------
    ## Launch the compilations
    echo " ------ Rewriting the options ..."
    export FLUPS_OPTS="-DMPI_BATCH_SEND=${batch_length[$idx]}"

    echo " ------ Compiling with ${FLUPS_OPTS} -- Executable will be there ${FLUPS_DIR}/samples/validation/"

    echo " ------ Compiling librairies ..."
    COMPILEJOB_ID=$(sbatch --parsable \
                           --job-name=flups_compile_MPI_${MPI_VERSION} \
                           --account=${ACCOUNT} \
                           --partition=${PARTITION} \
                           ${COMPILE_CLUSTER_SPEC} \
                           --nodes=${COMPILE_NNODE} \
                           --ntasks=${COMPILE_NTASK} \
                           --time=${COMPILE_TIME} \
                           ${FLUPS_DIR}/samples/validation/run/benchmark_compile.sh) 
    echo " ------ ... done ! "
done

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS

## 1 Node == 128 CPUS
export NPROC_X=4
export NPROC_Y=4
export NPROC_Z=8

export NRES=1

echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for i in {1..7}
do
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/ ($NPROC_NODES) ))
    
    #---------------------------------------------------------------------------
    export NGLOB_X=$(( ($NPROC_X)*($NPCPUS) ))
    export NGLOB_Y=$(( ($NPROC_Y)*($NPCPUS) ))
    export NGLOB_Z=$(( ($NPROC_Z)*($NPCPUS) ))
    export L_X=$(( ($NPROC_X) ))
    export L_Y=$(( ($NPROC_Y) ))
    export L_Z=$(( ($NPROC_Z) ))
    
    #---------------------------------------------------------------------------
    # Loop on the provided version 
    #---------------------------------------------------------------------------
    export MYNAME=flups_MPI${MPI_VERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
    echo "    sbatch -d afterok:${COMPILEJOB_ID} \
           --job-name=${MYNAME} \
           --account=${ACCOUNT} \
           ${KERNEL_CLUSTER_SPEC} \
           --partition=${PARTITION} \
           --nodes=${NNODE} \
           --ntasks-per-node=${NPROC_NODES} \
           --time=${KERNEL_TIME} \
           ${FLUPS_DIR}/samples/validation/run/benchmark_kernel_xplore.sh"
    echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    
    sbatch -d afterok:${COMPILEJOB_ID} \
	       --job-name=${MYNAME} \
           --account=${ACCOUNT} \
           ${KERNEL_CLUSTER_SPEC} \
           --partition=${PARTITION} \
           --nodes=${NNODE} \
           --ntasks-per-node=${NPROC_NODES} \
           --time=${KERNEL_TIME} \
           ${FLUPS_DIR}/samples/validation/run/benchmark_kernel_xplore.sh
    #---------------------------------------------------------------------------
    
    if [ $(($i%3)) -eq 0 ]
    then
        NPROC_Z=$((2*$NPROC_Z))
    fi
    if [ $((($i)%3)) -eq 1 ]
    then
        NPROC_Y=$((2*$NPROC_Y))
    fi
    if [ $(($i%3)) -eq 2 ]
    then
        NPROC_X=$((2*$NPROC_X))
    fi
done 


