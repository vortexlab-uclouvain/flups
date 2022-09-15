#!/bin/bash -l
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=weak_scaling_flups-${MPI_VERSION}-${TAG}

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
echo "------ Copying what's needed ..."
export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
export FLUPS_DIR=${SCRATCH_DIR}/flups/
mkdir -p ${H3LPR_DIR}
mkdir -p ${FLUPS_DIR}

cd ${SCRATCH_DIR}
rsync -r ${HOME_FLUPS} ${FLUPS_DIR}
rsync -r ${HOME_H3LPR} ${H3LPR_DIR}
echo " ------ ... done ! "

#-------------------------------------------------------------------------------
## Launch the compilations
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

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS
export ARR_NPROC_X=(4 4 8 8  8  16 16 16 32 32 32)
export ARR_NPROC_Y=(4 8 8 8  16 16 16 32 32 32 40)
export ARR_NPROC_Z=(8 8 8 16 16 16 32 32 32 48 40)

export NRES=1

export CODE_BCS='4,4,4,4,4,4 
                 3,3,3,3,3,3'

echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for idx in "${!ARR_NPROC_X[@]}";
do    
    export NPROC_X=${ARR_NPROC_X[$idx]}
    export NPROC_Y=${ARR_NPROC_Y[$idx]}
    export NPROC_Z=${ARR_NPROC_Z[$idx]}
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
    echo "sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/${CLUSTER}_kernel_valid.sh "
    echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    
    sbatch -d afterok:${COMPILEJOB_ID} \
           --job-name=${MYNAME} \
           --account=${ACCOUNT} \
           ${KERNEL_CLUSTER_SPEC} \
           --partition=${PARTITION} \
           --nodes=${NNODE} \
           --ntasks-per-node=${NPROC_NODES} \
           --time=${KERNEL_TIME} \
           ${FLUPS_DIR}/samples/validation/run/benchmark_kernel_scaling.sh
    #---------------------------------------------------------------------------
done 


