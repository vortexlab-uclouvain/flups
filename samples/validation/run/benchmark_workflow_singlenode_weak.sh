#!/bin/bash -l
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=single_weak_scaling_flups-${MPI_VERSION}-${TAG}

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
## LAUNCH THE JOBS - One a single node 

echo " ------ Submitting Job scripts"
export NNODE=1
    

    
#---------------------------------------------------------------------------
# Loop on the provided version 
#---------------------------------------------------------------------------
export MYNAME=flups_MPI${MPI_VERSION}_single_weak_scaling_flups

sbatch -d afterok:${COMPILEJOB_ID} \
       --job-name=${MYNAME} \
       --account=${ACCOUNT} \
       ${KERNEL_CLUSTER_SPEC} \
       --partition=${PARTITION} \
       --nodes=${NNODE} \
       --ntasks-per-node=${NPROC_NODES} \
       --time=${KERNEL_TIME} \
       ${FLUPS_DIR}/samples/validation/run/benchmark_kernel_singlenode_scaling.sh
#---------------------------------------------------------------------------
