#!/bin/bash -l
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=compare_flups_accfft-${MPI_VERSION}-${TAG}

#-------------------------------------------------------------------------------
## Ceation of the scratch directory
export SCRATCH_DIR=${BASE_SCRATCHDIR}/${SUBMISSION_NAME}  
echo "scratch directory = ${SCRATCH_DIR}"

#-------------------------------------------------------------------------------
## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
export SCRIPT_MODULE=${HOME_FLUPS}/samples/compareACCFFT/run/${CLUSTER}_modules.sh

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
                       ${FLUPS_DIR}/samples/compareACCFFT/run/benchmark_compile.sh) 
echo " ------ ... done ! "

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS

## 1 Node == 128 CPUS
export NPROC_X=1
export NPROC_Y=16
export NPROC_Z=16

export L_X=1 #$(( ${NPROC_X} ))
export L_Y=1 #$(( ${NPROC_Y} ))
export L_Z=1 #$(( ${NPROC_Z} ))

npcpu_list=(48)


export NGLOB_X=$(( ${NPROC_Z}*48 )) # $(( ${NPROC_X}* ${NPCPU} ))

echo " ------ Submitting Job scripts"

# Loop on the number of node needed for the test
for i in {0..2}
do
    for npcpu in ${npcpu_list[@]}
    do
        export NGLOB_Y=$(( ${NPROC_Y}* ${npcpu} )) # $(( ${NPROC_Y}* ${NPCPU} ))
        export NGLOB_Z=$(( ${NPROC_Z}* ${npcpu} )) # $(( ${NPROC_Z}* ${NPCPU} ))
        export NNODE=$(( (${NPROC_X} * ${NPROC_Y} * ${NPROC_Z})/128 ))

        export NPCPU=${npcpu}
        #---------------------------------------------------------------------------
        # Loop on the provided version 
        #---------------------------------------------------------------------------
        
        export MYNAME=flups_MPI${MPI_VERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}_NPCPU${npcpu}
        echo "sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/compareACCFFT/run/${CLUSTER}_kernel_compare.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} \
               --job-name=${MYNAME} \
               --account=${ACCOUNT} \
               ${KERNEL_CLUSTER_SPEC} \
               --partition=${PARTITION} \
               --nodes=${NNODE} \
               --ntasks-per-node=${NPROC_NODES} \
               --time=${KERNEL_TIME} \
               ${FLUPS_DIR}/samples/compareACCFFT/run/benchmark_kernel_compare.sh
        #---------------------------------------------------------------------------
        
        NPROC_Y=$((2*$NPROC_Y))
        NPROC_Z=$((2*$NPROC_Z))
        
        L_Y=$((2*$L_Y))
        L_Z=$((2*$L_Z))
    done
done 


