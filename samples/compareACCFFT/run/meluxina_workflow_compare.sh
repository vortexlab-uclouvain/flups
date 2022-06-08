#!/bin/bash -l
## This script is used to submit weak scaling job on the meluxina cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

OMPIVERSION=4.1.3
#CODE_VERSION='dprec_nb dprec_a2a nb a2a'
CODE_VERSION='nb isr'

##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=flups_vs_accfft_${TAG}

## Ceation of the scratch directory
SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  
echo "scratch file = ${SCRATCH_DIR}"

HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/

## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/compareACCFFT/run/meluxina_modules.sh

##-------------------------------------------------------------------------------------------------------------
## Go to the scratch directory and copy what's needed
echo "------ Copying what's needed ..."
export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
export FLUPS_DIR=${SCRATCH_DIR}/flups/
export ACCFFT_DIR=${HOME}/lib-OpenMPI-${OMPIVERSION}
mkdir -p ${H3LPR_DIR}
mkdir -p ${FLUPS_DIR}

cd ${SCRATCH_DIR}
rsync -r ${HOME_FLUPS} ${FLUPS_DIR}
rsync -r ${HOME_H3LPR} ${H3LPR_DIR}
echo " ------ ... done ! "


## Launch the compilations
echo " ------ Compiling librairies ..."
COMPILEJOB_ID=$(sbatch --parsable \
                       --job-name=compile \
                       --export=ALL,OMPIVERSION=${OMPIVERSION},MODULES=${SCRIPT_MODULE} \
                       ${FLUPS_DIR}/samples/compareACCFFT/run/meluxina_compile.sh) 
echo " ------ ... done ! "

##-------------------------------------------------------------------------------------------------------------
## LAUNCH THE JOBS

## 1 Node == 128 CPUS -> 1*64 = 64
export NPROC_X=1
export NPROC_Y=16
export NPROC_Z=16
#export NPCPU=16
npcpu_list=(48)

echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for i in {0..1}
do
for npcpu in ${npcpu_list[@]}
do
    export NGLOB_X=$(( ${NPROC_Z}* ${npcpu} )) # $(( ${NPROC_X}* ${NPCPU} ))
    export NGLOB_Y=$(( ${NPROC_Z}* ${npcpu} )) # $(( ${NPROC_Y}* ${NPCPU} ))
    export NGLOB_Z=$(( ${NPROC_Z}* ${npcpu} )) # $(( ${NPROC_Z}* ${NPCPU} ))
    export L_X=1 #$(( ${NPROC_X} ))
    export L_Y=1 #$(( ${NPROC_Y} ))
    export L_Z=1 #$(( ${NPROC_Z} ))
    export NNODE=$(( (${NPROC_X} * ${NPROC_Y} * ${NPROC_Z})/128 ))
    
    # -------------------------------------------
    # Loop on the provided version 
    # -------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_vs_accfft_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}_NPCPU${npcpu}/
        export MYNAME=flups_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}_NPCPU${npcpu}
        export MODULES=${SCRIPT_MODULE}
        export OMPIVERSION=${OMPIVERSION}
        mkdir -p ${SCRATCH_FLUPS}/prof/
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/compareACCFFT/run/meluxina_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/compareACCFFT/run/meluxina_kernel_valid.sh
    done
done
    NPROC_Y=$((2*$NPROC_Y))
    NPROC_Z=$((2*$NPROC_Z))
    #if [ $(($i%2)) -eq 0 ]
    #then
    #    NPROC_Z=$((2*$NPROC_Z))
    #fi
    #if [ $((($i)%2)) -eq 1 ]
    #then
    #    NPROC_Y=$((2*$NPROC_Y))
    #fi
done 


