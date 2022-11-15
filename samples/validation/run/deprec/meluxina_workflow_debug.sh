#!/bin/sh
## This script is used to submit weak scaling job on the meluxina cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

OMPIVERSION=4.1.3
CODE_VERSION='nb dprec_nb a2a dprec_a2a'

##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=FLUPS_DEBUG_OMPI${OMPIVERSION}_${TAG}

## Ceation of the scratch directory
SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  
echo "scratch file = ${SCRATCH_DIR}"

HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/

## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/meluxina_modules.sh

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


## Launch the compilations
echo " ------ Compiling librairies ..."
COMPILEJOB_ID=$(sbatch --parsable \
                       --job-name=flups_compile_OMPI${OMPIVERSION} \
                       --export=ALL,OMPIVERSION=${OMPIVERSION},MODULES=${SCRIPT_MODULE} \
                       ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 
echo " ------ ... done ! "

##-------------------------------------------------------------------------------------------------------------
## LAUNCH THE JOB

## 1 Node == 128 CPUS 
echo " ------ Submitting Job scripts"
for version in ${CODE_VERSION}
do 
    export EXEC_FLUPS=flups_validation_${version}
    export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
    export MYNAME=flups_${version}_OMPI${OMPIVERSION}
    export MODULES=${SCRIPT_MODULE}
    export OMPIVERSION=${OMPIVERSION}
    mkdir -p ${SCRATCH_FLUPS}/prof/
    echo "Submitting job with command: sbatch -d afterok:${COMPILEJOB_ID} --nodes=4 --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_debug.sh "
    sbatch -d afterok:${COMPILEJOB_ID} --nodes=4 --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_debug.sh
done 
