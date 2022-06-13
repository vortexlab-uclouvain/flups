#!/bin/sh
## This script is used to submit weak scaling job on the meluxina cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

OMPIVERSION=4.1.3
CODE_VERSION='a2a nb isr'

##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=weak_scaling_ompi_${OMPIVERSION}_${TAG}

## Ceation of the scratch directory
export SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  
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

## Create what's needed for the scaling
# echo " ------ Creating directories ..."
# for version in ${CODE_VERSION}
# do
#   echo " directory: ${SCRATCH_DIR}/simulations_${version}/prof created! "
#   mkdir -p ${SCRATCH_DIR}/simulations_${version}/prof 
# done 
# echo " ------ ... done ! "

##-------------------------------------------------------------------------------------------------------------
## LAUNCH THE JOBS
### Strong scaling information


## 1 Node == 256 CPUS -> 1*64 = 64
export NPROC_X=4
export NPROC_Y=4
export NPROC_Z=4

echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for i in {0..7}
    do
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


    export NGLOB_X=$(( ($NPROC_X)*96 ))
    export NGLOB_Y=$(( ($NPROC_Y)*96 ))
    export NGLOB_Z=$(( ($NPROC_Z)*96 ))
    export L_X=$(( ($NPROC_X) ))
    export L_Y=$(( ($NPROC_Y) ))
    export L_Z=$(( ($NPROC_Z) ))
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/128 ))
    
    # -------------------------------------------
    # Loop on the provided version 
    # -------------------------------------------
    export CODE_VERSION=${CODE_VERSION}
    export OMPIVERSION=${OMPIVERSION}
    export MODULES=${SCRIPT_MODULE}
    export MYNAME=flups_weak_scaling_OMPI${OMPIVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
    
    echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
    echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
done 


