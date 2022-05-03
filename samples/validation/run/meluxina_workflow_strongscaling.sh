#!/bin/sh
## This script is used to submit strong scaling job on the meluxina cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

OMPIVERSION=4.1.3
# CODE_VERSION='nb a2a'
CODE_VERSION='nb'

##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=CONSTANT_NNODES_OMPI_${OMPIVERSION}_${TAG}
## Ceation of the scratch directory
SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  
echo "scratch file = ${SCRATCH_DIR}"

HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/

## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/meluxina_modules.sh

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

## Launch the compilations of the librairies (H3LPR and FLUPS)
echo " ------ Compiling librairies ..."
COMPILEJOB_ID=$(sbatch --parsable \
                       --job-name=flups_compile_OMPI${OMPIVERSION} \
                       --export=ALL,OMPIVERSION=${OMPIVERSION},MODULES=${SCRIPT_MODULE} \
                       ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 
echo " ------ ... done ! "

## Create what's needed for the scaling
echo " ------ Creating directories ..."
for version in ${CODE_VERSION}
do
  echo " directory: ${SCRATCH_DIR}/simulations_${version}/prof created! "
  mkdir -p ${SCRATCH_DIR}/simulations_${version}/prof 
done 
echo " ------ ... done ! "


##-------------------------------------------------------------------------------------------------------------
## LAUNCH THE JOBS
### Strong scaling information
export NPROC_X=16
export NPROC_Y=16
export NPROC_Z=16

export NNODE=16


echo " ------ Submitting Job scripts"
## Loop on the number of point we want per direction
# ntot_dir=(4 8 12 16 24 32 44 64)
ntot_dir=(4 8 12 16 24 32)
for n in ${ntot_dir[@]}
do  
    export NGLOB_X=$(echo "$(( $n*64 ))")
    export NGLOB_Y=$(echo "$(( $n*64 ))")
    export NGLOB_Z=$(echo "$(( $n*64 ))")
    
    export L_X=$(echo "$(( ($n) ))")
	  export L_Y=$(echo "$(( ($n) ))")
    export L_Z=$(echo "$(( ($n) ))")
    # -------------------------------------------
    # Loop on the provided version 
    # -------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}/
        export MYNAME=flups_${version}_OMPI${OMPIVERSION}_N${NGLOB_X}x${NGLOB_Y}x${NGLOB_Z}
        export MODULES=${SCRIPT_MODULE}
        export OMPIVERSION=${OMPIVERSION}
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
    done
done 
