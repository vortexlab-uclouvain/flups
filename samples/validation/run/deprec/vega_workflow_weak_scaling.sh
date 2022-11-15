#!/bin/sh
## This script is used to submit weak scaling job on the vega cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

OMPIVERSION=4.1.3
CODE_VERSION='a2a'

##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=weak-ompi-${OMPIVERSION}-${TAG}

#-------------------------------------------------------------------------------

## Ceation of the scratch directory
SCRATCH_DIR=/ceph/hpc/scratch/user/euthomasg/${SUBMISSION_NAME}  
HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/
export FFTW_DIR=${PREFIX}
export HDF5_DIR=${PREFIX}

#-------------------------------------------------------------------------------
echo "scratch file = ${SCRATCH_DIR}"

## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/vega_modules.sh

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
                       ${FLUPS_DIR}/samples/validation/run/vega_compile.sh) 
echo " ------ ... done ! "

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS
### Strong scaling information


## 1 Node == 128 CPUS
export NPROC_X=4
export NPROC_Y=4
export NPROC_Z=8

export NPCPU=96

echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for i in {0..0}
do
    #---------------------------------------------------------------------------
    export NGLOB_X=$(( ($NPROC_X)*${NPCPU} ))
    export NGLOB_Y=$(( ($NPROC_Y)*${NPCPU} ))
    export NGLOB_Z=$(( ($NPROC_Z)*${NPCPU} ))
    export L_X=$(( ($NPROC_X) ))
    export L_Y=$(( ($NPROC_Y) ))
    export L_Z=$(( ($NPROC_Z) ))
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/128 ))
    
    #---------------------------------------------------------------------------
    # Loop on the provided version 
    #---------------------------------------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
        export MYNAME=flups_${version}_OMPI${OMPIVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
        export MODULES=${SCRIPT_MODULE}
        export OMPIVERSION=${OMPIVERSION}
        mkdir -p ${SCRATCH_FLUPS}/prof/
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/vega_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/vega_kernel_valid.sh
    done
    
    #---------------------------------------------------------------------------
    # update the proc distribution
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
    #---------------------------------------------------------------------------
done 


