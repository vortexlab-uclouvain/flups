#!/bin/sh
## This script is used to submit weak scaling job on the vega cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

CMPICHVERSION=8.1.12
#CODE_VERSION='a2a dprec_a2a nb dprec_nb'
CODE_VERSION='a2a dprec_a2a nb'

##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=debug_flups-${CMPICHVERSION}-${TAG}

#-------------------------------------------------------------------------------

## Ceation of the scratch directory
SCRATCH_DIR=/scratch/project_465000098/${SUBMISSION_NAME}  
HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/
export FFTW_DIR=${PREFIX}
export HDF5_DIR=${PREFIX}

#-------------------------------------------------------------------------------
echo "scratch file = ${SCRATCH_DIR}"

## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/lumi_modules.sh

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
                       --job-name=flups_compile_CRAYMPICH${CMPICHVERSION} \
                       --export=ALL,CMPICH_VERSION=${CMPICHVERSION},MODULES=${SCRIPT_MODULE} \
                       ${FLUPS_DIR}/samples/validation/run/lumi_compile.sh) 
echo " ------ ... done ! "

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS
### Strong scaling information


## 1 Node == 128 CPUS
export NPROC_X=128
export NPROC_Y=2
export NPROC_Z=1

echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for i in {0..0}
do
    #---------------------------------------------------------------------------
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/128 ))
    export NGLOB_X=384
    export NGLOB_Y=384
    export NGLOB_Z=1536
    export L_X=1
    export L_Y=1
    export L_Z=4
    export FLUPS_BC="4,4,4,4,3,3"
    
    #---------------------------------------------------------------------------
    # Loop on the provided version 
    #---------------------------------------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
        export MYNAME=flups_${version}_CMPICH${CMPICHVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
        export MODULES=${SCRIPT_MODULE}
        export CMPICH_VERSION=${CMPICHVERSION}
        mkdir -p ${SCRATCH_FLUPS}/prof/
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/vega_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/lumi_kernel_debug.sh
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


