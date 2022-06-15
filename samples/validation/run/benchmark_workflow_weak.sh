#!/bin/sh
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=weak_scaling_flups-${MPI_VERSION}-${TAG}

#-------------------------------------------------------------------------------
## Ceation of the scratch directory
SCRATCH_DIR=${BASE_SCRATCHDIR}/${SUBMISSION_NAME}  
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
## Compile the librairies  
export H3LPR_CXXFLAGS='-g -O3 -march=native -fopenmp -DNDEBUG'
export H3LPR_LDFLAGS='-fopenmp'

export FLUPS_CCFLAGS='-g -O3 -std=c99 -DNDEBUG'
export FLUPS_CXXFLAGS='-g -O3 -std=c++17 -DNDEBUG'
export FLUPS_LDFLAGS='-fopenmp '
export FLUPS_OPTS=''

export FFTW_DIR
export HDF5_DIR
## Launch the compilations
echo " ------ Compiling librairies ..."
COMPILEJOB_ID=$(sbatch --parsable \
                       --job-name=flups_compile_CRAYMPICH${MPI_VERSION} \
                       ${FLUPS_DIR}/samples/validation/run/${CLUSTER}_compile.sh) 
echo " ------ ... done ! "

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS

## 1 Node == 128 CPUS
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
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/128 ))
    
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
    export MYNAME=flups_${version}_MPI${MPI_VERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
    echo "sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/${CLUSTER}_kernel_valid.sh "
    echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/${CLUSTER}_kernel_valid.sh
    #---------------------------------------------------------------------------
done 


