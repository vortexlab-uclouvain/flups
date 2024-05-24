#!/bin/bash -l
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=validation_biot_savart_flups-${MPI_VERSION}-${TAG}

#-------------------------------------------------------------------------------
## Ceation of the scratch directory
export SCRATCH_DIR=${BASE_SCRATCHDIR}/${SUBMISSION_NAME}  
echo "scratch directory = ${SCRATCH_DIR}"

#-------------------------------------------------------------------------------
## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
export SCRIPT_MODULE=${HOME_FLUPS}/samples/solve_vtube/run/${CLUSTER}_modules.sh

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
DEPJOB_ID=$(sbatch --parsable \
                       --job-name=flups_compile_MPI_${MPI_VERSION} \
                       --account=${ACCOUNT} \
                       --partition=${PARTITION} \
                       ${COMPILE_CLUSTER_SPEC} \
                       --nodes=${COMPILE_NNODE} \
                       --ntasks=${COMPILE_NTASK} \
                       --time=${COMPILE_TIME} \
                       ${FLUPS_DIR}/samples/solve_vtube/run/benchmark_compile.sh) 
echo " ------ ... done ! "

#-------------------------------------------------------------------------------
## LAUNCH THE JOBS

## 1 Node == 128 CPUS
export NGLOB=(32 64 128 256 512 1024 2048)

## Number of proc per dimension for each case
export ARR_NPROC_X=(4 4 4 4 4 8 16)
export ARR_NPROC_Y=(4 4 4 4 4 8 16)
export ARR_NPROC_Z=(8 8 8 8 8 8 16)


#export NGLOB=(1024 2048)

### Number of proc per dimension for each case
#export ARR_NPROC_X=(16 16)
#export ARR_NPROC_Y=(16 16)
#export ARR_NPROC_Z=(16 32)

## Number of resolution to be tested inside a single job
export NRES=1

## Order of the differentiation
# 1 - spectral like
# 2 - FD2 
# 4 - FD4
# 6 - FD6 
export ORDER_DIFF='1 2 4 6'

## Test case used
# 0 - vortex tube
# 1 - Vortex Ring /!\ the analytical solution is not correct
# 4 - Fully periodic case
export FIELD_CASE='0'

## Is the field compact (options used with the tube or the ring)
#export FIELD_COMPACT='--compact'
export FIELD_COMPACT=''


echo " ------ Submitting Job scripts"
# Loop on the number of node needed for the test
for idx in "${!NGLOB[@]}";
do    
    export NPROC_X=${ARR_NPROC_X[$idx]}
    export NPROC_Y=${ARR_NPROC_Y[$idx]}
    export NPROC_Z=${ARR_NPROC_Z[$idx]}

    #---------------------------------------------------------------------------
    
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/ ($NPROC_NODES) ))
    export NGLOB_X=${NGLOB[$idx]}
    export NGLOB_Y=${NGLOB[$idx]}
    export NGLOB_Z=${NGLOB[$idx]}
    export L_X=1
    export L_Y=1
    export L_Z=1
    
    #---------------------------------------------------------------------------
    # Loop on the provided version 
    #---------------------------------------------------------------------------
    export MYNAME=flups_MPI${MPI_VERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
    echo "sbatch -d afterok:${DEPJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/solve_vtube/run/benchmark_kernel_tube.sh"
    echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    
    DEPJOB_ID=$(sbatch --parsable \
                       -d afterok:${DEPJOB_ID} \
                       --job-name=${MYNAME} \
                       --account=${ACCOUNT} \
                       ${KERNEL_CLUSTER_SPEC} \
                       --partition=${PARTITION} \
                       --nodes=${NNODE} \
                       --ntasks-per-node=${NPROC_NODES} \
                       --time=${KERNEL_TIME} \
                       ${FLUPS_DIR}/samples/solve_vtube/run/benchmark_kernel_tube.sh)
    #---------------------------------------------------------------------------
done 


