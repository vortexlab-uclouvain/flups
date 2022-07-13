#!/bin/bash -l
##-------------------------------------------------------------------------------------------------------------
## BUILD EVERYTHING AND COMPILE
## Compilation has to be done on the compute nodes (the logging node don't have any module information)

## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=xplore_order_list_flups-${MPI_VERSION}-${TAG}

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
export COMPILE_OPTIONS=(" " "-DMPI_DEFAULT_ORDER")
export COMPILE_PREFIXES=("priority_order" "default_order")
for idx in "${!COMPILE_OPTIONS[@]}";
do
    echo "------ Copying what's needed ..."
    export H3LPR_DIR=${SCRATCH_DIR}/${COMPILE_PREFIXES[$idx]}/h3lpr/
    export FLUPS_DIR=${SCRATCH_DIR}/${COMPILE_PREFIXES[$idx]}/flups/

    mkdir -p ${H3LPR_DIR}
    mkdir -p ${FLUPS_DIR}


    cd ${SCRATCH_DIR}
    rsync -r ${HOME_FLUPS} ${FLUPS_DIR}
    rsync -r ${HOME_H3LPR} ${H3LPR_DIR}
    echo " ------ ... done ! "

    #-------------------------------------------------------------------------------
    ## Launch the compilations
    echo " ------ Rewriting the options ..."
    export FLUPS_OPTS=${COMPILE_OPTIONS[$idx]}
    echo " ------ Compiling with ${FLUPS_OPTS} -- Executable will be there ${FLUPS_DIR}/samples/validation/"

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
done
#-------------------------------------------------------------------------------
## LAUNCH THE JOBS

## 1 Node == 128 CPUS
# export NGLOB=(32 64 128 256 512 1024 2048)

# export ARR_NPROC_X=(4 4 4 4 4 8 16)
# export ARR_NPROC_Y=(4 4 4 4 4 8 16)
# export ARR_NPROC_Z=(8 8 8 8 8 8 16)

# export CODE_BCS='4,4,4,4,4,4 
#                  0,0,1,0,3,3 
#                  4,0,4,4,1,4
#                  3,3,4,4,4,4'
# export NRES=1

# echo " ------ Submitting Job scripts"
# # Loop on the number of node needed for the test
# for idx in "${!NGLOB[@]}";
# do    
#     export NPROC_X=${ARR_NPROC_X[$idx]}
#     export NPROC_Y=${ARR_NPROC_Y[$idx]}
#     export NPROC_Z=${ARR_NPROC_Z[$idx]}

#     #---------------------------------------------------------------------------
#     export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/ ($NPROC_NODES) ))
#     export NGLOB_X=${NGLOB[$idx]}
#     export NGLOB_Y=${NGLOB[$idx]}
#     export NGLOB_Z=${NGLOB[$idx]}
#     export L_X=1
#     export L_Y=1
#     export L_Z=1
    
#     #---------------------------------------------------------------------------
#     # Loop on the provided version 
#     #---------------------------------------------------------------------------
#     export MYNAME=flups_MPI${MPI_VERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
#     echo "sbatch -d afterok:${DEPJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/${CLUSTER}_kernel_valid.sh "
#     echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    
#     DEPJOB_ID=$(sbatch --parsable \
#                        -d afterok:${DEPJOB_ID} \
#                        --job-name=${MYNAME} \
#                        --account=${ACCOUNT} \
#                        ${KERNEL_CLUSTER_SPEC} \
#                        --partition=${PARTITION} \
#                        --nodes=${NNODE} \
#                        --ntasks-per-node=${NPROC_NODES} \
#                        --time=${KERNEL_TIME} \
#                        ${FLUPS_DIR}/samples/validation/run/benchmark_kernel_valid.sh)
#     #---------------------------------------------------------------------------
# done 


