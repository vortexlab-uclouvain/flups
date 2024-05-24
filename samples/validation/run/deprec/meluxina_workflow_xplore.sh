#!/bin/sh
## This script is used to submit weak scaling job on the meluxina cluster. 
## In the first part of the script, you have to choose the version of OpenMPI you want to use 
## as well as the code version (all to all or non-blocking). The running directory is then named 
## after the version of OMPI

CODE_VERSION='nb'

#-------------------------------------------------------------------------------------------------------------
# CLUSTER SPECIFIC INFO
#-------------------------------------------------------------------------------------------------------------
export OMPIVERSION=4.1.3
export HOME_FLUPS=${HOME}/flups/
export HOME_H3LPR=${HOME}/h3lpr/
#export FFTW_DIR=${EBROOTFFTW}
#export HDF5_DIR=${EBROOTHDF5}

# indicate wherre the scratch should be put
ROOT_SCRATCH=/project/scratch/p200053

#-------------------------------------------------------------------------------------------------------------
## Definition of the directories 
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=xplore-ompi-${OMPIVERSION}-${TAG}

#-------------------------------------------------------------------------------
## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/meluxina_modules.sh


#-------------------------------------------------------------------------------
# INT_MAX = 2147483647
#batch_length=(1 4 16 32 64 128 256 512 2147483647)
batch_length=(1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192)

# get the number of cpus, one node = 128
export NPROC_X=16
export NPROC_Y=16
export NPROC_Z=16
# here choosing 4096 cpus = 32 nodes

export NPCPU=96

echo "================================================================================="
echo " MAIN SCRATCH DIR = ${ROOT_SCRATCH}/${SUBMISSION_NAME}"
echo "================================================================================="

#-------------------------------------------------------------------------------
## Launch the compilations
echo " ------ Compiling librairies ..."
export MODULES=${SCRIPT_MODULE}

for bl in ${batch_length[@]}
do
    #---------------------------------------------------------------------------
    # with new balance
    export COMPILE_OPT="-DMPI_BATCH_SEND=${bl}"
    export COMPILE_SUFFIX="batch${bl}"

    #---------------------------------------------------------------------------
    SCRATCH_DIR=${ROOT_SCRATCH}/${SUBMISSION_NAME}/${COMPILE_SUFFIX}
    echo "scratch file = ${SCRATCH_DIR}"

    #---------------------------------------------------------------------------
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

    #---------------------------------------------------------------------------
    COMPILEJOB_ID=$(sbatch --parsable \
                       --job-name=flups_compile_OMPI${OMPIVERSION} \
                       ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 
    
    ## with old balance
    #export COMPILE_OPT="-DMPI_BATCH_SEND=${bl} -DBALANCE_DPREC"
    #export COMPILE_SUFFIX="batch${bl}_bal0"
    #COMPILEJOB_ID=$(sbatch --parsable \
    #                   --job-name=flups_compile_OMPI${OMPIVERSION} \
    #                   ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 

    #---------------------------------------------------------------------------
    export NGLOB_X=$(( ($NPROC_X)*${NPCPU} ))
    export NGLOB_Y=$(( ($NPROC_Y)*${NPCPU} ))
    export NGLOB_Z=$(( ($NPROC_Z)*${NPCPU} ))
    export L_X=$(( ($NPROC_X) ))
    export L_Y=$(( ($NPROC_Y) ))
    export L_Z=$(( ($NPROC_Z) ))
    export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/128 ))
   
    #---------------------------------------------------------------------------
    # start for the requested versions  - new balance
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
        export MYNAME=flups_${version}_OMPI${OMPIVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
        export MODULES=${SCRIPT_MODULE}
        export OMPIVERSION=${OMPIVERSION}
        mkdir -p ${SCRATCH_FLUPS}/prof/
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
    done

    ##---------------------------------------------------------------------------
    ## start for the requested versions  - old balance
    #for version in ${CODE_VERSION}
    #do 
    #    export EXEC_FLUPS=flups_validation_${version}
    #    export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
    #    export MYNAME=flups_${version}_OMPI${OMPIVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
    #    export MODULES=${SCRIPT_MODULE}
    #    export OMPIVERSION=${OMPIVERSION}
    #    mkdir -p ${SCRATCH_FLUPS}/prof/
    #    echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
    #    echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
    #    sbatch -d afterok:${COMPILEJOB_ID} --nodes=${NNODE} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
    #done
    
    
    #---------------------------------------------------------------------------
done
echo " ------ ... done ! "

echo "================================================================================="
echo " MAIN SCRATCH DIR = ${ROOT_SCRATCH}/${SUBMISSION_NAME}"
echo "================================================================================="


