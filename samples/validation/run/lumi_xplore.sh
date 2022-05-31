#!/bin/bash -l
#SBATCH --account=project_465000098
#SBATCH --partition=standard
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=128
#SBATCH --time=04:00:00
#SBATCH --hint=nomultithread

export CMPICHVERSION=8.1.12

#-------------------------------------------------------------------------------
# USER SPECIFIC INFO
#-------------------------------------------------------------------------------
CODE_VERSION='nb'
HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/

# indicate wherre the scratch should be put
ROOT_SCRATCH=/scratch/project_465000098/

#-------------------------------------------------------------------------------
## create the main dir
#-------------------------------------------------------------------------------
# get a unique id
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=xplore-${TAG}

#get the scratch
SCRATCH_MAIN=${ROOT_SCRATCH}/${SUBMISSION_NAME}
mkdir -p ${SCRATCH_MAIN}

# get flups + h3lpr
MAIN_H3LPR=${SCRATCH_MAIN}/h3lpr/
MAIN_FLUPS=${SCRATCH_MAIN}/flups/
mkdir -p ${MAIN_H3LPR}
mkdir -p ${MAIN_FLUPS}
rsync -r ${HOME_FLUPS} ${MAIN_FLUPS}
rsync -r ${HOME_H3LPR} ${MAIN_H3LPR}

#-------------------------------------------------------------------------------
## We use a script to load the different module. This script takes as an argument, the version of OMPI choosen at 
## the beginning of the file
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/lumi_modules.sh
source ${SCRIPT_MODULE} ${OMPIVERSION}

#-------------------------------------------------------------------------------
# INT_MAX = 2147483647
size_length=(32 48 64 96 128)
batch_length=(1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192)

# get the number of cpus, one node = 128
# here choosing 512 cpus = 4 nodes
NPROC_X=16
NPROC_Y=16
NPROC_Z=16
#NPCPU=96

echo "================================================================================="
echo " MAIN SCRATCH DIR = ${SCRATCH_MAIN}"
echo "================================================================================="

for NPCPU in ${size_length[@]}
do
    for bl in ${batch_length[@]}
    do
        #-----------------------------------------------------------------------
        # with new balance
        export COMPILE_OPT="-DMPI_BATCH_SEND=${bl}"
        export COMPILE_SUFFIX="batch${bl}_npcpu${NPCPU}"
    
        #-----------------------------------------------------------------------
        SCRATCH_DIR=${SCRATCH_MAIN}/${COMPILE_SUFFIX}
        mkdir -p ${SCRATCH_DIR}
    
        echo "================================================================================="
        echo " BATCH LENGTH = ${bl} - NPCPU = ${NPCPU}"
        echo "scratch file = ${SCRATCH_DIR}"
        echo "================================================================================="
    
        #-----------------------------------------------------------------------
        ## Go to the scratch directory and copy what's needed
        export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
        export FLUPS_DIR=${SCRATCH_DIR}/flups/
        mkdir -p ${H3LPR_DIR}
        mkdir -p ${FLUPS_DIR}
        rsync -r ${MAIN_FLUPS} ${FLUPS_DIR}
        rsync -r ${MAIN_H3LPR} ${H3LPR_DIR}
    
        #-----------------------------------------------------------------------
        # do the compilation 
        ${FLUPS_DIR}/samples/validation/run/lumi_compile.sh
        
        #-----------------------------------------------------------------------
        export NGLOB_X=$(( ($NPROC_X)*${NPCPU} ))
        export NGLOB_Y=$(( ($NPROC_Y)*${NPCPU} ))
        export NGLOB_Z=$(( ($NPROC_Z)*${NPCPU} ))
        export L_X=$(( ($NPROC_X) ))
        export L_Y=$(( ($NPROC_Y) ))
        export L_Z=$(( ($NPROC_Z) ))
        export NNODE=$(( ($NPROC_X * $NPROC_Y * $NPROC_Z)/128 ))
       
        #-----------------------------------------------------------------------
        # start for the requested versions  - new balance
        for version in ${CODE_VERSION}
        do 
            export EXEC_FLUPS=flups_validation_${version}
            export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
            export MYNAME=flups_${version}_OMPI${OMPIVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
    	    #-------------------------------------------------------------------
    	    mkdir -p ${SCRATCH_FLUPS}
            mkdir -p ${SCRATCH_FLUPS}/prof/
    	    cd ${SCRATCH_FLUPS}
    	    cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}
    
    	    echo "================================================================================="
    	    echo " STARTING JOB - ${version} - batch length ${bl} - npcpu = ${NPCPU}"
    	    echo "================================================================================="
    	    OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=3,3,3,3,3,3 --nres=1 --nsolve=100 --kernel=0 > ${MYNAME}.log
    
    	    # go back
    	    cd -
    	    #-------------------------------------------------------------------
        done

        #-----------------------------------------------------------------------
        # removes the flups and h3lpr folder
        rm -rf ${FLUPS_DIR} ${H3LPR_DIR}

        #-----------------------------------------------------------------------
    done
done

#-------------------------------------------------------------------------------
# remove the flups and h3lpr folder
rm -rf ${MAIN_FLUPS} ${MAIN_H3LPR}

#-------------------------------------------------------------------------------
echo "================================================================================="
echo " MAIN SCRATCH DIR = ${ROOT_SCRATCH}/${SUBMISSION_NAME}"
echo "================================================================================="


