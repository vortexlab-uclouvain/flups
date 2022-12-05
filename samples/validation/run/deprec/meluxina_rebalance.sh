#!/bin/bash -l
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=32
#SBATCH --time=04:00:00
#SBATCH --hint=nomultithread

export OMPIVERSION=4.1.3

#-------------------------------------------------------------------------------
# USER SPECIFIC INFO
#-------------------------------------------------------------------------------
CODE_VERSION='nb a2a'
HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/

# indicate wherre the scratch should be put
ROOT_SCRATCH=/project/scratch/p200053

#-------------------------------------------------------------------------------
## create the main dir
#-------------------------------------------------------------------------------
# get a unique id
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SUBMISSION_NAME=rebalance-${TAG}

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
SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/meluxina_modules.sh
source ${SCRIPT_MODULE} ${OMPIVERSION}

#-------------------------------------------------------------------------------
# INT_MAX = 2147483647
config_list=(0 1)
balance_list=(0 1)

# get the number of cpus, one node = 128
# here choosing 4096 cpus = 32 nodes
NPROC_X=16
NPROC_Y=16
NPROC_Z=16
#NPCPU=96

echo "================================================================================="
echo " MAIN SCRATCH DIR = ${SCRATCH_MAIN}"
echo "================================================================================="

for conf in ${config_list[@]}
do
for bl in ${balance_list[@]}
do
    #---------------------------------------------------------------------------
    # with new balance
    if [ "${bl}" -eq "0" ]; then
    	export COMPILE_OPT="-DBALANCE_DPREC";
    else
    	export COMPILE_OPT=" ";
    fi
    export COMPILE_SUFFIX="newbalance${bl}_config${conf}"

    #---------------------------------------------------------------------------
    SCRATCH_DIR=${SCRATCH_MAIN}/${COMPILE_SUFFIX}
    mkdir -p ${SCRATCH_DIR}

    echo "================================================================================="
    echo " BALANCE ${bl} - CONFIG = ${conf}"
    echo " scratch file = ${SCRATCH_DIR}"
    echo "================================================================================="

    #---------------------------------------------------------------------------
    ## Go to the scratch directory and copy what's needed
    export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
    export FLUPS_DIR=${SCRATCH_DIR}/flups/
    mkdir -p ${H3LPR_DIR}
    mkdir -p ${FLUPS_DIR}
    rsync -r ${MAIN_FLUPS} ${FLUPS_DIR}
    rsync -r ${MAIN_H3LPR} ${H3LPR_DIR}

    #---------------------------------------------------------------------------
    # do the compilation 
    ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh
    
    #---------------------------------------------------------------------------
    if [ "${conf}" -eq "0" ]; then
        export NPROC_X=1
        export NPROC_Y=256
        export NPROC_Z=16
        export NGLOB_X=1088
        export NGLOB_Y=1088
        export NGLOB_Z=1088
        export L_X=1
        export L_Y=1
        export L_Z=1
        export NNODE=32 
    elif [ "${conf}" -eq "1" ]; then
        export NPROC_X=1
        export NPROC_Y=128
        export NPROC_Z=32
        export NGLOB_X=1088
        export NGLOB_Y=1088
        export NGLOB_Z=1088
        export L_X=1
        export L_Y=1
        export L_Z=1
        export NNODE=32 
    fi
   
    #---------------------------------------------------------------------------
    # start for the requested versions  - new balance
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}/
        export MYNAME=flups_${version}_OMPI${OMPIVERSION}_N${NPROC_X}x${NPROC_Y}x${NPROC_Z}
	    #---------------------------------------------------------------------
	    mkdir -p ${SCRATCH_FLUPS}
        mkdir -p ${SCRATCH_FLUPS}/prof/
	    cd ${SCRATCH_FLUPS}
	    cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}

	    echo "================================================================================="
	    echo " STARTING JOB - ${version} - balance ${bl} - config = ${conf}"
	    echo "================================================================================="
	    OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --np=${NPROC_X},${NPROC_Y},${NPROC_Z} --res=${NGLOB_X},${NGLOB_Y},${NGLOB_Z} --dom=${L_X},${L_Y},${L_Z} --bc=3,3,3,3,3,3 --nres=1 --nsolve=100 --kernel=0 > ${MYNAME}.log

	    # go back
	    cd -
	    #------------------------------------------------------
    done
    #---------------------------------------------------------------------------
done
done


echo "================================================================================="
echo " MAIN SCRATCH DIR = ${ROOT_SCRATCH}/${SUBMISSION_NAME}"
echo "================================================================================="


