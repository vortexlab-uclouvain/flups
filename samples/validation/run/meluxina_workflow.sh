#!/bin/sh

## ID of the submission -- Needed to differentiate the submissions 
SUBMISSION_NAME=weak_scal_test

## Definition of the directories 
SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  
HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/


## Go to the scratch directory and copy what's needed
mkdir -p ${SCRATCH_DIR}
cd ${SCRATCH_DIR}
cp -r ${HOME_FLUPS} ${SCRATCH_DIR}
cp -r ${HOME_H3LPR} ${SCRATCH_DIR}

## Launch the compilations
export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
export FLUPS_DIR=${SCRATCH_DIR}/flups/
COMPILEJOB_ID=$(sbatch --parsable ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 

## Create what's needed for the scaling
SCRATCH_A2A=${SCRATCH_DIR}/flups_a2a_simulations/
mkdir -p ${SCRATCH_A2A}/prof 

SCRATCH_NB=${SCRATCH_DIR}/flups_nb_simulations/
mkdir -p ${SCRATCH_NB}/prof 

SCRATCH_AGG=${SCRATCH_DIR}/flups_agg_simulations/
mkdir -p ${SCRATCH_AGG}/prof 


## 1 Node == 256 CPUS -> 1*64 = 64
export NPROC_Y=8
export NPROC_Z=8

export NGLOB_Y=32
export NGLOB_Z=32

export L_Y = 2
export L_Z = 2

i=1
while [ $i -le 3 ]
do
    export NPROC_X=$(echo "$(( ($i)*4  ))")
    export NGLOB_X=$(echo "$(( ($i)*32 ))")
    export L_X=$(echo "$(( ($i) ))")
	  
    # -------------------------------------------
    # All 2 All 
    # ------------------------------------------- 
    export EXEC_FLUPS=flups_validation_a2a
    export SCRATCH_FLUPS=${SCRATCH_A2A}
    export MYNAME=flups_a2a_N${i}
    echo "Submitting job ${MYNAME} -- Nnodes=${i} "
    sbatch -d afterok:COMPILEJOB_ID --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/run/meluxina_kernel_valid.sh

    # -------------------------------------------
    # Non Blocking 
    # ------------------------------------------- 
    export EXEC_FLUPS=flups_validation_nb
    export SCRATCH_FLUPS=${SCRATCH_NB}
    export MYNAME=flups_nb_N${i}
    echo "Submitting job ${MYNAME} -- Nnodes=${i} "
    sbatch -d afterok:COMPILEJOB_ID --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/run/meluxina_kernel_valid.sh


    # -------------------------------------------
    # Non Blocking Aggressive 
    # ------------------------------------------- 
    export EXEC_FLUPS=flups_validation_agg
    export SCRATCH_FLUPS=${SCRATCH_AGG}
    export MYNAME=flups_agg_N${i}
    echo "Submitting job ${MYNAME} -- Nnodes=${i} "
    sbatch -d afterok:COMPILEJOB_ID --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/run/meluxina_kernel_valid.sh

    # increment the counter
  	((i++))
done 


