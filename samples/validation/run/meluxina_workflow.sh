#!/bin/sh

## ID of the submission -- Needed to differentiate the submissions 
SUBMISSION_NAME=weak_scal_test
CODE_VERSION='a2a nb deprec_a2a deprec_nb'

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
echo " ------ Compiling librairies ..."
export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
export FLUPS_DIR=${SCRATCH_DIR}/flups/
COMPILEJOB_ID=$(sbatch --parsable ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 
echo " ------ ... done ! "

## Create what's needed for the scaling
echo " ------ Creating directories ..."
for version in ${CODE_VERSION}
do
  echo " directory: ${SCRATCH_DIR}/simulations_${version}/prof created! "
  mkdir -p ${SCRATCH_DIR}/simulations_${version}/prof 
done 
echo " ------ ... done ! "



## 1 Node == 256 CPUS -> 1*64 = 64
export NPROC_Y=8
export NPROC_Z=8

export NGLOB_Y=32
export NGLOB_Z=32

export L_Y=2
export L_Z=2

echo " ------ Submitting Job scripts"
## Loop on the number of node needed for the test
i=1
while [ $i -le 1 ]
do
    export NPROC_X=$(echo "$(( ($i)*4  ))")
    export NGLOB_X=$(echo "$(( ($i)*32 ))")
    export L_X=$(echo "$(( ($i) ))")
	  
    # -------------------------------------------
    # All 2 All 
    # -------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}/
        export MYNAME=flups_${version}_N${i}
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
    done
    # # -------------------------------------------
    # # Non Blocking 
    # # ------------------------------------------- 
    # export EXEC_FLUPS=flups_validation_nb
    # export SCRATCH_FLUPS=${SCRATCH_NB}
    # export MYNAME=flups_nb_N${i}
    # echo "Submitting job ${MYNAME} -- Nnodes=${i} "
    # sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh


    # # -------------------------------------------
    # # Non Blocking Aggressive 
    # # ------------------------------------------- 
    # export EXEC_FLUPS=flups_validation_agg
    # export SCRATCH_FLUPS=${SCRATCH_AGG}
    # export MYNAME=flups_agg_N${i}
    # echo "Submitting job ${MYNAME} -- Nnodes=${i} "
    # sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh

    # increment the counter
  	((i++))
done 


