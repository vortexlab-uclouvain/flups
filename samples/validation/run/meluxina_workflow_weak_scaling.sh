#!/bin/sh

## ID of the submission -- Needed to differentiate the submissions 
SUBMISSION_NAME=weak_scal_test
#CODE_VERSION='a2a nb dprec_a2a dprec_nb'
CODE_VERSION='dprec_nb'

## Definition of the directories 
SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  
HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/


## Go to the scratch directory and copy what's needed
export H3LPR_DIR=${SCRATCH_DIR}/h3lpr/
export FLUPS_DIR=${SCRATCH_DIR}/flups/
mkdir -p ${H3LPR_DIR}
mkdir -p ${FLUPS_DIR}

cd ${SCRATCH_DIR}
rsync -r --exclude '.git' ${HOME_FLUPS} ${FLUPS_DIR}
rsync -r --exclude '.git' ${HOME_H3LPR} ${H3LPR_DIR}

## Launch the compilations
echo " ------ Compiling librairies ..."
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

export NGLOB_Y=512
export NGLOB_Z=512

export L_Y=2
export L_Z=2

echo " ------ Submitting Job scripts"
## Loop on the number of node needed for the test
i=1
while [ $i -le 2 ]
do
    export NPROC_X=$(echo "$(( ($i)*4  ))")
    export NGLOB_X=$(echo "$(( ($NPROC_X)*64 ))")
    export L_X=$(echo "$(( ($i) ))")
	  
    # -------------------------------------------
    # Loop on the provided version 
    # -------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}/
        export MYNAME=flups_${version}_N${i}
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
    done

    # increment the counter
  	((i++))
done 


