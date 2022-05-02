#!/bin/sh

## ID of the submission -- Needed to differentiate the submissions 
OMPIVERSION=4.1.3
CODE_VERSION='nb a2a'

## Definition of the directories 
SUBMISSION_NAME=strong_scaling_OMPI_${OMPIVERSION}

HOME_FLUPS=${HOME}/flups/
HOME_H3LPR=${HOME}/h3lpr/

SCRIPT_MODULE=${HOME_FLUPS}/samples/validation/run/meluxina_modules.sh

SCRATCH_DIR=/project/scratch/p200053/${SUBMISSION_NAME}  


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
COMPILEJOB_ID=$(sbatch --parsable \
                       --export=ALL,OMPIVERSION=${OMPIVERSION},MODULES=${SCRIPT_MODULE} \
                       ${FLUPS_DIR}/samples/validation/run/meluxina_compile.sh) 
echo " ------ ... done ! "

## Create what's needed for the scaling
echo " ------ Creating directories ..."
for version in ${CODE_VERSION}
do
  echo " directory: ${SCRATCH_DIR}/simulations_${version}/prof created! "
  mkdir -p ${SCRATCH_DIR}/simulations_${version}/prof 
done 
echo " ------ ... done ! "



### Strong scaling information
export NPROC_X=16
export NPROC_Y=16
export NPROC_Z=16




echo " ------ Submitting Job scripts"
## Loop on the number of point we want per direction
ntot_dir=(1 2 4 8 12 16 24 32 44 64)
for n in ${ntot_dir[@]}
do  
    export NGLOB_X=$(echo "$(( $n*64 ))")
    export NGLOB_Y=$(echo "$(( $n*64 ))")
    export NGLOB_Z=$(echo "$(( $n*64 ))")
    
    export L_X=$(echo "$(( ($n) ))")
	  export L_Y=$(echo "$(( ($n) ))")
    export L_Z=$(echo "$(( ($n) ))")
    # -------------------------------------------
    # Loop on the provided version 
    # -------------------------------------------
    for version in ${CODE_VERSION}
    do 
        export EXEC_FLUPS=flups_validation_${version}
        export SCRATCH_FLUPS=${SCRATCH_DIR}/simulations_${version}/
        export MYNAME=flups_${version}_N${NGLOB_X}x${NGLOB_Y}x${NGLOB_Z}
        export MODULES=${SCRIPT_MODULE}
        export OMPIVERSION=${OMPIVERSION}
        echo "Submitting job with command:  sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh "
        echo "NGLOB = ${NGLOB_X} ${NGLOB_Y} ${NGLOB_Z} -- NPROC = ${NPROC_X} ${NPROC_Y} ${NPROC_Z} -- L = ${L_X} ${L_Y} ${L_Z}"
        sbatch -d afterok:${COMPILEJOB_ID} --nodes=${i} --job-name=${MYNAME} ${FLUPS_DIR}/samples/validation/run/meluxina_kernel_valid.sh
    done

    # increment the counter
  	((i++))
done 


