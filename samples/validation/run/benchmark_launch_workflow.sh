#!/bin/bash -l
#---------------------------------------------------------------------------------------
#               Case dependent info
#---------------------------------------------------------------------------------------
export CODE_VERSION='dprec_nb dprec_a2a a2a nb isr'
export NPCPUS=96

## Compile the librairies  
export H3LPR_CXXFLAGS='-g -O3 -march=native -fopenmp -DNDEBUG'
export H3LPR_LDFLAGS='-fopenmp -rdynamic -ldl'

export FLUPS_CCFLAGS='-g -O3 -std=c99'
export FLUPS_CXXFLAGS='-g -O3 -std=c++17'
export FLUPS_LDFLAGS='-fopenmp '
export FLUPS_OPTS=''

# ..................................................................
## Compile bash options
export COMPILE_NNODE=1
export COMPILE_NTASK=4
export COMPILE_TIME='00:10:00'

# ..................................................................
## Kernel bash options
export KERNEL_TIME='00:30:00'

# ..................................................................
export CLUSTER='vega'

# ..................................................................
#---------------------------------------------------------------------------------------
#               MELUXINA
#---------------------------------------------------------------------------------------
if [[ ${CLUSTER} == "meluxina" ]]; then
    export BASE_SCRATCHDIR=/project/scratch/p200053
    ## .................................
    ## NEEDED dir
    export HOME_FLUPS=${HOME}/flups/
    export HOME_H3LPR=${HOME}/h3lpr/

    export FFTW_DIR=${EBROOTFFTW}
    export HDF5_DIR=${EBROOTHDF5}

    ## .................................
    ## MPI information
    export MPI_VERSION=4.1.3
    export MPICC='mpicc'
    export MPICXX='mpic++'
    export NPROC_NODES=128

    ## BASH OPTIONS -- GENERAL
    export PARTITION='cpu'
    export ACCOUNT='p200053'

    ## BASH OPTIONS -- Compilation job 
    export COMPILE_CLUSTER_SPEC='--qos=short'

    ## .................................
    ## BASH OPTIONS -- kernel job 
    export KERNEL_CLUSTER_SPEC='--qos=default'

    ## .................................
    export LCOMMAND='srun'
fi


#---------------------------------------------------------------------------------------
#               VEGA
#---------------------------------------------------------------------------------------
if [[ ${CLUSTER} == "vega" ]]; then
    export BASE_SCRATCHDIR=/ceph/hpc/data/b2203-024-users/
    ## .................................
    ## NEEDED dir
    export HOME_FLUPS=${HOME}/flups/
    export HOME_H3LPR=${HOME}/h3lpr/

    export FFTW_DIR=${EBROOTFFTW}
    export HDF5_DIR=${EBROOTHDF5}

    ## .................................
    ## MPI information
    export MPI_VERSION=4.1.3
    export MPICC='mpicc'
    export MPICXX='mpic++'
    export NPROC_NODES=128

    ## BASH OPTIONS -- GENERAL
    export PARTITION='cpu'
    export ACCOUNT='b2203-024-users'

    ## BASH OPTIONS -- Compilation job 
    export COMPILE_CLUSTER_SPEC=''

    ## .................................
    ## BASH OPTIONS -- kernel job 
    export KERNEL_CLUSTER_SPEC=''

    ## .................................
    export LCOMMAND='mpirun'
fi


# ("${HOME_FLUPS}/samples/validation/run/benchmark_workflow_weak.sh")
("${HOME_FLUPS}/samples/validation/run/benchmark_workflow_strong.sh")