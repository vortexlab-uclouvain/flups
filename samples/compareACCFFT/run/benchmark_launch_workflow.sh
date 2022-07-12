#!/bin/bash -l
#---------------------------------------------------------------------------------------
#               Case dependent info
#---------------------------------------------------------------------------------------
## Version: 
# All to all                 = a2a 
# Non blocking               = nb 
# MPI_datatypes              = isr
# First version all to all   = dprec_a2a 
# First version non blocking = dprec_nb
export CODE_VERSION='a2a nb isr dprec_nb dprec_a2a'

## Compile the librairies  
export H3LPR_CXXFLAGS='-g -O3 -march=native -fopenmp -DNDEBUG -flto'
export H3LPR_LDFLAGS='-fopenmp -rdynamic -ldl -flto'

export FLUPS_CCFLAGS='-g -O3 -std=c99   -DMPI_BATCH_SEND=1 -DNDEBUG -march=native -flto'
export FLUPS_CXXFLAGS='-g -O3 -std=c++17 -DMPI_BATCH_SEND=1 -DNDEBUG -march=native -flto'
export FLUPS_LDFLAGS='-fopenmp -flto'
export FLUPS_OPTS=''

# ..................................................................
## Compile bash options
export COMPILE_NNODE=1
export COMPILE_NTASK=4
export COMPILE_TIME='05:00:00'

# ..................................................................
## Kernel bash options
export KERNEL_TIME='01:00:00'

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

    # export UCX_TLS=self,shm,rc,ud
    # export UCX_TLS=self,shm,rc,ud
    # export UCX_UD_TX_QUEUE_LEN=4096

fi


#---------------------------------------------------------------------------------------
#               VEGA
#---------------------------------------------------------------------------------------
if [[ ${CLUSTER} == "vega" ]]; then
    export BASE_SCRATCHDIR=/ceph/hpc/data/b2203-024-users/
    ## .................................
    ## MPI information
    export MPI_VERSION=4.1.3
    export MPICC='mpicc'
    export MPICXX='mpic++'
    
    ## .................................
    ## NEEDED dir
    export HOME_FLUPS=${HOME}/flups/
    export HOME_H3LPR=${HOME}/h3lpr/

    export FFTW_DIR=${EBROOTFFTW}
    export HDF5_DIR=${EBROOTHDF5}

    export ACCFFT_DIR=${HOME}/lib-OpenMPI-${MPI_VERSION}

    ## .................................
    ## Cluster information
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
    # export UCX_TLS=ud,sm
    export UCX_TLS=self,shm,rc,dc
    export OMPI_MCA_pml=ucx
    export OMPI_MCA_osc=ucx
fi


# Performance 
# ("${HOME_FLUPS}/samples/validation/run/benchmark_workflow_weak.sh")
# ("${HOME_FLUPS}/samples/validation/run/benchmark_workflow_strong.sh")

# Validation 
("${HOME_FLUPS}/samples/compareACCFFT/run/benchmark_workflow_compare.sh")
