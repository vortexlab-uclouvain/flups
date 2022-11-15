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

export FLUPS_AR='gcc-ar'
export FLUPS_CCFLAGS='-g -O3 -std=c99   -DMPI_BATCH_SEND=1 -march=native -flto'
export FLUPS_CXXFLAGS='-g -O3 -std=c++17 -DMPI_BATCH_SEND=1  -march=native -flto'
export FLUPS_LDFLAGS='-fopenmp -flto'
export FLUPS_OPTS='-DNDEBUG -DFFTW_FLAG=FFTW_MEASURE'

# ..................................................................
## Compile bash options
export COMPILE_NNODE=1
export COMPILE_NTASK=4
export COMPILE_TIME='05:00:00'

# ..................................................................
## Kernel bash options
export KERNEL_TIME='01:00:00'

# ..................................................................
#export CLUSTER='meluxina'
export CLUSTER='vega'
# ..................................................................
#---------------------------------------------------------------------------------------
#               MELUXINA
#---------------------------------------------------------------------------------------
if [[ ${CLUSTER} == "meluxina" ]]; then
    export BASE_SCRATCHDIR=/project/scratch/p200053
    # export BASE_SCRATCHDIR=/project/scratch/p200067
    
    ## .................................
    ## NEEDED dir
    export HOME_FLUPS=${HOME}/flups/
    export HOME_H3LPR=${HOME}/h3lpr/

    #export FFTW_DIR=${EBROOTFFTW}
    #export HDF5_DIR=${EBROOTHDF5}
    export FFTW_DIR=${HOME}/lib-MPICH-4.1a1-UCX-1.13.1-fast
    export HDF5_DIR=${HOME}/lib-MPICH-4.1a1-UCX-1.13.1-fast

    ## .................................
    ## MPI information
    #export MPI_VERSION=4.1.4
    export MPI_VERSION=4.1a1
    export MPICC='mpicc'
    export MPICXX='mpic++'
    export NPROC_NODES=128
    
    #export ACCFFT_DIR=${HOME}/lib-OpenMPI-${MPI_VERSION}
    export ACCFFT_DIR=${HOME}/lib-MPICH-4.1a1-UCX-1.13.1-fast

    ## .................................

    ## BASH OPTIONS -- GENERAL
    export PARTITION='cpu'
    export ACCOUNT='p200053'
    #export ACCOUNT='p200067'

    ## BASH OPTIONS -- Compilation job 
    export COMPILE_CLUSTER_SPEC='--qos=short --mem=491520'
#    export COMPILE_CLUSTER_SPEC='--res verylargecpu --qos=large --mem=491520'
    ## .................................
    ## BASH OPTIONS -- kernel job 
    export KERNEL_CLUSTER_SPEC='--qos=default --mem=491520'
    #export KERNEL_CLUSTER_SPEC='--res verylargecpu --qos=large --mem=491520'
    
    export UCX_TLS=self,shm,dc
    export UCX_DC_MLX5_NUM_DCI=16
    # export UCX_TLS=self,shm,rc,ud
    # export UCX_TLS=self,shm,rc,ud
    # export UCX_UD_TX_QUEUE_LEN=4096

    #export LAUNCH_COMMAND='mpirun'
    export LAUNCH_COMMAND='mpiexec -bind-to core'
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

    #export FFTW_DIR=${EBROOTFFTW}
    #export HDF5_DIR=${EBROOTHDF5}
    export FFTW_DIR=${HOME}/lib-MPICH-4.1a1
    export HDF5_DIR=${HOME}/lib-MPICH-4.1a1

    #export ACCFFT_DIR=${HOME}/lib-OpenMPI-${MPI_VERSION}
    export ACCFFT_DIR=${HOME}/lib-MPI-4.1a1/

    ## .................................
    ## MPI information
    #export MPI_VERSION=4.1.3
    export MPI_VERSION=4.1a1
    export MPICC='mpicc'
    export MPICXX='mpic++'
    
    ## .................................
    ## Cluster information
    export NPROC_NODES=128

    ## BASH OPTIONS -- GENERAL
    export PARTITION='cpu'
    export ACCOUNT='b2203-024-users'

    ## BASH OPTIONS -- Compilation job 
    export COMPILE_CLUSTER_SPEC='--mem=256000 --reservation=Benchmark-2301733'

    ## .................................
    ## BASH OPTIONS -- kernel job 
    export KERNEL_CLUSTER_SPEC='--mem=256000 --reservation=Benchmark-2301733'

    ## .................................
    # export UCX_TLS=ud,sm
    export UCX_TLS=self,shm,dc
    #export OMPI_MCA_pml=ucx
    #export OMPI_MCA_osc=ucx

    export LAUNCH_COMMAND='mpiexec -bind-to core'
fi


("${HOME_FLUPS}/samples/compareACCFFT/run/benchmark_workflow_compare.sh")
