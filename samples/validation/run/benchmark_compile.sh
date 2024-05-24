#!/bin/bash -l
#SBATCH --job-name="Flups Compilation"
#SBATCH --output=flups_compile%j.out
#SBATCH --error=flups_compile%j.err

echo "sourcing the modules -- ${SCRIPT_MODULE} ${MPI_VERSION}"
source ${SCRIPT_MODULE} ${MPI_VERSION}

## Compile H3LPR
echo "Compiling h3lpr"
cd ${H3LPR_DIR}
touch make_arch/make.${CLUSTER}
echo " ## This make_arch has been created by the run/benchmark_scripts architecture" > make_arch/make.${CLUSTER} 
echo "CXX=${MPICXX} " >> make_arch/make.${CLUSTER} 
echo "CC=${MPICC} " >> make_arch/make.${CLUSTER}
echo "CXXFLAGS :=  ${H3LPR_CXXFLAGS}" >> make_arch/make.${CLUSTER}
echo "LDFLAGS  := ${H3LPR_LDFLAGS}" >> make_arch/make.${CLUSTER}
ARCH_FILE=make_arch/make.${CLUSTER} make clean 
PREFIX=. ARCH_FILE=make_arch/make.${CLUSTER} make install -j4 
cd - 


## Compile Flups 
echo "Compiling flups -- with options OPTS=${FLUPS_OPTS}"
## Warning -- Only the static librairies are installed as we move the executable of place.. 
cd ${FLUPS_DIR} 
echo " ## This make_arch has been created by the run/benchmark_scripts architecture" > make_arch/make.${CLUSTER} 
echo "CXX=${MPICXX}" >> make_arch/make.${CLUSTER}
echo "CC=${MPICC} " >> make_arch/make.${CLUSTER}
echo "AR=${FLUPS_AR} " >> make_arch/make.${CLUSTER}

echo "CXXFLAGS := ${FLUPS_CXXFLAGS} " >> make_arch/make.${CLUSTER}
echo "CCFLAGS  := ${FLUPS_CCFLAGS}" >> make_arch/make.${CLUSTER}
echo "LDFLAGS  := ${FLUPS_LDFLAGS}"  >> make_arch/make.${CLUSTER}
echo "OPTS  := ${FLUPS_OPTS}"  >> make_arch/make.${CLUSTER}

echo "#---------------------------------------------------------   " >> make_arch/make.${CLUSTER}
echo "# DEPENDENCES DIRECTORIES  " >> make_arch/make.${CLUSTER}
echo "#---------------------------------------------------------  " >> make_arch/make.${CLUSTER}
echo "  " >> make_arch/make.${CLUSTER}
echo "## FFTW3  " >> make_arch/make.${CLUSTER}
echo "FFTW_DIR  := ${FFTW_DIR}  " >> make_arch/make.${CLUSTER}
echo "FFTW_LIB := ${FFTW_DIR}/lib  " >> make_arch/make.${CLUSTER}
echo "FFTW_INC := ${FFTW_DIR}/include  " >> make_arch/make.${CLUSTER}
echo "  " >> make_arch/make.${CLUSTER}
echo "## HDF5  " >> make_arch/make.${CLUSTER}
echo "HDF5_DIR := ${HDF5_DIR}  " >> make_arch/make.${CLUSTER}
echo "HDF5_LIB := ${HDF5_DIR}/lib  " >> make_arch/make.${CLUSTER}
echo "HDF5_INC := ${HDF5_DIR}/include  " >> make_arch/make.${CLUSTER}
echo "  " >> make_arch/make.${CLUSTER}
echo "##H3LPR  " >> make_arch/make.${CLUSTER}
echo "H3LPR_INC := ${H3LPR_DIR}/include/  " >> make_arch/make.${CLUSTER}
echo "H3LPR_LIB := ${H3LPR_DIR}/lib/  " >> make_arch/make.${CLUSTER}

ARCH_FILE=make_arch/make.${CLUSTER} make clean
PREFIX=. ARCH_FILE=make_arch/make.${CLUSTER} make install_static -j4
cd - 


## Compile the validation testcase
echo "Compiling the validation client"
cd ${FLUPS_DIR}/samples/validation
ARCH_FILE=make_arch/make.${CLUSTER} make clean
ARCH_FILE=make_arch/make.${CLUSTER} make all -j4
cd -

