#!/bin/bash -l
#SBATCH --job-name="Flups Compilation"
#SBATCH --output=flups_compile%j.out
#SBATCH --error=flups_compile%j.err

echo "sourcing the modules"
source ${SCRIPT_MODULE} ${MPI_VERSION}

## Compile H3LPR
cd ${H3LPR_DIR}
touch make_arch/make.${CLUSTER}
echo "CXX=${MPICXX} " > make_arch/make.${CLUSTER} 
echo "CC=${MPICC} " >> make_arch/make.${CLUSTER}
echo "CXXFLAGS :=  ${H3LPR_CXXFLAGS}" >> make_arch/make.${CLUSTER}
echo "LDFLAGS  := ${H3LPR_LDFLAGS}" >> make_arch/make.${CLUSTER}
ARCH_FILE=make_arch/make.${CLUSTER} make clean
PREFIX=. ARCH_FILE=make_arch/make.${CLUSTER} make install -j4
cd - 


## Compile Flups 
## Warning -- Only the static librairies are installed as we move the executable of place.. 
cd ${FLUPS_DIR} 
echo "CXX=${MPICXX}" > make_arch/make.${CLUSTER}
echo "CC=${MPICC} " >> make_arch/make.${CLUSTER}

echo "CXXFLAGS := ${FLUPS_CXXFLAGS} " >> make_arch/make.${CLUSTER}
echo "CCFLAGS  := ${FLUPS_CCFLAGS}" >> make_arch/make.${CLUSTER}

echo "LDFLAGS  := ${FLUPS_LDFLAGS}"  >> make_arch/make.${CLUSTER}
echo "AR := ${FLUPS_AR} " >> make_arch/make.${CLUSTER}

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
echo "  " >> make_arch/make.${CLUSTER}
echo "##ACCFFT  " >> make_arch/make.${CLUSTER}
echo "ACCFFT_DIR := ${ACCFFT_DIR}" >> make_arch/make.${CLUSTER}

ARCH_FILE=make_arch/make.${CLUSTER} make clean
PREFIX=. ARCH_FILE=make_arch/make.${CLUSTER} make install -j4

cd - 


## Compile the validation testcase
cd ${FLUPS_DIR}/samples/compareACCFFT
ARCH_FILE=make_arch/make.${CLUSTER} make clean
ARCH_FILE=make_arch/make.${CLUSTER} make all -j4
cd -

