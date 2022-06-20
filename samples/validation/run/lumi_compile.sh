#!/bin/bash -l
#SBATCH --job-name="Flups Compilation"
#SBATCH --account=project_465000098
#SBATCH --partition=small
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --output=flups_compile%j.out
#SBATCH --error=flups_compile%j.err


source ${MODULES} ${CMPICHVERSION}

## Compile H3LPR
cd ${H3LPR_DIR}
touch make_arch/make.lumi
echo "CXX=CC " >> make_arch/make.lumi 
echo "CC=cc " > make_arch/make.lumi
echo "CXXFLAGS := -g -O3 -march=native -fopenmp -DNDEBUG " >> make_arch/make.lumi
echo "LDFLAGS  := -fopenmp " >> make_arch/make.lumi
ARCH_FILE=make_arch/make.lumi make install -j4 
cd - 


## Compile Flups 
## Warning -- Only the static librairies are installed as we move the executable of place.. 
cd ${FLUPS_DIR} 
echo "CXX=CC" > make_arch/make.lumi
echo "CC=cc " >> make_arch/make.lumi
echo "CXXFLAGS := -g -O3 -std=c++17 -DNDEBUG " >> make_arch/make.lumi
echo "CCFLAGS  := -g -O3 -std=c99 -DNDEBUG" >> make_arch/make.lumi
echo "LDFLAGS  := -fopenmp #-lstc++ "  >> make_arch/make.lumi

echo "#---------------------------------------------------------   " >> make_arch/make.lumi
echo "# DEPENDENCES DIRECTORIES  " >> make_arch/make.lumi
echo "#---------------------------------------------------------  " >> make_arch/make.lumi
echo "  " >> make_arch/make.lumi
echo "## FFTW3  " >> make_arch/make.lumi
echo "FFTW_DIR  := ${FFTW_DIR}  " >> make_arch/make.lumi
echo "FFTW_LIB := ${FFTW_DIR}/lib  " >> make_arch/make.lumi
echo "FFTW_INC := ${FFTW_DIR}/include  " >> make_arch/make.lumi
echo "  " >> make_arch/make.lumi
echo "## HDF5  " >> make_arch/make.lumi
echo "HDF5_DIR := ${HDF5_DIR}  " >> make_arch/make.lumi
echo "HDF5_LIB := ${HDF5_DIR}/lib  " >> make_arch/make.lumi
echo "HDF5_INC := ${HDF5_DIR}/include  " >> make_arch/make.lumi
echo "  " >> make_arch/make.lumi
echo "##H3LPR  " >> make_arch/make.lumi
echo "H3LPR_INC := ${H3LPR_DIR}/include/  " >> make_arch/make.lumi
echo "H3LPR_LIB := ${H3LPR_DIR}/lib/  " >> make_arch/make.lumi

ARCH_FILE=make_arch/make.lumi make install_static -j4
cd - 


## Compile the validation testcase
cd ${FLUPS_DIR}/samples/validation
ARCH_FILE=make_arch/make.lumi make all -j4
cd -

