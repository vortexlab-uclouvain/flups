#!/bin/bash -l
#SBATCH --job-name="Flups Compilation"
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=%x.out
#SBATCH --error=%x.err

## Load modules 
# module use /apps/USE/easybuild/staging/2021.1/modules/all
## First load the librairies relying on the old version of OpenMPI
# module load HDF5/1.12.1-gompi-2021a
# module load FFTW/3.3.10-gompi-2021a 
## Then load the latest version of OpenMPI
#module load OpenMPI/5.0.0-GCC-10.3.0
# module load OpenMPI/4.1.3-GCC-10.3.0 
#module load OpenMPI/4.1.1-GCC-10.3.0

source ${MODULES} ${OMPIVERSION}

## Compile H3LPR
cd ${H3LPR_DIR}
touch make_arch/make.meluxina
echo "CC = mpicc " > make_arch/make.meluxina
echo "CXX = mpic++ " >> make_arch/make.meluxina 
echo "CXXFLAGS := -g -O3 -march=native -fopenmp -DNDEBUG " >> make_arch/make.meluxina
echo "LDFLAGS := -fopenmp " >> make_arch/make.meluxina
ARCH_FILE=make_arch/make.meluxina make install -j 
cd - 


## Compile Flups 
## Warning -- Only the static librairies are installed as we move the executable of place.. 
cd ${FLUPS_DIR} 
echo "CXX = ${EBROOTOPENMPI}/bin/mpic++ " > make_arch/make.meluxina
echo "CC = ${EBROOTOPENMPI}/bin/mpicc " >> make_arch/make.meluxina
echo "## For production: " >> make_arch/make.meluxina
echo "CXXFLAGS := -O3 -std=c++17 -DNDEBUG -Wno-variadic-macros " >> make_arch/make.meluxina
echo "CCFLAGS := -O3 -std=c99 -DNDEBUG -Wno-variadic-macros " >> make_arch/make.meluxina
# echo "CXXFLAGS := -O0 -std=c++17 " >> make_arch/make.meluxina
# echo "CCFLAGS := -O0 -std=c99 " >> make_arch/make.meluxina
echo "LDFLAGS = -fopenmp #-lstc++ " >> make_arch/make.meluxina
echo "#---------------------------------------------------------   " >> make_arch/make.meluxina
echo "# DEPENDENCES DIRECTORIES  " >> make_arch/make.meluxina
echo "#---------------------------------------------------------  " >> make_arch/make.meluxina
echo "  " >> make_arch/make.meluxina
echo "## FFTW3  " >> make_arch/make.meluxina
echo "FFTW_DIR  := ${EBROOTFFTW}  " >> make_arch/make.meluxina
echo "FFTW_LIB := ${EBROOTFFTW}/lib  " >> make_arch/make.meluxina
echo "FFTW_INC := ${EBROOTFFTW}/include  " >> make_arch/make.meluxina
echo "  " >> make_arch/make.meluxina
echo "## HDF5  " >> make_arch/make.meluxina
echo "HDF5_DIR  := ${EBROOTHDF5}  " >> make_arch/make.meluxina
echo "HDF5_LIB := ${EBROOTHDF5}/lib  " >> make_arch/make.meluxina
echo "HDF5_INC := ${EBROOTHDF5}/include  " >> make_arch/make.meluxina
echo "  " >> make_arch/make.meluxina
echo "##H3LPR  " >> make_arch/make.meluxina
echo "H3LPR_INC := ${H3LPR_DIR}/include/  " >> make_arch/make.meluxina
echo "H3LPR_LIB := ${H3LPR_DIR}/lib/  " >> make_arch/make.meluxina

ARCH_FILE=make_arch/make.meluxina make install_static -j
cd - 


## Compile the validation testcase
cd ${FLUPS_DIR}/samples/validation
ARCH_FILE=make_arch/make.meluxina make all -j 
cd -
