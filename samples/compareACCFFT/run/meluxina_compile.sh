#!/bin/bash -l
#SBATCH --job-name="Flups Compilation"
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00

source ${MODULES} ${OMPIVERSION}

## Compile H3LPR
cd ${H3LPR_DIR}
touch make_arch/make.meluxina
echo "CC = mpicc " > make_arch/make.meluxina
echo "CXX = mpic++ " >> make_arch/make.meluxina 
echo "CXXFLAGS := -g -O3 -march=native -fopenmp -DNDEBUG " >> make_arch/make.meluxina
echo "LDFLAGS := -fopenmp -ldl " >> make_arch/make.meluxina
ARCH_FILE=make_arch/make.meluxina make install -j 
cd - 


## Compile Flups 
## Warning -- Only the static librairies are installed as we move the executable of place.. 
cd ${FLUPS_DIR} 
echo "CXX = ${EBROOTOPENMPI}/bin/mpic++ " > make_arch/make.meluxina
echo "CC = ${EBROOTOPENMPI}/bin/mpicc " >> make_arch/make.meluxina
#echo "CXXFLAGS := -g -O3 -std=c++17 -DMPI_BATCH_SEND=2147483647 -DBALANCE_DPREC -DFFTW_FLAG=FFTW_MEASURE -DNDEBUG -march=native #-DNO_PROF" >> make_arch/make.meluxina
#echo "CCFLAGS :=  -g -O3 -std=c99   -DMPI_BATCH_SEND=2147483647 -DBALANCE_DPREC -DFFTW_FLAG=FFTW_MEASURE -DNDEBUG -march=native #-DNO_PROF" >> make_arch/make.meluxina
#echo "CXXFLAGS := -g -O3 -std=c++17 -DMPI_BATCH_SEND=2147483647 -DFFTW_FLAG=FFTW_MEASURE -DNDEBUG -march=native #-DNO_PROF" >> make_arch/make.meluxina
#echo "CCFLAGS :=  -g -O3 -std=c99   -DMPI_BATCH_SEND=2147483647 -DFFTW_FLAG=FFTW_MEASURE -DNDEBUG -march=native #-DNO_PROF" >> make_arch/make.meluxina
echo "CXXFLAGS := -g -O3 -std=c++17 -DMPI_BATCH_SEND=2147483647 -DNDEBUG -march=native #-DNO_PROF" >> make_arch/make.meluxina
echo "CCFLAGS :=  -g -O3 -std=c99   -DMPI_BATCH_SEND=2147483647 -DNDEBUG -march=native #-DNO_PROF" >> make_arch/make.meluxina

# Working configuration
# echo "CXXFLAGS := -g -Wall -O0 --debug -std=c++11" >> make_arch/make.meluxina
# echo "CCFLAGS := -g -Wall -O0 --debug -std=c99" >> make_arch/make.meluxina

echo "LDFLAGS = -fopenmp  #-lstc++ " >> make_arch/make.meluxina
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
echo "  " >> make_arch/make.meluxina
echo "##ACCFFT  " >> make_arch/make.meluxina
echo "ACCFFT_DIR := ${ACCFFT_DIR}" >> make_arch/make.meluxina

ARCH_FILE=make_arch/make.meluxina make install -j
cd - 


## Compile the validation testcase
cd ${FLUPS_DIR}/samples/compareACCFFT
ARCH_FILE=make_arch/make.meluxina make all -j 
cd -
