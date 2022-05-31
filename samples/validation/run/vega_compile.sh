#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --account=d2202-040-users
#SBATCH --job-name="flups Compilation"
#SBATCH --ntasks=4
#SBATCH --time=00:10:00

source ${MODULES} ${OMPIVERSION}

## Compile H3LPR
cd ${H3LPR_DIR}
touch make_arch/make.vega
echo "CC = mpicc " > make_arch/make.vega
echo "CXX = mpic++ " >> make_arch/make.vega 
echo "CXXFLAGS := -g -O3 -march=native -fopenmp -DNDEBUG " >> make_arch/make.vega
echo "LDFLAGS := -fopenmp " >> make_arch/make.vega
PREFIX=. ARCH_FILE=make_arch/make.vega make install -j4
cd - 


## Compile Flups 
## Warning -- Only the static librairies are installed as we move the executable of place.. 
cd ${FLUPS_DIR} 
echo "CXX = mpic++ " > make_arch/make.vega
echo "CC = mpicc " >> make_arch/make.vega
echo "## For production: " >> make_arch/make.vega
echo "CXXFLAGS := -O3 -std=c++17 -DNDEBUG -march=native -DHAVE_WISDOM=\\\"${HOME}/flups/fftw-wisdom/wisdom/vega.wsdm\\\" ${COMPILE_OPT}" >> make_arch/make.vega
echo "CCFLAGS := -O3 -std=c99 -DNDEBUG " >> make_arch/make.vega
echo "LDFLAGS := -fopenmp" >> make_arch/make.vega
#echo "CXXFLAGS := -O3 -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -std=c++17 " >> make_arch/make.vega
#echo "CCFLAGS :=  -O3 -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -std=c99 " >> make_arch/make.vega
#echo "LDFLAGS =       -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -fopenmp #-lstc++ " >> make_arch/make.vega
echo "#---------------------------------------------------------   " >> make_arch/make.vega
echo "# DEPENDENCES DIRECTORIES  " >> make_arch/make.vega
echo "#---------------------------------------------------------  " >> make_arch/make.vega
echo "  " >> make_arch/make.vega
echo "## FFTW3  " >> make_arch/make.vega
#echo "FFTW_DIR  := ${EBROOTFFTW}  " >> make_arch/make.vega
echo "FFTW_LIB := ${FFTW_DIR}/lib  " >> make_arch/make.vega
echo "FFTW_INC := ${FFTW_DIR}/include  " >> make_arch/make.vega
echo "  " >> make_arch/make.vega
echo "## HDF5  " >> make_arch/make.vega
#echo "HDF5_DIR  := ${EBROOTHDF5}  " >> make_arch/make.vega
echo "HDF5_LIB := ${HDF5_DIR}/lib  " >> make_arch/make.vega
echo "HDF5_INC := ${HDF5_DIR}/include  " >> make_arch/make.vega
echo "  " >> make_arch/make.vega
echo "##H3LPR  " >> make_arch/make.vega
echo "H3LPR_INC := ${H3LPR_DIR}/include/  " >> make_arch/make.vega
echo "H3LPR_LIB := ${H3LPR_DIR}/lib/  " >> make_arch/make.vega

PREFIX=. ARCH_FILE=make_arch/make.vega make clean
PREFIX=. ARCH_FILE=make_arch/make.vega make install_static -j4
cd - 

echo "FLUPS is done, doing the client now"
## Compile the validation testcase
cd ${FLUPS_DIR}/samples/validation
ARCH_FILE=make_arch/make.vega make clean
ARCH_FILE=make_arch/make.vega make all -j4

#if [-z ${COMPILE_SUFFIX+x}];
#then mv flups_validation_a2a flups_validation_a2a_${COMPILE_SUFFIX} && mv flups_validation_nb flups_validation_nb_${COMPILE_SUFFIX};
#fi

cd -
