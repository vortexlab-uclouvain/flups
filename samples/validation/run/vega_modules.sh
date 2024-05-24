#module purge
#module load HDF5/1.10.7-gompi-2021a
#module load FFTW/3.3.9-gompi-2021a
#module load OpenMPI/${1}-GCC-11.2.0

#module load GCC/10.3.0
#module load Automake/1.16.3-GCCcore-10.3.0
#module load Autoconf/2.71-GCCcore-10.3.0
#module load libtool/2.4.6-GCCcore-10.3.0
#module load CMake/3.20.1-GCCcore-10.3.0

export PATH=${HOME}/lib-MPICH-4.1a1-UCX-1.13.1-fast/bin/:${PATH}
export LD_LIBRARY_PATH=${HOME}/lib-MPICH-4.1a1-UCX-1.13.1-fast/lib:${LD_LIBRARY_PATH}
