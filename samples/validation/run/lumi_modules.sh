## First load the librairies relying on the old version of OpenMPI
#module load cray-mpich/${1}
#module load cray-hdf5-parallel/1.12.0.7
#module load cray-fftw/3.3.8.12
#module load CrayEnv
module load gcc/11.2.0
module load libfabric/1.15.0.0
export PATH=${HOME}/lib-MPICH-4.1a1/bin/:${PATH}
export LD_LIBRARY_PATH=${HOME}/lib-MPICH-4.1a1/lib:${LD_LIBRARY_PATH}
