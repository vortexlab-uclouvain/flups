
## Change the list of module we use 
#module use /apps/USE/easybuild/staging/2021.5/modules/all
#
### First load the librairies relying on the old version of OpenMPI
#module load HDF5/1.12.1-gompi-2021a
#module load FFTW/3.3.10-gompi-2021a 
#
### Then load the latest version of OpenMPI
##module load OpenMPI/5.0.0-GCC-10.3.0
#module load OpenMPI/${1}-GCC-10.3.0
##module load OpenMPI/4.1.1-GCC-10.3.0

module use /apps/USE/easybuild/release/2021.5/modules/all
module load HDF5/1.12.1-gompi-2021a
module load FFTW/3.3.10-gompi-2021a
module load OpenMPI/4.1.4-GCC-10.3.0
module load CMake/3.20.1-GCCcore-10.3.0
