#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=/home/acad/ucl-tfl/tgillis/flups/samples/validation
EXEC_FLUPS=flups_validation_nb

######### STRONG -> increase the number of CPU, same size
SCRATCH=/SCRATCH/acad/examples/tgillis/flups_strong_nb

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/zenobe_kernel_nb.sh $SCRATCH

cd $SCRATCH

#####################################################################
export MY_RES=1024
#####################################################################

## 128 PROCS
qsub -q large -v MY_NX=8,MY_NY=4,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=32:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=4,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=2,L_X=1,L_Y=1,L_Z=1 -l select=32:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=4,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=4,L_X=1,L_Y=1,L_Z=1 -l select=32:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

## 256 PROCS
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=64:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=2,L_X=1,L_Y=1,L_Z=1 -l select=64:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=4,L_X=1,L_Y=1,L_Z=1 -l select=64:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

## 512 PROCS
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=8,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=128:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=2,L_X=1,L_Y=1,L_Z=1 -l select=128:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=4,L_X=1,L_Y=1,L_Z=1 -l select=128:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

## 1024 PROCS
qsub -q large -v MY_NX=16,MY_NY=8,MY_NZ=8,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=256:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=8,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=2,L_X=1,L_Y=1,L_Z=1 -l select=256:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=8,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=4,L_X=1,L_Y=1,L_Z=1 -l select=256:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

## 2048 PROCS
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=8,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=512:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=2,L_X=1,L_Y=1,L_Z=1 -l select=512:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=4,L_X=1,L_Y=1,L_Z=1 -l select=512:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

## 4096 PROCS
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=16,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=8,MY_SIZE=${MY_RES},MY_NTH=2,L_X=1,L_Y=1,L_Z=1 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=4,L_X=1,L_Y=1,L_Z=1 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

#end of file
