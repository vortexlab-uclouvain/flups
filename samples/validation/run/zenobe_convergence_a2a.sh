#!/bin/sh
## RM the previous validation dir

#HOME_FLUPS=/home/acad/ucl-tfl/dcaprace/FLUPS/flups_green/samples/validation
HOME_FLUPS=/home/acad/ucl-tfl/tgillis/flups/samples/validation
EXEC_FLUPS=flups_validation_a2a

#SCRATCH=/SCRATCH/acad/examples/dcaprace/flups_convergence_a2a
SCRATCH=/SCRATCH/acad/examples/tgillis/flups_convergence_a2a

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/zenobe_kernelConv_a2a.sh $SCRATCH

cd $SCRATCH

#####################################################################

## 64 
export MY_RES=64
#qsub -q main -v MY_NX=1,MY_NY=2,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a.sh
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a_3d.sh

## 128
export MY_RES=128
#qsub -q main -v MY_NX=1,MY_NY=2,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a.sh
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a_3d.sh

## 256 
export MY_RES=256
#qsub -q main -v MY_NX=1,MY_NY=2,MY_NZ=2,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a.sh
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=1:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a_3d.sh

## 512 
export MY_RES=512
#qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=4,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=4:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a.sh
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=4:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernelConv_a2a_3d.sh

## 1024 
export MY_RES=1024
#qsub -q large -v MY_NX=4,MY_NY=6,MY_NZ=8,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=8:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelConv_a2a.sh
qsub -q large -v MY_NX=24,MY_NY=16,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=16:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelConv_a2a_3d.sh

## 2048 
export MY_RES=2048
#qsub -q large -v MY_NX=8,MY_NY=12,MY_NZ=16,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=64:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelConv_a2a.sh
qsub -q large -v MY_NX=48,MY_NY=48,MY_NZ=1,MY_SIZE=${MY_RES},MY_NTH=1,L_X=1,L_Y=1,L_Z=1 -l select=96:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelConv_a2a_3d.sh

#end of file
