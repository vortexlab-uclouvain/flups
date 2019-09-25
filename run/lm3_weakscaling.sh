#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=/home/ucl/tfl/tgillis/flups
EXEC_FLUPS=flups_validation

######### WEAK -> keep the size/proc constant
SCRATCH=$GLOBALSCRATCH/flups_weak

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/lm3_kernel_valid.sh $SCRATCH

cd $SCRATCH

export SIZE_PER_PROC=64

############################################
## 1 PROC -> 1*64 = 64
export MY_NX=1
export MY_NY=1
export MY_SIZE=$(($SIZE_PER_PROC * MY_NX))
#-- 1 thread
export MY_NZ=1
export MY_NTHREADS=1
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh

############################################
## 2 PROC -> 2*64 = 128
export MY_NX=2
export MY_NY=2
export MY_SIZE=$(($SIZE_PER_PROC * MY_NX))
#-- 1 thread
export MY_NZ=2
export MY_NTHREADS=1
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
#-- 2 thread
export MY_NZ=1
export MY_NTHREADS=2
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh


############################################
## 4 PROC -> 4*64 = 256
export MY_NX=4
export MY_NY=4
export MY_SIZE=$(($SIZE_PER_PROC * MY_NX))
#-- 1 thread
export MY_NZ=4
export MY_NTHREADS=1
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
#-- 2 thread
export MY_NZ=2
export MY_NTHREADS=2
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh


############################################
## 8 PROC -> 8*64 = 512
export MY_NX=8
export MY_NY=8
export MY_SIZE=$(($SIZE_PER_PROC * MY_NX))
#-- 1 thread
export MY_NZ=8
export MY_NTHREADS=1
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
#-- 2 thread
export MY_NZ=4
export MY_NTHREADS=2
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh


