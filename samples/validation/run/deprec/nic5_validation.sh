#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=/home/ucl/tfl/pbalty/flups/samples/validation/
EXEC_FLUPS=flups_validation_nb

######### WEAK -> keep the size/proc constant
SCRATCH=$GLOBALSCRATCH/flups_validation/vector

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/nic5_kernel_valid.sh $SCRATCH

cd $SCRATCH

export SIZE_PER_PROC=64

############################################
## 64 PROC -> 4*4*4 = 64
export MY_NX=4
export MY_NY=4
export MY_NZ=4
export MY_SIZEX=32
export MY_SIZEY=32
export MY_SIZEZ=32
export MY_NTHREADS=1
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
export MY_NRES=5
export MY_NSOLVE=1
#export MY_KERNEL=0


for MYKERNEL in 0 1 2 3 4 5 6
do
        export MY_KERNEL=$((MYKERNEL))
        #sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
	export MY_BC="0 0 1 0 3 3"
        sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./nic5_kernel_valid.sh
        
        export MY_BC="4 4 4 4 4 4"
        sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./nic5_kernel_valid.sh
        
        export MY_BC="4 0 4 4 1 4"
        sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./nic5_kernel_valid.sh
done

#export MY_BC="0 0 1 0 3 3"
#sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./nic5_kernel_valid.sh
#
#export MY_BC="4 4 4 4 4 4"
#sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./nic5_kernel_valid.sh
#
#export MY_BC="4 0 4 4 1 4"
#sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./nic5_kernel_valid.sh
