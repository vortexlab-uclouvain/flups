#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=/home/ucl/tfl/tgillis/flups
EXEC_FLUPS=flups_validation

######### STRONG -> increase the number of CPU, same size
SCRATCH=$GLOBALSCRATCH/flups_strong

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/lm3_kernel_valid.sh $SCRATCH

cd $SCRATCH

export MY_SIZE=512

### 8 PROCS (max 24)
#export MY_NX=2
#export MY_NY=2
#export MY_NZ=2
#export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
#
#for MYTHREADS in 1 2 3
#do		
#	export MY_NTHREADS=$((MYTHREADS))
#	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
#done
#
#
### 16 PROCS (max 48)
#export MY_NX=4
#export MY_NY=2
#export MY_NZ=2
#export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))
#
#for MYTHREADS in 1 2 3
#do		
#	export MY_NTHREADS=$((MYTHREADS))
#	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
#done
#
## 32 PROCS (max 96)
export MY_NX=4
export MY_NY=4
export MY_NZ=2
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))

for MYTHREADS in 1 2 3
do		
	export MY_NTHREADS=$((MYTHREADS))
	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
done

## 64 PROCS (max 192)
export MY_NX=4
export MY_NY=4
export MY_NZ=4
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))

for MYTHREADS in 1 2 3
do		
	export MY_NTHREADS=$((MYTHREADS))
	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
done

## 128 PROCS (max 384)
export MY_NX=8
export MY_NY=4
export MY_NZ=4
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))

for MYTHREADS in 1 2 3
do		
	export MY_NTHREADS=$((MYTHREADS))
	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
done

## 256 PROCS (max 768)
export MY_NX=8
export MY_NY=8
export MY_NZ=4
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))

for MYTHREADS in 1 2 3
do		
	export MY_NTHREADS=$((MYTHREADS))
	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
done


## 512 PROCS (max 1536)
export MY_NX=8
export MY_NY=8
export MY_NZ=8
export MY_NTASKS=$(($MY_NX*$MY_NY*$MY_NZ))

for MYTHREADS in 1 2 3
do		
	export MY_NTHREADS=$((MYTHREADS))
	sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} ./lm3_kernel_valid.sh
done

