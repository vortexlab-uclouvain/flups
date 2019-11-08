#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=/home/pr1ekp00/$(whoami)/flups/samples/validation


## fixed parameters
export SIZE_PER_PROC=128

export ver=$1

if [[ -z $1 ]];
then
  echo "you must specify a version (small/large) as an argument"
  exit 1
else
  echo "starting as $1"
fi


if [ "$ver" = "small" ]; then

############################################################
############################################################
############################################################
#                           SMALL
############################################################
############################################################
############################################################


############################################################
# ALL TO ALL
#-----------------------------------------------------------
export EXEC_FLUPS=flups_validation_a2a

SCRATCH=/gpfs/scratch/pr1ekp00/$(whoami)/flups_weak_a2a_${ver}

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof

# copy the needed info
cp $HOME_FLUPS/${EXEC_FLUPS} $SCRATCH
cp $HOME_FLUPS/run/marenostrum_kernel_valid.sh $SCRATCH
# go to it
cd $SCRATCH


#================== 1152 CPU's ================
#-- requested walltime
export WT='00:10:00'
#-- proc domain
export MY_NX=8
export MY_NY=12
export MY_NZ=12
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 4 thread
export WT='00:20:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} ")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} / 2")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#================== 2304 CPU's ================
#-- requested walltime
export WT='00:15:00'
#-- proc domain
export MY_NX=16
export MY_NY=12
export MY_NZ=12
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 4 thread
export WT='00:30:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} / 2")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} ")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh


#================== 4608 CPU's ================
#-- requested walltime
export WT='00:20:00'
#-- proc domain
export MY_NX=16
export MY_NY=24
export MY_NZ=12
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 4 thread
export WT='00:40:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} / 2")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} ")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh





############################################################
# NON-BLOCKING
#-----------------------------------------------------------
export EXEC_FLUPS=flups_validation_nb

SCRATCH=/gpfs/scratch/pr1ekp00/$(whoami)/flups_weak_nb_${ver}

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof

# copy the needed info
cp $HOME_FLUPS/${EXEC_FLUPS} $SCRATCH
cp $HOME_FLUPS/run/marenostrum_kernel_valid.sh $SCRATCH
# go to it
cd $SCRATCH


#================== 1152 CPU's ================
#-- requested walltime
export WT='00:10:00'
#-- proc domain
export MY_NX=8
export MY_NY=12
export MY_NZ=12
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 8 thread
export WT='00:20:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} ")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} / 2")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#================== 2304 CPU's ================
#-- requested walltime
export WT='00:15:00'
#-- proc domain
export MY_NX=16
export MY_NY=12
export MY_NZ=12
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 4 thread
export WT='00:30:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} / 2")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} ")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#================== 4608 CPU's ================
#-- requested walltime
export WT='00:20:00'
#-- proc domain
export MY_NX=16
export MY_NY=24
export MY_NZ=12
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 4 thread
export WT='00:40:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} / 2")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} ")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh





elif [ "$ver" = "small" ]; then

############################################################
############################################################
############################################################
#                           LARGE
############################################################
############################################################
############################################################

exit 0

#================== 9216 CPU's ================
#-- requested walltime
export WT='00:25:00'
#-- proc domain
export MY_NX=16
export MY_NY=24
export MY_NZ=24
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- requested walltime

#-- 4 thread
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} ")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} / 2")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh



############################################################
# NON-BLOCKING
#-----------------------------------------------------------



#================== 9216 CPU's ================
#-- requested walltime
export WT='00:25:00'
#-- proc domain
export MY_NX=16
export MY_NY=24
export MY_NZ=24
#-- domain length
export L_X=1.0
export L_Y=$(bc<<< "scale=6 ; $MY_NY / $MY_NX")
export L_Z=$(bc<<< "scale=6 ; $MY_NZ / $MY_NX")
#-- global size
export SIZE_X=$(($SIZE_PER_PROC*$MY_NX))
export SIZE_Y=$(($SIZE_PER_PROC*$MY_NY))
export SIZE_Z=$(($SIZE_PER_PROC*$MY_NZ))
#-- 1 thread
export MY_NTHREADS=1
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh

#-- 4 thread
export WT='00:40:00'
export MY_NX=$(bc<<< "scale=6 ; ${MY_NX} ")
export MY_NY=$(bc<<< "scale=6 ; ${MY_NY} / 2")
export MY_NZ=$(bc<<< "scale=6 ; ${MY_NZ} / 2")
export MY_NTHREADS=4
export MY_NTASKS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./marenostrum_kernel_valid.sh


    echo 'still to do'



fi