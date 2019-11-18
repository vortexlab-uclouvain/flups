#!/bin/sh

HOME_FLUPS=/p/project/prpa79/flups/samples/validation
KERNEL=juwels_kernel_valid.sh

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

export nPerSwitch=1152 #number of process per switch
export SW_TIMEOUT=1440 #minutes to wait before releasing the constraint on switches


if [ "$ver" = "small" ]; then
############################################################
############################################################
############################################################
#                           SMALL (<=4k)
############################################################
############################################################
############################################################


############################################################
# ALL TO ALL
#-----------------------------------------------------------
export EXEC_FLUPS=flups_validation_a2a

SCRATCH=/p/scratch/prpa79/$(whoami)/flups_weak_a2a_${ver}

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof

# copy the needed info
cp $HOME_FLUPS/${EXEC_FLUPS} $SCRATCH
cp $HOME_FLUPS/run/${KERNEL} $SCRATCH
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


#================== 4608 CPU's ================
#-- requested walltime
export WT='00:10:00'
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}




############################################################
# NON-BLOCKING
#-----------------------------------------------------------
export EXEC_FLUPS=flups_validation_nb

SCRATCH=/p/scratch/prpa79/$(whoami)/flups_weak_nb_${ver}

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof

# copy the needed info
cp $HOME_FLUPS/${EXEC_FLUPS} $SCRATCH
cp $HOME_FLUPS/run/${KERNEL} $SCRATCH
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


#================== 4608 CPU's ================
#-- requested walltime
export WT='00:10:00'
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


elif [ "$ver" = "large" ]; then

############################################################
############################################################
############################################################
#                           LARGE (>4k, <=18k)
############################################################
############################################################
############################################################



# ############################################################
# # ALL TO ALL
# #-----------------------------------------------------------
export EXEC_FLUPS=flups_validation_a2a

SCRATCH=/p/scratch/prpa79/$(whoami)/flups_weak_a2a_${ver}

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof

# copy the needed info
cp $HOME_FLUPS/${EXEC_FLUPS} $SCRATCH
cp $HOME_FLUPS/run/${KERNEL} $SCRATCH
# go to it
cd $SCRATCH


#================== 9216 CPU's ================
#-- requested walltime
export WT='00:20:00'
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


#================== 18432 CPU's ================
#-- requested walltime
export WT='00:20:00'
#-- proc domain
export MY_NX=32
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}




############################################################
# NON-BLOCKING
#-----------------------------------------------------------
export EXEC_FLUPS=flups_validation_nb

SCRATCH=/p/scratch/prpa79/$(whoami)/flups_weak_nb_${ver}

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof

# copy the needed info
cp $HOME_FLUPS/${EXEC_FLUPS} $SCRATCH
cp $HOME_FLUPS/run/${KERNEL} $SCRATCH
# go to it
cd $SCRATCH


#================== 9216 CPU's ================
#-- requested walltime
export WT='00:20:00'
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


#================== 18432 CPU's ================
#-- requested walltime
export WT='00:20:00'
#-- proc domain
export MY_NX=32
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
export N_SWITCH=$(bc<<< "scale=0 ; $MY_NTASKS / $nPerSwitch") #CAUTION, THIS WORKS ONLY BECAUSE WE ALWAYS HAVE A MULTIPLE OF 24 SWITCHES
# echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT}  -- ./${KERNEL}"
# sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} --switches=${N_SWITCH}@${SW_TIMEOUT} ./${KERNEL}
echo "sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT} ./${KERNEL}"
sbatch --ntasks=${MY_NTASKS} --cpus-per-task=${MY_NTHREADS} --time=${WT}  ./${KERNEL}


elif [ "$ver" = "Xlarge" ]; then

############################################################
############################################################
############################################################
#                     EXTRA  LARGE (>18k)
############################################################
############################################################
############################################################

echo "must be done"

fi
