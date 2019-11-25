#!/bin/sh
## RM the previous validation dir

export HOME_FLUPS=/home/acad/ucl-tfl/tgillis/flups/samples/validation

## fixed parameters
export SIZE_TOT=2048

################################################################################################################
#           NON-BLOCKING COMMUNICATION
################################################################################################################
export EXEC_FLUPS=flups_validation_nb
SCRATCH=/SCRATCH/acad/examples/tgillis/flups_strong_nb

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/zenobe_kernelScaling.sh $SCRATCH
cd $SCRATCH

#================== 288 CPU's ================
#-- proc domain
export NX=6
export NY=6
export NZ=8
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 576 CPU's ================
#-- proc domain
export NX=12
export NY=6
export NZ=8
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 1152 CPU's ================
#-- proc domain
export NX=12
export NY=12
export NZ=8
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY/2")
export NZ2=$NZ
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 2304 CPU's ================
#-- proc domain
export NX=12
export NY=12
export NZ=16
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 4608 CPU's ================
#-- proc domain
export NX=12
export NY=24
export NZ=16
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX")
export NY2=$(bc<<< "scale=0 ; $NY/2")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 6912 CPU's ================
#-- proc domain
export NX=18
export NY=24
export NZ=16
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY/2")
export NZ2=$(bc<<< "scale=0 ; $NZ")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh


################################################################################################################
#           ALL TO ALL COMMUNICATION
################################################################################################################
export EXEC_FLUPS=flups_validation_a2a
SCRATCH=/SCRATCH/acad/examples/tgillis/flups_strong_a2a

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/zenobe_kernelScaling.sh $SCRATCH
cd $SCRATCH

#================== 288 CPU's ================
#-- proc domain
export NX=6
export NY=6
export NZ=8
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 576 CPU's ================
#-- proc domain
export NX=12
export NY=6
export NZ=8
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 1152 CPU's ================
#-- proc domain
export NX=12
export NY=12
export NZ=8
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY/2")
export NZ2=$(bc<<< "scale=0 ; $NZ")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 2304 CPU's ================
#-- proc domain
export NX=12
export NY=12
export NZ=16
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 4608 CPU's ================
#-- proc domain
export NX=12
export NY=24
export NZ=16
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX")
export NY2=$(bc<<< "scale=0 ; $NY/2")
export NZ2=$(bc<<< "scale=0 ; $NZ/2")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh

#================== 6912 CPU's ================
#-- proc domain
export NX=18
export NY=24
export NZ=16
#-- domain length
export LX=1.0
export LY=1.0
export LZ=1.0
#-- global size
export SIZE_X=$SIZE_TOT
export SIZE_Y=$SIZE_TOT
export SIZE_Z=$SIZE_TOT
#-- 1 thread
export N_NODE=$(bc<<< "scale=0 ; ($NX*$NY*$NZ)/24")
echo "qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX},MY_NY=${NY},MY_NZ=${NZ},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=1,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernelScaling.sh
#-- 4 thread
export NX2=$(bc<<< "scale=0 ; $NX/2")
export NY2=$(bc<<< "scale=0 ; $NY/2")
export NZ2=$(bc<<< "scale=0 ; $NZ")
echo "qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=${MY_NTHREADS},L_X=${LX},L_Y=${LY},L_Z=${LZ} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./EXEC_FLUPS.sh"
qsub -q R4661004 -v MY_NX=${NX2},MY_NY=${NY2},MY_NZ=${NZ2},MY_SIZE_X=${SIZE_X},MY_SIZE_Y=${SIZE_Y},MY_SIZE_Z=${SIZE_Z},MY_NTH=4,L_X=${LX},L_Y=${LY},L_Z=${LZ},EXEC_FLUPS=${EXEC_FLUPS} -l select=${N_NODE}:ncpus=24:mem=63000mb:mpiprocs=6:ompthreads=4 ./zenobe_kernelScaling.sh


#end of file
