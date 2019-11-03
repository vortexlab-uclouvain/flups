#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=${PWD}/../


VER=a2a

EXE=flups_vs_p3dfft++_${VER}

######### WEAK -> increase the number of CPU and the size
SCRATCH=/SCRATCH/acad/examples/dcaprace/flupsVSp3dfft3_weak_$VER

# clean the validation dir
# rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXE $SCRATCH
cp $HOME_FLUPS/${EXE}_noP3D $SCRATCH
cp $HOME_FLUPS/run/zenobe_kernel.sh $SCRATCH

cd $SCRATCH

#####################   size = 64^3/proc  #################################
# cpu = 64
# same on large
# qsub -q large -v EXE=${EXE},MY_NY=8,MY_NZ=16,LX=4,LY=4,LZ=8,MY_SIZE=64,MY_NTH=1, -l select=32:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

# cpu = 256
# qsub -q large -v EXE=${EXE},MY_NY=16,MY_NZ=16,LX=4,LY=8,LZ=8,MY_SIZE=64,MY_NTH=1, -l select=64:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

# cpu = 512
# qsub -q large -v EXE=${EXE},MY_NY=16,MY_NZ=32,LX=8,LY=8,LZ=8,MY_SIZE=64,MY_NTH=1, -l select=128:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

# cpu = 1024
# qsub -q large -v EXE=${EXE},MY_NY=32,MY_NZ=32,LX=8,LY=8,LZ=16,MY_SIZE=64,MY_NTH=1, -l select=256:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

#####################   size = 128^3/proc  #################################
# cpu = 64
# same on large
qsub -q large -v EXE=${EXE},MY_NY=8,MY_NZ=16,LX=4,LY=4,LZ=8,MY_SIZE=128,MY_NTH=1, -l select=32:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

# cpu = 256
qsub -q large -v EXE=${EXE},MY_NY=16,MY_NZ=16,LX=4,LY=8,LZ=8,MY_SIZE=128,MY_NTH=1, -l select=64:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

# cpu = 512
qsub -q large -v EXE=${EXE},MY_NY=16,MY_NZ=32,LX=8,LY=8,LZ=8,MY_SIZE=128,MY_NTH=1, -l select=128:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh

# cpu = 1024
qsub -q large -v EXE=${EXE},MY_NY=32,MY_NZ=32,LX=8,LY=8,LZ=16,MY_SIZE=128,MY_NTH=1, -l select=256:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel.sh


#end of file
