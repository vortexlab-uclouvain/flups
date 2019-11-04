#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=${PWD}/../


VER=a2a

EXE=flups_vs_p3dfft++_${VER}

######### WEAK -> increase the number of CPU and the size
SCRATCH=/SCRATCH/acad/examples/dcaprace/flupsVSp3dfft3_weak_${VER}_V3

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
## CANNOT DO cpu=96,192,384... due to P3D !

# cpu = 128 (-> actually allocating 144)
# same on large
qsub -q large -v EXE=${EXE},MY_NY=8,MY_NZ=16,LX=4,LY=4,LZ=8,MY_SIZE=128,MY_NTH=1, -l select=6:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernel.sh

# cpu = 256 (->264)
qsub -q large -v EXE=${EXE},MY_NY=16,MY_NZ=16,LX=4,LY=8,LZ=8,MY_SIZE=128,MY_NTH=1, -l select=11:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernel.sh

# cpu = 512 (->528)
qsub -q large -v EXE=${EXE},MY_NY=16,MY_NZ=32,LX=8,LY=8,LZ=8,MY_SIZE=128,MY_NTH=1, -l select=22:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernel.sh

# cpu = 1024 (->1032)
qsub -q large -v EXE=${EXE},MY_NY=32,MY_NZ=32,LX=8,LY=8,LZ=16,MY_SIZE=128,MY_NTH=1, -l select=43:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernel.sh

# cpu = 2048 (->2064)
qsub -q large -v EXE=${EXE},MY_NY=32,MY_NZ=64,LX=8,LY=16,LZ=16,MY_SIZE=128,MY_NTH=1, -l select=86:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 ./zenobe_kernel.sh

#end of file
