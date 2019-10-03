#!/bin/sh
## RM the previous validation dir

HOME_FLUPS=/home/acad/ucl-tfl/tgillis/flups
EXEC_FLUPS=flups_validation_nb

######### STRONG -> increase the number of CPU, same size
SCRATCH=/SCRATCH/acad/examples/tgillis/flups_weak_nb

# clean the validation dir
rm -rf ${SCRATCH}
mkdir -p $SCRATCH
mkdir -p $SCRATCH/data
mkdir -p $SCRATCH/prof
# copy the needed info
cp $HOME_FLUPS/$EXEC_FLUPS $SCRATCH
cp $HOME_FLUPS/run/zenobe_kernel_nb.sh $SCRATCH

cd $SCRATCH

######################################################
# 64 N per CPU
# n = 1 - size = 64 - cpu = 1
qsub -q main -v MY_NX=1,MY_NY=1,MY_NZ=1,MY_SIZE=64,MY_NTH=1 -l select=1:ncpus=1:mem=10500mb:mpiprocs=1:ompthreads=1 ./zenobe_kernel_nb.sh
# n = 2 - size = 128 - cpu = 8
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=2,MY_SIZE=128,MY_NTH=1 -l select=2:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=1,MY_SIZE=128,MY_NTH=2 -l select=2:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
# n = 3 - size = 256 - cpu = 64
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=4,MY_SIZE=256,MY_NTH=1 -l select=16:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=2,MY_SIZE=256,MY_NTH=2 -l select=16:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=1,MY_SIZE=256,MY_NTH=4 -l select=16:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh
# n = 4 - size = 512 - cpu = 512
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=8,MY_SIZE=512,MY_NTH=1 -l select=128:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=4,MY_SIZE=512,MY_NTH=2 -l select=128:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=2,MY_SIZE=512,MY_NTH=4 -l select=128:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh
# n = 5 - size = 1024 - cpu = 4096
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=16,MY_SIZE=1024,MY_NTH=1 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=8,MY_SIZE=1024,MY_NTH=2 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=4,MY_SIZE=1024,MY_NTH=4 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

######################################################
# 128 N per CPU
# n = 1 - size = 128 - cpu = 1
qsub -q main -v MY_NX=1,MY_NY=1,MY_NZ=1,MY_SIZE=128,MY_NTH=1 -l select=1:ncpus=1:mem=10500mb:mpiprocs=1:ompthreads=1 ./zenobe_kernel_nb.sh
# n = 2 - size = 256 - cpu = 8
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=2,MY_SIZE=256,MY_NTH=1 -l select=2:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=1,MY_SIZE=256,MY_NTH=2 -l select=2:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
# n = 3 - size = 512 - cpu = 64
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=4,MY_SIZE=512,MY_NTH=1 -l select=16:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=2,MY_SIZE=512,MY_NTH=2 -l select=16:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=1,MY_SIZE=512,MY_NTH=4 -l select=16:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh
# n = 4 - size = 1024 - cpu = 512
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=8,MY_SIZE=1024,MY_NTH=1 -l select=128:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=4,MY_SIZE=1024,MY_NTH=2 -l select=128:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=2,MY_SIZE=1024,MY_NTH=4 -l select=128:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh
# n = 5 - size = 2048 - cpu = 4096
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=16,MY_SIZE=2048,MY_NTH=1 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=8,MY_SIZE=2048,MY_NTH=2 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=4,MY_SIZE=2048,MY_NTH=4 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

#######################################################
## 192 N per CPU
## n = 1 - size = 128 - cpu = 1
#qsub -q main -v MY_NX=1,MY_NY=1,MY_NZ=1,MY_SIZE=192,MY_NTH=1 -l select=1:ncpus=1:mem=10500mb:mpiprocs=1:ompthreads=1 ./zenobe_kernel_nb.sh
## n = 2 - size = 256 - cpu = 8
#qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=2,MY_SIZE=384,MY_NTH=1 -l select=2:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
#qsub -q main -v MY_NX=2,MY_NY=2,MY_NZ=1,MY_SIZE=384,MY_NTH=2 -l select=2:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
## n = 3 - size = 512 - cpu = 64
#qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=4,MY_SIZE=768,MY_NTH=1 -l select=16:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
#qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=2,MY_SIZE=768,MY_NTH=2 -l select=16:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
#qsub -q main -v MY_NX=4,MY_NY=4,MY_NZ=1,MY_SIZE=768,MY_NTH=4 -l select=16:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh
## n = 4 - size = 1024 - cpu = 512
#qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=8,MY_SIZE=1536,MY_NTH=1 -l select=128:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
#qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=4,MY_SIZE=1536,MY_NTH=2 -l select=128:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
#qsub -q large -v MY_NX=8,MY_NY=8,MY_NZ=2,MY_SIZE=1536,MY_NTH=4 -l select=128:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh
## n = 5 - size = 2048 - cpu = 4096
#qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=16,MY_SIZE=3072,MY_NTH=1 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=4:ompthreads=1 ./zenobe_kernel_nb.sh
#qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=8,MY_SIZE=3072,MY_NTH=2 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=2:ompthreads=2 ./zenobe_kernel_nb.sh
#qsub -q large -v MY_NX=16,MY_NY=16,MY_NZ=4,MY_SIZE=3072,MY_NTH=4 -l select=1024:ncpus=4:mem=10500mb:mpiprocs=1:ompthreads=4 ./zenobe_kernel_nb.sh

#end of file
