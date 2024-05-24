#!/bin/bash
# Submission script for JUWELS
#SBATCH --account=prpa79
#SBATCH --job-name=scaling
#
#SBATCH --output=flups_%j.out
#SBATCH --error=flups_%j.err

export OMP_NUM_THREADS=${MY_NTHREADS}

export I_MPI_DEBUG=+5

echo "----------------- Load modules -----------------"
module purge
#module load Intel/2019.3.199-GCC-8.3.0
#module load IntelMPI/2018.5.288
#module load IntelMPI/2019.3.199

#module load intel-para/2019a
##CHANGING COMMUNICATION METHOD, OTHERWISE MPI USE TOO MUCH MEMORY !
#export PSP_UCP=1
#export UCX_TLS=ud_mlx5,self,sm
##export PSI_LOGGERDEBUG=1
##export PSI_FORWARDERDEBUG=1

#OLD DEVELOPMENT IMPI
# module use /p/software/juwels/otherstages/
# module load Stages/Devel-2019a
# module load Intel
# module load IntelMPI/2019.6.RC20191024
#----------------------
module use $OTHERSTAGES
module load Stages/Devel-2019a
module load Intel
module load IntelMPI/2019.6.154

module load FFTW/3.3.8
module load HDF5/1.10.5
module load METIS/5.1.0
module list


echo "----------------- launching job -----------------"
echo "launch command: srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${SIZE_X} ${SIZE_Y} ${SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0"

srun --label ${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${SIZE_X} ${SIZE_Y} ${SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0

scontrol show job ${SLURM_JOB_ID}

sacct --format="JobID,NCPUS,NNodes,Elapsed,MaxRSS,MaxVMSize,ExitCode" | grep "$SLURM_JOB_ID"