#!/bin/bash
# Submission script for Zenobe 
#PBS -N scalability
#PBS -r y 
#PBS -W group_list=examples
#PBS -l walltime=00:15:00 

exec > ${PBS_O_WORKDIR}/${PBS_JOBNAME}_${PBS_JOBID}.log 
echo "------------------ Work dir --------------------" 
cd ${PBS_O_WORKDIR} && echo ${PBS_O_WORKDIR} 
echo "------------------ Job Info --------------------" 
echo "jobid : $PBS_JOBID" 
echo "jobname : $PBS_JOBNAME" 
echo "job type : $PBS_ENVIRONMENT" 
echo "submit dir : $PBS_O_WORKDIR" 
echo "queue : $PBS_O_QUEUE" 
echo "user : $PBS_O_LOGNAME" 
echo "threads : $OMP_NUM_THREADS -> ${MY_NTH}" 
NCPUS=$((1*${MY_NY}*${MY_NZ}))
echo "CPU         :  $NCPUS"

module purge
module load compiler/intel/comp_and_lib/2017.4.196 intelmpi/2017.3.196/64 compiler/gcc/7.2.0
module load metis/5.1.0/32/intel/2016.2.181
#module load devtoolset/8


echo
echo "----------------- Start computation -------------"
echo "Starting time : " $(date)

MY_SIZE_X=$((${MY_SIZE}*${LX}))
MY_SIZE_Y=$((${MY_SIZE}*${LY}))
MY_SIZE_Z=$((${MY_SIZE}*${LZ}))

echo "launching  mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} ./${EXE} -np ${MY_NY} ${MY_NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -ni 100 >> stdout_${PBS_JOBID}"
mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} ./${EXE} -np ${MY_NY} ${MY_NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -ni 100 >> stdout_${PBS_JOBID}

################## 
echo "End time : " $(date)
echo "----------------- Computation over, bye bye! ----"

qstat -f $PBS_JOBID
