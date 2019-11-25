#!/bin/bash
# Submission script for Zenobe 
#PBS -N scalability
#PBS -r y 
#PBS -W group_list=examples
#PBS -l walltime=01:00:00 

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
echo "threads : $OMP_NUM_THREADS" 
NCPUS=$((${MY_NX}*${MY_NY}*${MY_NZ}))
echo "CPU         :  $NCPUS"

module purge
module load compiler/gcc/7.2.0
module load compiler/intel/comp_and_lib/2017.4.196
module load intelmpi/2017.3.196/64
module load devtoolset/8

echo
echo "----------------- Start computation -------------"
echo "Starting time : " $(date)

NCPUS=$(bc<<< "scale=0 ; ($MY_NX*$MY_NY*$MY_NZ)/1")

echo "launching mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} ./${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0 >> stdout_${PBS_JOBID}"
mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} ./${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -L ${L_X} ${L_Y} ${L_Z} -nres 1 -ns 20 -k 0 >> stdout_${PBS_JOBID}


################## 
echo "End time : " $(date)
echo "----------------- Computation over, bye bye! ----"




qstat -f $PBS_JOBID
