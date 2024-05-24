#!/bin/bash
# Submission script for Zenobe 
#PBS -N scalability
#PBS -q large
#PBS -r y 
#PBS -W group_list=examples
#PBS -l walltime=00:30:00 

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
NCPUS="$( cat $PBS_NODEFILE | wc -l )"
echo "CPU         :  $NCPUS"

module purge
module load compiler/gcc/7.2.0
module load compiler/intel/comp_and_lib/2019.3.199
module load intelmpi/5.1.3.181/64



EXEC_FLUPS=flups_validation

echo
echo "----------------- Start computation -------------"
echo "Starting time : " $(date)

echo "launching mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} valgrind --leak-check=full ./${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZE} -nres 1 -ns 100 -k 0 >> stdout_${PBS_JOBID} "
#mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} valgrind --leak-check=full ./${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZE} -nres 1 -ns 100 -k 0 >> stdout_${PBS_JOBID}
mpirun -n ${NCPUS} -genv OMP_NUM_THREADS=${MY_NTH} ./${EXEC_FLUPS} -np ${MY_NX} ${MY_NY} ${MY_NZ} -res ${MY_SIZE} -nres 1 -ns 250 -k 0 >> stdout_${PBS_JOBID}


################## 
echo "End time : " $(date)
echo "----------------- Computation over, bye bye! ----"




qstat -f $PBS_JOBID
