#!/bin/bash
# Submission script for Zenobe 
#PBS -N conv3D
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
NCPUS=$((${MY_NX}*${MY_NY}*${MY_NZ}))
echo "CPU         :  $NCPUS"

#module purge
#module load compiler/gcc/7.2.0
#module load compiler/intel/comp_and_lib/2017.4.196
#module load intelmpi/2017.3.196/64
#module load devtoolset/8

module load OpenMPI/3.1.4-GCC-8.3.0
module load HDF5/1.10.5-gompi-2019b
module load FFTW/3.3.8-gompi-2019b

EXEC_FLUPS=flups_validation_a2a

echo
echo "----------------- Start computation -------------"
echo "Starting time : " $(date)

MY_SIZE_X=$((${MY_SIZE} * ${L_X}/${L_X}))
MY_SIZE_Y=$((${MY_SIZE} * ${L_Y}/${L_X}))
MY_SIZE_Z=$((${MY_SIZE} * ${L_Z}/${L_X}))

########################## --bc=4,4,4,4,4,4 ###########################
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=7 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=7 --bc=4,4,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}



######################### --bc=0,0,1,0,3,3 ###########################
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=0,0,1,0,3,3 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}


######################### --bc=4,0,4,4,1,4 ###########################
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"    
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}                               
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"    
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}                              
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=7 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"    
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=7 --bc=4,0,4,4,1,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}     

######################### --bc=3,3,4,4,4,4 ###########################
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=0 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
# echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
# env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=1 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=2 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=3 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=4 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=5 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}
echo "launching  env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}"
env OMP_NUM_THREADS=${MY_NTH} mpirun -n ${NCPUS} -s   ./${EXEC_FLUPS} --np=${MY_NX},${MY_NY},${MY_NZ} --res=${MY_SIZE_X},${MY_SIZE_Y},${MY_SIZE_Z} --dom=${L_X},${L_Y},${L_Z} --nres=1 --kernel=6 --bc=3,3,4,4,4,4 --center=0 --outdir=${PBS_O_WORKDIR} >> stdout_${PBS_JOBID}

################## 
echo "End time : " $(date)
echo "----------------- Computation over, bye bye! ----"




qstat -f $PBS_JOBID
