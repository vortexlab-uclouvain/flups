#!/bin/bash

##########################################################
#       EXAMPLE SCRIPT TO RUN FLUPS FOR BENCHMARKING
#       
# note: All 3 calls below MUST be executed.
##########################################################

#number of MPI ranks in each direction. The product is equal to the total number of ranks.
NX=8
NY=8
NZ=8
NCPUS=$((${NX}*${NY}*${NZ}))

#resolution per MPI rank (recommanded 128, can be reduced to 64 if not enough memory)
res=128

#number of repetitions of the solve operation (for more reliable statistics of the execution time)
NS=100

EXEC_FLUPS=flups_validation_a2a

###########################################################

export OMP_NUM_THREADS=1 

echo "----------------- Start computation -------------"
echo "Starting time : " $(date)

MY_SIZE_X=$((${NX} * ${res}))
MY_SIZE_Y=$((${NY} * ${res}))
MY_SIZE_Z=$((${NZ} * ${res}))

mpirun -n ${NCPUS} ./${EXEC_FLUPS} -np ${NX} ${NY} ${NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -ns $NS -bc 1 1 1 1 1 1 >> stdout_fullOdd
mpirun -n ${NCPUS} ./${EXEC_FLUPS} -np ${NX} ${NY} ${NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -ns $NS -bc 0 0 0 0 0 0 >> stdout_fullEven
mpirun -n ${NCPUS} ./${EXEC_FLUPS} -np ${NX} ${NY} ${NZ} -res ${MY_SIZE_X} ${MY_SIZE_Y} ${MY_SIZE_Z} -ns $NS -bc 1 4 4 0 4 4 >> stdout_UnbOdd
#to reduce the pressure on the memory, you might consider reducing res by a factor 2 (only for this last test)

echo "End time : " $(date)
echo "----------------- Computation over, bye bye! ----"

