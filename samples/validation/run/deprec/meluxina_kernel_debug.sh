#!/bin/bash -l
#SBATCH --account=p200053
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --ntasks-per-node=128
#SBATCH --time=00:30:00
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=%x_%j_%4t.out
#SBATCH --error=%x_%j_%4t.err

#--------------------------------------------------------
# module use /apps/USE/easybuild/staging/2021.1/modules/all
source ${MODULES} ${OMPIVERSION}

module load Valgrind

cd ${SCRATCH_FLUPS}
cp ${FLUPS_DIR}/samples/validation/${EXEC_FLUPS} ${SCRATCH_FLUPS}


echo "----------------- launching job -----------------"
echo "OMP_NUM_THREADS=1 srun ./${EXEC_FLUPS} --np=128,4,1 --bc=4,4,4,4,3,3 --center=1"

OMP_NUM_THREADS=1 mpirun -np 256 --tag-output --output-filename log0 valgrind --leak-check=full --track-origins=yes --show-reachable=yes ./${EXEC_FLUPS} --np=128,2,1 --bc=4,4,4,4,3,3 --res=384,384,1536 --dom=1,1,4 --center=1 --kernel=0 --nres=1 --nsolve=20
OMP_NUM_THREADS=1 mpirun -np 256 --tag-output --output-filename log1 valgrind --leak-check=full --track-origins=yes --show-reachable=yes ./${EXEC_FLUPS} --np=128,2,1 --bc=4,4,4,4,3,3 --res=384,384,1536 --dom=1,1,4 --center=0 --kernel=0 --nres=1 --nsolve=20

OMP_NUM_THREADS=1 mpirun -np 256 --tag-output --output-filename log2 ./${EXEC_FLUPS} --np=128,2,1 --bc=4,4,4,4,3,3 --res=384,384,1536 --dom=1,1,4 --center=1 --kernel=0 --nres=1 --nsolve=20
OMP_NUM_THREADS=1 mpirun -np 256 --tag-output --output-filename log3 ./${EXEC_FLUPS} --np=128,2,1 --bc=4,4,4,4,3,3 --res=384,384,1536 --dom=1,1,4 --center=0 --kernel=0 --nres=1 --nsolve=20
cd -
