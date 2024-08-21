# Running FLUPS on GPUs

FLUPS can be run on Nvidia GPUs using CUDA. AMD GPUs support is WIP, and should be available soon. To pick the backend you want to run on (CUDA, HIP or FFTW), you can use the `-DBACKEND_*` variable (default is FFTW for CPUs).

## GPU-aware MPI

To run on GPU, a GPU-aware MPI library is required (this allows the MPI to automatically handle GPU-allocated memory without having to explicitely copying it back to the CPU). We have tested FLUPS extensively with MPICH, and thus recommend you to use MPICH.

To install a CUDA-aware version of MPICH, you can follow this [DBS recipe for Lucia](config/lucia/mpich_ofi_cuda_opt.arch). (Note that this recipe is for the Lucia cluster, and will need to be adapted for your own cluster).

As of the time writing these lines, there is a [known issue with MPICH that can cause the GPU memory to overflow](https://github.com/pmodels/mpich/pull/7055). To avoid this, it is recommended to set `MPIR_CVAR_NOLOCAL=1` as to bypass MPICH's shared memory optimizations and use the underlying provider's instead (note: you can also set `MPIR_CVAR_CH4_IPC_GPU_HANDLE_CACHE=disabled` not to use the GPU handle cache that's overflowing).

## Example scripts

#### Lucia

```bash
#SBATCH --job-name=gpu_job
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=48:00:00
#SBATCH --account=examples
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

# DBS-compiled MPICH
export PREFIX=${HOME}/soft/lib-MPICH-4.2.2-OFI-1.21.0-CUDA-12.2.0
export ARCH_FILE=make_arch/make.dbs
module load EasyBuild/2023a
module load CUDA/12.2.0
module load GCC/12.3.0
module load GDRCopy/2.3.1-GCCcore-12.3.0.lua

# output dirs
flups_dir=${HOME}/Codes/flups

# compile the code
export MPICH_CC=nvcc
export MPICH_CXX=nvcc
cd ${flups_dir}
make clean && make install -j

# compile the client
unset MPICH_CC MPICH_CXX
cd ${flups_dir}/test
make clean && make -j

# Export some stuff
export FI_PROVIDER="verbs,shm"
export MPIR_CVAR_NOLOCAL=1

# lezgo
${PREFIX}/bin/mpirun -l --bind-to core ./flups_test_a2a
```
