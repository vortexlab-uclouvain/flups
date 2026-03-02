# PyFLUPS Examples

This directory contains example scripts demonstrating various features of PyFLUPS.

## Examples

### 1. simple_poisson.py
Basic usage of PyFLUPS to solve a simple 3D Poisson equation with a point source.

**Features demonstrated:**
- Creating topology
- Setting up solver
- Solving standard Poisson equation
- Computing statistics

**Run:**
```bash
mpirun -np 4 python simple_poisson.py
```

### 2. vector_field.py
Solving for vector fields using the rotational (Biot-Savart) formulation.

**Features demonstrated:**
- 3-component vector fields (lda=3)
- ROT solver type
- Spectral derivatives
- Creating vortex ring vorticity distribution
- Computing velocity from vorticity

**Run:**
```bash
mpirun -np 4 python vector_field.py
```

### 3. advanced_usage.py
Advanced features and low-level operations.

**Features demonstrated:**
- Using profiler for timing
- Getting spectral information
- Accessing internal buffers
- Manual FFT operations
- Comparing different Green functions
- Memory layout and index computation

**Run:**
```bash
mpirun -np 4 python advanced_usage.py
```

## Requirements

All examples require:
- FLUPS library compiled and accessible
- mpi4py
- numpy
- MPI runtime (e.g., OpenMPI, MPICH)

Set `FLUPS_LIB_PATH` environment variable:
```bash
export FLUPS_LIB_PATH=/path/to/flups/lib
```

## Notes

- Adjust the number of MPI processes (`-np`) based on your system
- Some examples may generate output files in a `data/` directory
- Examples use periodic boundary conditions by default
- Grid sizes can be adjusted in each script
