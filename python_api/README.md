# PyFLUPS - Python API for FLUPS

Python bindings for FLUPS (Fourier-based Library of Unbounded Poisson Solvers).

## Overview

FLUPS is a high-performance library for solving Poisson equations using FFT-based methods. PyFLUPS provides a Python interface to this library, enabling easy integration with Python scientific computing workflows.

## Features

- **Object-oriented Python API** for FLUPS functionality
- **Full MPI support** via mpi4py
- **NumPy integration** for array operations
- Support for various **boundary conditions**: EVEN, ODD, PERiodic, UNBounded
- Multiple **Green's function types**: LGF, HEJ, CHAT, MEHR
- Both **standard** (STD) and **rotational** (ROT) Poisson solvers
- **HDF5 export** capabilities with XDMF support

## Installation

### Prerequisites

1. **FLUPS library**: You must first build the FLUPS C++ library
2. **Python dependencies**:
   - Python >= 3.7
   - numpy >= 1.18.0
   - mpi4py >= 3.0.0

### Build FLUPS

```bash
# Navigate to the FLUPS root directory
cd /path/to/flups

# Build the library (example with default architecture)
make ARCH_FILE=make_arch/make.default
```

### Install PyFLUPS

```bash
# Navigate to the python_api directory
cd python_api

# Install in development mode (recommended for testing)
pip install -e .

# Or install normally
pip install .
```

### Set Library Path

Set the `FLUPS_LIB_PATH` environment variable to point to the compiled FLUPS library:

```bash
export FLUPS_LIB_PATH=/path/to/flups/lib
```

Or add it to your `.bashrc` or `.zshrc` file.

## Quick Start

### Basic Example

```python
import numpy as np
from mpi4py import MPI
import pyflups as pf

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Domain parameters
nglob = [64, 64, 64]  # Global grid size
L = [1.0, 1.0, 1.0]   # Domain size
h = [L[i]/nglob[i] for i in range(3)]  # Grid spacing

# Process decomposition
nproc = [2, 2, 1]  # Number of processes per dimension

# Boundary conditions
bc = [
    [pf.BoundaryType.PER, pf.BoundaryType.PER],  # x: periodic
    [pf.BoundaryType.PER, pf.BoundaryType.PER],  # y: periodic
    [pf.BoundaryType.EVEN, pf.BoundaryType.EVEN] # z: zero flux
]

# Create topology
topo = pf.Topology(
    axis=0,
    lda=1,  # Scalar field
    nglob=nglob,
    nproc=nproc,
    is_complex=False,
    comm=comm
)

# Create solver
solver = pf.Solver(
    topology=topo,
    boundary_conditions=bc,
    h=h,
    L=L,
    order_diff=pf.DiffType.NOD
)

# Set Green's function type
solver.set_green_type(pf.GreenType.LGF_2)

# Setup solver
solver.setup()

# Prepare arrays - use nmem (padded memory size) for allocation
nloc = [topo.get_nloc(i) for i in range(3)]
nmem = [topo.get_nmem(i) for i in range(3)]
size = nmem[0] * nmem[1] * nmem[2]

rhs = np.zeros(size, dtype=np.float64)
field = np.zeros(size, dtype=np.float64)

# Initialize RHS (example: point source at center)
center_idx = pf.locID(0, nloc[0]//2, nloc[1]//2, nloc[2]//2, 0, 
                      topo.axis, nmem, 1)
if rank == 0:
    rhs[center_idx] = 1.0

# Solve
solver.solve(field, rhs, pf.SolverType.STD)

# Export results
if rank == 0:
    print("Solution computed successfully")
```

## API Documentation

### Main Classes

#### `Topology`
Manages domain decomposition and memory layout for parallel FFT operations.

**Constructor:**
```python
Topology(axis, lda, nglob, nproc, is_complex=False, 
         axproc=None, alignment=16, comm=None)
```

#### `Solver`
Main class for solving Poisson equations.

**Constructor:**
```python
Solver(topology, boundary_conditions, h, L, 
       order_diff=DiffType.NOD, center_type=None, profiler=None)
```

**Key Methods:**
- `set_green_type(green_type)`: Set Green's function type
- `set_alpha(alpha)`: Set alpha for Hejlesen kernels
- `setup(change_comm=False)`: Setup and allocate memory
- `solve(field, rhs, solver_type=STD)`: Solve the equation

#### `Profiler`
Timing measurements for performance analysis.

**Constructor:**
```python
Profiler(name="default")
```

### Enumerations

- `BoundaryType`: EVEN, ODD, PER, UNB, NONE
- `GreenType`: CHAT_2, LGF_2/4/6/8, HEJ_0/2/4/6/8/10, MEHR_4L/6L/4F/6F
- `SolverType`: STD (standard), ROT (rotational/Biot-Savart)
- `DiffType`: NOD, SPE, FD2, FD4, FD6
- `CenterType`: NODE_CENTER, CELL_CENTER

### Utility Functions

- `hdf5_dump(topology, filename, data)`: Export data to HDF5
- `print_data(topology, data)`: Print data to console
- `locID(axsrc, i0, i1, i2, lia, axtrg, size, nf)`: Compute memory index

## Advanced Usage

### Using Inner Buffers

```python
# Get FLUPS internal buffer
buffer = solver.get_inner_buffer()

# Get physical topology
topo_phys = solver.get_inner_topo_physical()
```

### Profiling

```python
# Create profiler
prof = pf.Profiler("my_solver")

# Pass to solver
solver = pf.Solver(..., profiler=prof)

# Display results
prof.display()
```

### Spectral Operations

```python
# Get spectral information
kfact, koffset, symstart = solver.get_spectral_info()

# Manual FFT operations
solver.do_copy(topology, data, pf.FLUPS_FORWARD)
solver.do_fft(data, pf.FLUPS_FORWARD)
solver.do_mult(data, pf.SolverType.STD)
solver.do_fft(data, pf.FLUPS_BACKWARD)
solver.do_copy(topology, data, pf.FLUPS_BACKWARD)
```

## Examples

Complete examples are provided in the `examples/` directory:

- `simple_poisson.py`: Basic Poisson solver
- `vector_field.py`: Solving for vector fields
- `advanced_usage.py`: Using advanced features

Run examples with MPI:
```bash
mpirun -np 4 python examples/simple_poisson.py
```

## License

Copyright © UCLouvain 2020

Licensed under the Apache License, Version 2.0. See LICENSE file for details.

## Citation

If you use FLUPS in your research, please cite:

```
@article{gillis2019,
  title={A parallel multigrid-based preconditioner for the 3D heterogeneous Helmholtz equation},
  author={Gillis, Thomas and others},
  journal={Journal of Computational Physics},
  year={2019}
}
```

## Support

- **Issues**: https://github.com/vortexlab-uclouvain/flups/issues
- **Documentation**: See `doc/` directory in main FLUPS repository

## Authors

- Thomas Gillis
- Denis-Gabriel Caprace
- Contributors from UCLouvain Vortex Lab
