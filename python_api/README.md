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

### 1. Prerequisites

First, ensure you have the FLUPS C++ library compiled:

```bash
cd /path/to/flups
ARCH_FILE=/your/arch/file make python
```

**Required dependencies:**
- Python >= 3.7
- numpy >= 1.18.0
- mpi4py >= 3.0.0

### 2. Set up Python environment (optional but recommended)

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Upgrade pip
pip install --upgrade pip
```

### 3. Install PyFLUPS

```bash
# Navigate to the python_api directory
cd python_api

# Install dependencies
pip install numpy mpi4py

# Install PyFLUPS in development mode (recommended for testing)
pip install -e .

# Or install normally
pip install .
```

### 4. Set library path

Add to your `.bashrc` or `.zshrc`:

```bash
export FLUPS_LIB_PATH=/path/to/flups/lib
```

Or set it temporarily:
```bash
export FLUPS_LIB_PATH=/path/to/flups/lib
```

`make python` builds one library per SwitchTopo communication backend (`libflups_a2a`, `libflups_nb`,
`libflups_isr`). PyFLUPS loads `libflups` if present, then `libflups_nb`, `libflups_a2a` and `libflups_isr`,
searching `FLUPS_LIB_PATH` first, then the `lib` folder of the FLUPS checkout it belongs to. Set `FLUPS_LIB` to
choose: either a backend name (`export FLUPS_LIB=a2a`) or the full path of a library.

## Testing the Installation

### Quick Test

Create a file `test_import.py`:

```python
import pyflups as pf
print("PyFLUPS version:", pf.__version__)
print("Available classes:", dir(pf))
```

Run:
```bash
python test_import.py
```

### MPI Test

```bash
cd examples
mpirun -np 4 python simple_poisson.py
```

## Troubleshooting

### Library not found

**Error:** `ImportError: Could not find the FLUPS library ...` when importing pyflups. The message lists the
library names and the folders that were searched.

**Solution:**
- Verify the library exists: `ls $FLUPS_LIB_PATH` should show `libflups_a2a.so`, `libflups_nb.so`, `libflups_isr.so`
  (`.dylib` on macOS)
- Ensure FLUPS is compiled with `make python`, e.g. `ARCH_FILE=make_arch/make.osx_gcc make python` (for macOS)
- Point `FLUPS_LIB_PATH` to the folder, or `FLUPS_LIB` to the library itself

### MPI errors

Make sure `mpi4py` is installed and compatible with your MPI installation.

## Verifying the Setup

Run all examples:
```bash
cd examples

# Simple Poisson
mpirun -np 4 python simple_poisson.py

# Vector field
mpirun -np 4 python vector_field.py

# Advanced usage
mpirun -np 4 python advanced_usage.py
```

Expected output: Each script should complete without errors and display results.

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
- `get_stencil(stencil_type=StencilType.LHS)`: Stencil inverted exactly by the lattice kernels, as `(offsets, weights)`
- `do_copy`, `do_fft`, `do_mult`, `get_wavenumbers`, `get_inner_topo_spectral`: the steps of `solve`, see [Advanced Usage](#advanced-usage)

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
- `StencilType`: LHS ($L_h$), RHS ($M_h$)

### Utility Functions

- `hdf5_dump(topology, filename, data)`: Export data to HDF5
- `print_data(topology, data)`: Print data to console
- `locID(axsrc, i0, i1, i2, lia, axtrg, size, nf)`: Compute memory index
- `field_view(topology, data)`: View of `data` in the layout of `topology`, indexed `[component, i0, i1, i2]` (x, y, z order, padding dropped, complex if the topology is)

## Advanced Usage

### Profiling

```python
# Create profiler
prof = pf.Profiler("my_solver")

# Pass to solver
solver = pf.Solver(..., profiler=prof)

# Display results
prof.display()
```

### Step by Step Solve and Spectral Operations

`solve(field, rhs)` is a sequence of five steps on the inner buffer of FLUPS, which can be called one by one:

```python
solver.do_copy(topo, rhs, pf.FLUPS_FORWARD)        # rhs -> inner buffer (host -> device on a GPU build)
solver.do_fft(sign=pf.FLUPS_FORWARD)               # physical -> spectral, in place
solver.do_mult(solver_type=pf.SolverType.STD)      # Green's function and FFT normalisation
solver.do_fft(sign=pf.FLUPS_BACKWARD)              # spectral -> physical, in place
solver.do_copy(topo, field, pf.FLUPS_BACKWARD)     # inner buffer -> field
```

`do_fft` and `do_mult` work on the inner buffer by default. Between the two FFTs, the inner buffer holds the spectrum
in the layout of `solver.get_inner_topo_spectral()`: real if every direction is symmetric, complex as soon as one is
periodic. To modify it (filter, other Green's function, ...):

```python
topo_spec = solver.get_inner_topo_spectral()
spec = np.zeros(topo_spec.memsize)
solver.do_copy(topo_spec, spec, pf.FLUPS_BACKWARD)            # inner buffer -> spec
kx, ky, kz = solver.get_wavenumbers()                          # broadcastable to the local spectral grid
pf.field_view(topo_spec, spec)[...] *= np.exp(-0.5 * sigma**2 * (kx**2 + ky**2 + kz**2))
solver.do_copy(topo_spec, spec, pf.FLUPS_FORWARD)             # spec -> inner buffer
```

- The topology of the spectral data, `field_view` and `get_wavenumbers` are only valid between
  `do_fft(sign=FLUPS_FORWARD)` and `do_fft(sign=FLUPS_BACKWARD)` (`get_wavenumbers` raises otherwise).
- `get_wavenumbers` returns the magnitudes $|k_d|$ used by the Green's functions: fine for isotropic multipliers, not
  for derivatives (use `SolverType.ROT`).
- The FFTs are unnormalised. A multiplier applied after `do_mult` needs nothing more; a multiplier that replaces
  `do_mult` must include the normalisation factor (measured once by a forward + backward FFT of a constant field in
  `examples/manual_solve.py`).
- On an FFTW build, `pf.field_view(topo_spec, solver.get_inner_buffer())` is a view of the inner buffer itself and
  avoids the two copies. On a GPU build the inner buffer lives on the device: use `do_copy`.

See `examples/manual_solve.py` for a complete, validated example (step by step solve, spectral filter, screened Poisson
equation).

### Non-homogeneous Dirichlet / Neumann Walls

FLUPS solves homogeneous symmetric problems (ODD: $u = 0$, EVEN: $\partial_n u = 0$). The lattice kernels
(`LGF_*`, `MEHR_*`) invert a known stencil exactly, so a wall with data ($u = g$ or $\partial_n u = q$) is imposed by
correcting the right-hand side (derivation in `doc/nonhomogeneous_conditions.md`). `pyflups.bc` provides the steps for
global arrays indexed `[x, y, z]`, node-centred by default (pass `center=pf.CenterType.CELL_CENTER` to every call
for cell-centred data):

```python
from pyflups import bc as pbc

solver.set_green_type(pf.GreenType.LGF_4)
lhs = solver.get_stencil(pf.StencilType.LHS)        # L_h (weights in units of 1/h^2)
rhs_op = solver.get_stencil(pf.StencilType.RHS)     # M_h (identity for LGF_*, apply it yourself for MEHR_*L)
r = max(pbc.stencil_reach(lhs[0]), pbc.stencil_reach(rhs_op[0]))

# D = u - R(u) on the grid extended by r ghost layers, from the wall data: ghost(d, side, m) returns
# g (m = 0, ODD walls of node-centred data) or the ghost data of layer m, e.g. pbc.taylor_ghost(a, m, h, wall_type)
D = pbc.ghost_remainder(bc, n, r, ghost)
rhs = pbc.corrected_rhs(*lhs, pbc.apply_stencil(*rhs_op, f_ext, r), D, h, bc)

# ... scatter rhs to the FLUPS layout, solve (ODD for Dirichlet, EVEN for Neumann walls), gather v ...
u = pbc.recover(v, D, bc)                           # sets the Dirichlet wall nodes to g
```

See `examples/nonhomogeneous_bc.py` for a complete, validated example.

## Examples

Complete examples are provided in the `examples/` directory:

- `simple_poisson.py`: Basic Poisson solver
- `vector_field.py`: Solving for vector fields
- `advanced_usage.py`: Using advanced features
- `manual_solve.py`: The steps of `solve` one by one, with a spectral filter and a user Green's function in between (validated against analytic solutions)
- `nonhomogeneous_bc.py`: Non-homogeneous Dirichlet/Neumann walls through a right-hand-side correction (`--suite` runs the validation)

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
@article{caprace2021,
  title={{FLUPS}: A {F}ourier-Based Library of Unbounded {P}oisson Solvers},
  author={Caprace, Denis-Gabriel and Gillis, Thomas and Chatelain, Philippe},
  journal={SIAM Journal on Scientific Computing},
  volume={43},
  number={1},
  pages={C31--C60},
  year={2021},
  doi={10.1137/19M1303848}
}
```

and, for the parallel FFT library and its communication backends:

```
@article{balty2023,
  title={{FLUPS}: A Flexible and Performant Massively Parallel {F}ourier Transform Library},
  author={Balty, Pierre and Chatelain, Philippe and Gillis, Thomas},
  journal={IEEE Transactions on Parallel and Distributed Systems},
  volume={34},
  number={7},
  pages={2011--2024},
  year={2023},
  doi={10.1109/TPDS.2023.3254302}
}
```

## Support

- **Issues**: https://github.com/vortexlab-uclouvain/flups/issues
- **Documentation**: See `doc/` directory in main FLUPS repository

## Authors

- Pierre Balty
- Thomas Gillis
- Denis-Gabriel Caprace
- Contributors from UCLouvain Vortex Lab
