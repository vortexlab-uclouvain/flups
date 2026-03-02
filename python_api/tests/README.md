# Python API Tests for FLUPS

Test suite mirroring `samples/validation` with periodic Poisson solver validation and analytical solution comparison.

## Test Files

- **`test_validation_simple.py`**: pytest-based test (single rank only)
  - Solves 3D periodic Poisson with `sin(x)+sin(y)+sin(z)` RHS
  - Compares solver output to analytical solution `φ = -rhs`
  - Validates relative error < 1e-2

- **`mpi_validation_simple.py`**: Standalone MPI validation script (no pytest)
  - Automatically decomposes MPI ranks into 3D grid
  - Allocates with `nmem` (padded memory) for FLUPS compatibility
  - Uses `MPI.allreduce` to compute global max error
  - Exits with code 0 on success, 1 on failure

## Requirements

- `numpy` >= 1.20
- `mpi4py` >= 3.0
- `pytest` (for `test_validation_simple.py` only)
- Built FLUPS libraries: `make install_macos` or `make install_static`

## Quick Start

### Single Rank (pytest)

```bash
PYTHONPATH=$PWD/python_api:$PYTHONPATH pytest -v python_api/tests/test_validation_simple.py
```

### Single Rank (standalone)

```bash
PYTHONPATH=$PWD/python_api:$PYTHONPATH python python_api/tests/mpi_validation_simple.py
```

### Multiple Ranks (MPICH/Homebrew)

```bash
# 2 ranks
PYTHONPATH=$PWD/python_api:$PYTHONPATH mpirun -np 2 python python_api/tests/mpi_validation_simple.py

# 4 ranks
PYTHONPATH=$PWD/python_api:$PYTHONPATH mpirun -np 4 python python_api/tests/mpi_validation_simple.py

# 8 ranks
PYTHONPATH=$PWD/python_api:$PYTHONPATH mpirun -np 8 python python_api/tests/mpi_validation_simple.py
```

## Expected Output

```
MPI ranks: 4, nproc grid: [2, 2, 1], relative error: 8.036e-04
OK: validation passed
```

## Notes

- **Memory allocation**: Buffers are allocated using `nmem` (padded memory size) instead of `nloc` to respect FLUPS internal padding.
- **OpenMP warning**: If you see `libomp` duplicate library warnings, set:
  ```bash
  export KMP_DUPLICATE_LIB_OK=TRUE
  ```
- **MPI troubleshooting**: On macOS, ensure your `mpirun` matches the MPI used to compile FLUPS (MPICH from Homebrew recommended).
