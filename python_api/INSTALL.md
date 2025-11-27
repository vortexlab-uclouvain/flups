# PyFLUPS Installation and Testing Guide

## Installation

### 1. Prerequisites

First, ensure you have the FLUPS C++ library compiled:

```bash
cd /path/to/flups-dev-APIpython
make ARCH_FILE=make_arch/make.default
```

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
cd python_api

# Install dependencies
pip install numpy mpi4py

# Install PyFLUPS in development mode
pip install -e .
```

### 4. Set library path

Add to your `.bashrc` or `.zshrc`:

```bash
export FLUPS_LIB_PATH=/path/to/flups-dev-APIpython/lib
```

Or set it temporarily:
```bash
export FLUPS_LIB_PATH=/Users/pbalty/Documents/These/codes/flups-dev-APIpython/lib
```

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

**Error:** `Could not find FLUPS library`

**Solution:** 
- Verify library exists: `ls $FLUPS_LIB_PATH`
- Check for `libflups.so` (Linux) or `libflups.dylib` (macOS)
- Ensure FLUPS is compiled: `make ARCH_FILE=make_arch/make.osx_gcc` (for macOS)

### MPI errors

**Error:** `mpi4py` import fails

**Solution:**
```bash
# Install MPI first (macOS)
brew install openmpi

# Then install mpi4py
pip install mpi4py
```

### Import errors

**Error:** Module not found

**Solution:**
```bash
# Make sure you're in the right directory
cd python_api

# Reinstall
pip install -e .
```

## Development Setup

For development work:

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Run tests (if available)
pytest

# Check code style
flake8 pyflups/
```

## Platform-Specific Notes

### macOS

Use the appropriate make file:
```bash
make ARCH_FILE=make_arch/make.osx_gcc
```

Install MPI via Homebrew:
```bash
brew install open-mpi
```

### Linux

Ensure you have MPI installed:
```bash
# Ubuntu/Debian
sudo apt-get install libopenmpi-dev openmpi-bin

# Fedora/RHEL
sudo dnf install openmpi openmpi-devel
```

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
