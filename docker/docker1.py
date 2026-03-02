import hpccm
from hpccm.building_blocks import apt_get, python, gnu, llvm, openmpi, hdf5, fftw, ofed
from hpccm.primitives import baseimage

# Two stages: build and runtime to reduce final image size
# Build stage
Stage0 = hpccm.Stage()
stage0_image = 'ubuntu:24.04'
Stage0 += baseimage(image=stage0_image)

# ===========================
# Install system packages
apt = apt_get(ospackages=['git', 'wget', 'apt-transport-https', 'ca-certificates',
                          'autotools-dev', 'autoconf', 'libtool', 'automake',
                          'vim',
                          'make', 'cmake',
                          'valgrind',
                          'm4',
                          'libucx-dev',
                          'libstdc++-14-dev'])

Stage0 += apt

# ===========================
# GCC/GFortran 14 - Fortran only
fc = gnu(version='14', cc=False, cxx=False, fortran=True, environment=True, toolchain=False)
Stage0 += fc


# ===========================
# Clang 18 - Main C/C++ toolchain (stable for Ubuntu 24.04)
cc = llvm(version='18', openmp=True, extra_tools=True, toolset=True, environment=True, toolchain=True)
Stage0 += cc


# ===========================
# OpenMPI - use custom Clang toolchain with libc++
Stage0 += ofed(toolchain=cc.toolchain)
omp = openmpi(version='5.0.9', cuda=False, ucx=True, toolchain=cc.toolchain,
              configure_opts=['--disable-mpi-fortran'])
Stage0 += omp

# ===========================
# HDF5
Stage0 += hdf5(version='1.14.6', toolchain=omp.toolchain,
               configure_opts=['--disable-fortran','--enable-parallel','--enable-optimization=high',
                               '--enable-build-mode=production','--with-default-api-version=v110'])

# ===========================
# FFTW - with --enable-shared for .so compatibility
Stage0 += fftw(version='3.3.10', toolchain=omp.toolchain,
               configure_opts=['--disable-fortran','--enable-openmp','--enable-shared','--enable-static'])

# ===========================
# Python 3 only (install directly with apt)
Stage0 += apt_get(ospackages=['python3', 'python3-pip', 'python3-dev', 'python3-venv'])
# Create python command as symlink to python3
Stage0 += hpccm.primitives.shell(commands=[
    'ln -s /usr/bin/python3 /usr/bin/python',
    'pip3 install numpy --break-system-packages',
    'MPICC=/usr/local/openmpi/bin/mpicc pip3 install mpi4py --break-system-packages --no-cache-dir'
])

# ===========================
# Runtime stage - optimized final image
# Use .runtime() to copy only necessary libraries
# without compilation tools
Stage1 = hpccm.Stage()
Stage1 += baseimage(image=stage0_image)


# System packages
Stage1 += apt
Stage1 += cc
Stage1 += fc
Stage1 += omp

# OpenMPI, HDF5, FFTW, Python runtime copied from Stage0
Stage1 += hdf5(version='1.14.6').runtime(_from='0')
Stage1 += fftw(version='3.3.10').runtime(_from='0')

# Python 3 (install directly with apt)
Stage1 += apt_get(ospackages=['python3', 'python3-pip', 'python3-dev'])
# Remove PEP 668 restriction to allow pip install without --break-system-packages
Stage1 += hpccm.primitives.shell(commands=[
    'rm /usr/lib/python3.*/EXTERNALLY-MANAGED',
    'ln -s /usr/bin/python3 /usr/bin/python || true',
    'pip3 install numpy --break-system-packages',
    'MPICC=/usr/local/openmpi/bin/mpicc pip3 install mpi4py --break-system-packages --force-reinstall --no-cache-dir'
])


print(Stage0)
print(Stage1)
