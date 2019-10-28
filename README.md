# FLUPS - A Fourier-based Library of Unbounded Poisson Solvers

#### Licensing and authorship
Copyright © UCLouvain 2019

FLUPS is distributed under GPLv3.0 (or any later version) license.

The main authors are (by alphabetical order):
- Denis-Gabriel Caprace
- Thomas Gillis

For the list of all the contributors to the development of FLUPS, description and a complete License: see LICENSE file.

#### Citation information
If you use FLUPS, please cite it as follows in your publications:
- Caprace et al., **FLUPS - A Fourier-based Library of Unbounded Poisson Solvers**, SIAM Journal on Scientific Computing, 2019 (under review)


### Installation
#### 1. Dependencies
First, you need to install the dependencies, typically using the following configuration commands (for the intel compilers)
- FFTW (> v3.3.8) in the `fftw_prefix` dir:
```shell
CC=icc CXX=icpc FC=ifort ./configure --prefix=fftw_prefix --enable-mpi --enable-openmp --enable-avx2 --enable-shared
```
- HDF5 (> v1.10) in the `hdf5_prefix` dir:
```shell
CC=mpiicc CXX=mpiicpc FC=mpif90 ./configure --prefix=hdf5_prefix --enable-build-mode=production --enable-parallel
```

#### 2. The Library
You need now to create a architecture/compiler dependent file in `make_arch` to define `CXX`, `CXXFLAGS`, `FFTWDIR` and `HDF5DIR`.
For example:
```makefile
#---------------------------------------------------------
# COMPILERS
#---------------------------------------------------------
# specify the compiler (intel in this case, may aslo be gcc)
CXX = mpiicpc

# set the flag (optimisation or not)
CXXFLAGS := -O3 -g -DNDEBUG -stdc++11

#---------------------------------------------------------
# DEPENDENCES DIRECTORIES
#---------------------------------------------------------
FFTWDIR  := fftw_prefix
HDF5DIR  := hdf5_prefix
FFTW_LIB := ${FFTW_DIR}/lib
FFTW_INC := ${FFTW_DIR}/include
HDF5_LIB := ${HDF5_DIR}/lib
HDF5_INC := ${HDF5_DIR}/include
```
By default, the Makefile is looking for `-lfftw3_openmp -lfftw3` and `-lhdf5`. You can overwrite this by changing the variable `FFTW_LIBNAME` and `HDF5_LIBNAME` in your arch file.

Then you need to reference the created configuration file and the prefix you wish to :
```shell
export ARCH_FILE=make_arch/my_arch_dependent_file
export 
```
You can now 

Finally, go to the main folder and type the compilation command.
- Check the compilation details before doing the installation
```shell
make info
## or
ARCH_FILE=make_arch/my_arch_dependent_file PREFIX=/my/lib/prefix make info
```
- Install the library
```shell
make install
## or
ARCH_FILE=make_arch/my_arch_dependent_file PREFIX=/my/lib/prefix make install
```

#### 3. The documentation
To build the documentation, go to the `doc` subfolder and type `doxygen`.

#### 4. Compilation flags
Here is an exhautstive list of the compilation flags that can be used to change the behavior of the code. To use `MY_FLAG`, simply add `-DMY_FLAG` to the variable `CXXFLAGS` in your `make_arch`.
- `DUMP_DBG`: if specified, the solver will I/O fields using the HDF5 library.
- `COMM_NONBLOCK`: if specified, the code will use the non-blocking communication pattern instead of the all to all version.
- `PERF_VERBOSE`: requires an extensive I/O on the communication pattern used. For performance tuning and debugging purpose only.
- `NDEBUG`: use this flag to remove various checks inside the library
- `PROF`: allow you to use the build-in profiler to have a detailed view of the timing in each part of the solve

:warning: You may also change the memory alignement and the FFTW planner flag in the `flups.h` file.

### How to use a solver?

#### Detailed reference
For a detailed view of the API, have a look at @ref flups.h


#### FLUPS in a nutshell
To use the solver, you first need to create a topology
```cpp
int  axis      = 0;              // aligned along the first dimension
int  nglob[3]  = {64, 128, 64};  // global size of 64x64x64
int  nproc[3]  = {2, 1, 3};      // 6 procs; 2 x 1 x 3
bool isComplex = false;          // real data

// no specific alignement => we put a value of 1
Topology *topo = new Topology(axis, nglob, nproc, isComplex,NULL,1, MPI_COMM_WORLD);

// define additional quantities
double L = {1.0, 2.0, 1.0};
double h = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};
```

Then, you can define a new solver and it's boundary condition
```cpp
// define the solver
const BoundaryType mybc[3][2] = {{UNB, UNB}, {EVEN, ODD}, {UNB, EVEN}};  // BC in X,Y,Z
Solver *      mysolver   = new Solver(topo, mybc, h, L);

// setup the solver
mysolver->set_GreenType(HEJ2);
mysolver->setup(false);
```

To solve a field `rhs` that has been defined on the topology, use
```cpp
mysolver->solve(rhs, rhs, SRHS);
```

Then, destroy the solver
```
delete (mysolver);
```

#### Memory usage

For the recommanded configuration of 128^3 unknowns per processor in full unbounded, we have measured the memory usage of FLUPS on a 2000 cores run:
- the all to all version uses ~530Mb (O.253kB/unknown)
- the non-blocking version uses ~560Mb (O.267kB/unknown)

<!--
(1500/(560/128^3))^(1/3)
For 1.5Go, max 168
14*12
21*8 
7*24-->

**CAUTION**
FLUPS was nerver tested above 1024^3 unknowns per core.

### Implementation details and developers guide
#### C++ use
We use the C++ language in a very limited way, on purpose.
The features used are the object oriented layout and some usefull features of the standard library.

#### Conventions

- Put a ```BEGIN_FUNC;``` at the begining of each function
- how to name an action? ```action_mySuperFunction``` where ```action``` = ```set```, ```get```, ```execute```, ```switch```, ```cmpt```
- how to name a function? ```mySuperFunction```
- how to name an class? ```MyClass```
- how to name an type? ```MyType```

#### Format Guide
We follow the Google formating rules, see https://google.github.io/styleguide/cppguide.html for more details

To configure the auto-formatter in VsCode, search in the settings for `C_Cpp.clang_format_fallbackStyle`.

Set then the value:
```{ BasedOnStyle: Google, ColumnLimit: 0, IndentWidth: 4, AlignConsecutiveAssignments: true, AlignConsecutiveDeclarations: true }```.

Inspired from https://clang.llvm.org/docs/ClangFormatStyleOptions.html (*Configurable Format Style Options* section)

#### Memory layout
In this project we choose to handle the memory in a **Fortran** way of doing iven if we are in C/C++.
So, the memory is aligned as a single row of size `n[0] * n[1] * n[2]`.
The fastest rotating index is set to be `n[0]` then `n[1]` and finally `n[2]`.

We have chosen this way of doing to reuse the 3D code in a 2D framework.
Indeed having the last dimension in the slower rotating index does not penalize the loops writting.

As an example, we here is how we access the memory

```cpp
double* data =(double*) flups_malloc(n[0] * n[1] * n[2] * sizeof(double));

for(int iz=0; iz<n[2]; iz++){
    for(int iy=0; iy<n[1]; iy++){
        for(int ix=0; ix<n[0]; ix++){
            // n[0] is the fastest rotating index
            const int id = iz*n[1]*n[0] + iy * n[0] + ix;

            data[id] = 1.0 ;
        }
    }
}
```