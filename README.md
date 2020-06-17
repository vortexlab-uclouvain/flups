# FLUPS - A Fourier-based Library of Unbounded Poisson Solvers

#### Licensing and authorship
> FLUPS is distributed under Apache 2.0 license, copyright © UCLouvain 2020.

The main authors are (by alphabetical order):
- Denis-Gabriel Caprace
- Thomas Gillis

For the list of all the contributors to the development of FLUPS, description and a complete License: see LICENSE file.

#### Citation information
FLUPS is described in [this paper](https://arxiv.org/abs/2006.09300). If you use FLUPS, please cite it as follows in your publications:
- Caprace et al., **FLUPS - A Fourier-based Library of Unbounded Poisson Solvers**, SIAM Journal on Scientific Computing, 2019 (under review)


### Why should you use FLUPS?
- You can solve the Poisson on rectangular and uniform distributed 2D/3D grids;
- You can use any boundary conditions, including truly unbounded boundary conditions and semi-unbounded conditions
- You can solve may times the same Poisson problem at low cost using precomputed Green's function and communication patterns;
- You can use threads and/or MPI to fasten the execution;
- You can use the build-in profiler to optimize the execution speed;
- You can use any part of the library on its own, especially the pre-computed communications and the FFTs;
- You can apply filters or do any computation you want while in the Fourier space.

### Installation

FLUPS is a C++ library, with an API in C.
The compilation of FLUPS was tested with Intel compilers and GCC.,

#### Dependencies
First, you need to install the dependencies, typically using the following configuration commands (for the intel compilers)
- FFTW (> v3.3.8) in the `fftw_prefix` dir:
```shell
CC=icc CXX=icpc FC=ifort ./configure --prefix=fftw_prefix --enable-mpi --enable-openmp --enable-avx2 --enable-shared
```
- HDF5 (> v1.10) in the `hdf5_prefix` dir:
```shell
CC=mpiicc CXX=mpiicpc FC=mpif90 ./configure --prefix=hdf5_prefix --enable-build-mode=production --enable-parallel
```

#### Compilation
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
For example:
```makefile
FFTW_LIBNAME := -lfftw3_omp -lfftw3
HDF5_LIBNAME := -lhdf5_openmpi
```

Then you need to reference the created configuration file (using `ARCH_FILE`) and the prefix in you wish to install the library (using `PREFIX`).
You can either `export` the variables or reference them later while calling the Makefile.
If no prefix is given, `make install` uses the current working directory to install the library

Finally, go to the main folder and type the compilation command.
- Check the compilation details before doing the installation\
```shell
export ARCH_FILE=make_arch/my_arch_dependent_file
export PREFIX=/my/lib/prefix
make info
## or
ARCH_FILE=make_arch/my_arch_dependent_file PREFIX=/my/lib/prefix make info
```
- Install the library (to the PREFIX location, or by default in ./lib and ./include )
```shell
make install
## or
ARCH_FILE=make_arch/my_arch_dependent_file PREFIX=/my/lib/prefix make install
```

:warning: you must **install** the library. Indeed, we copy some data required by the solver.
If you wish to keep everything local, simply do not give a prefix and the current directory will be selected.

#### Documentation

The documentation is built using Doxygen.
To build the documentation, go to the `./doc` subfolder and type `doxygen`.

#### Available compilation flags
Here is an exhautstive list of the compilation flags that can be used to change the behavior of the code. To use `MY_FLAG`, simply add `-DMY_FLAG` to the variable `CXXFLAGS` in your `make_arch`.
- `DUMP_DBG`: if specified, the solver will I/O fields using the HDF5 library.
- `COMM_NONBLOCK`: if specified, the code will use the non-blocking communication pattern instead of the all to all version.
- `PERF_VERBOSE`: requires an extensive I/O on the communication pattern used. For performance tuning and debugging purpose only.
- `NDEBUG`: use this flag to bypass various checks inside the library
- `PROF`: allow you to use the build-in profiler to have a detailed view of the timing in each part of the solve. Make sure you have created a folder ```./prof``` next to your executable.
- `REORDER_RANKS`: try to reorder the MPI ranks based on the precomputed communication graph, using call to MPI_Dist_graph. We recommend the use of this feature when the number of processes > 128 and the nodes are allocated exclusive for your application, especially on fully unbounded domains.
- `HAVE_METIS`: in combination with REORDER_RANKS, use METIS instead of MPI_Dist_graph to partition the call graph based on the allocated ressources. You must hence install metis for this functionality.

:warning: You may also change the memory alignement and the FFTW planner flag in the `flups.h` file.

### How to use a solver?

#### Detailed reference
The scientific background of the library is explained in "Caprace et al., **FLUPS - A Fourier-based Library of Unbounded Poisson Solvers**, SIAM Journal on Scientific Computing, 2019 (under review)".

A detailed description of the API is provided in the documentation (@ref flups.h), as well as many implementation details.

#### Memory layout
In this project we choose to handle the memory in a **Fortran** way of doing even if we are in C/C++.
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

flups_free(data);
```

#### FLUPS in a nutshell
To use the solver, you first need to create a topology
```cpp
int  axis      = 0;              // aligned along the first dimension
int  lda       = 1;              // scalar field
int  nglob[3]  = {64, 128, 64};  // global size of 64x64x64
int  nproc[3]  = {2, 1, 3};      // 6 procs; 2 x 1 x 3
bool isComplex = false;          // real data

// no specific alignement => we put a value of 1
FLUPS_Topology *topo = flups_topo_new(axis, lda, nglob, nproc, isComplex, NULL, 1, MPI_COMM_WORLD);

// define additional quantities
double L = {1.0, 2.0, 1.0};
double h = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};
```

Then, you can define a new solver and its boundary condition
```cpp
// define the solver
const FLUPS_BoundaryType mybc[3][2] = {{UNB, UNB}, {EVEN, ODD}, {UNB, EVEN}};  // BC in X,Y,Z
FLUPS_Solver *mysolver = flups_init(topo, mybc, h, L,prof);

// setup the solver
flups_set_greenType(mysolver,typeGreen);
flups_setup(mysolver,false);
```

To solve a field `rhs` that has been defined on the topology, use
```cpp
flups_solve(mysolver,rhs, rhs);
```

Then, destroy the solver and the created topology
```
flups_cleanup(mysolver);
flups_topo_free(topo);
```

#### Advanced usage
Examples of usage of FLUPS in C programs are provided in the `./sample` subfolder.

#### Memory footprint
For the recommanded configuration of 128^3 unknowns per processor in full unbounded, we have measured the memory usage of FLUPS on a 2000 cores run:
- the all to all version uses ~530Mb (O.253kB/unknown)
- the non-blocking version uses ~560Mb (O.267kB/unknown)

<!--
(1500/(560/128^3))^(1/3)
For 1.5Go, max 168
14*12
21*8 
7*24-->

:warning: FLUPS was nerver tested above 1024^3 unknowns per core.

### Implementation details and developers guide
#### C++ use
We use the C++ language in a very limited way, on purpose.
The features used are the object oriented layout and some usefull features of the standard library.

#### Conventions
- Put a ```BEGIN_FUNC;``` at the begining and a ```END_FUNC;``` at the end of each function
- Use ```FLUPS_INFO``` for verbosity (several levels available), ```FLUPS_CHECK``` for assertions and ```FLUPS_ERROR``` for error management
- Use ```flups_malloc``` and ```flups_free``` function to allocate/free memory
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

#### Debugging

FLUPS can be compiled with different levels of verbosity. The following compilation flags are accepted:
- ```-DVERBOSE(=1)``` provides basic output with essential information
- ```-DVERBOSE=2``` generates an output at the beginning and at the end of each function call. If the flag ```-DPROF``` is also defined, the execution of each function call is timed and displayed when exiting the function.
- ```-DVERBOSE=3``` or ```-DVERBOSE=4``` adds even more debugging information