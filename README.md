# FLUPS - A Fourier-based Library of Unbounded Poisson Solvers

## Licensing and authorship
> FLUPS is distributed under BSD-3 clause license, copyright (c) UCLouvain 2022.

The main authors are (by alphabetical order):
- Pierre Balty (`v2.x`)
- Denis-Gabriel Caprace (`v1.x`)
- Thomas Gillis (`v1.x` and `v2.x`)

For the list of all the contributors to the development of FLUPS, description and a complete License: see `LICENSE` file.


### Citation information

FLUPS' design, implementation, and performances are described in two papers.

If you use FLUPS, please cite them in your publications:
- [Balty et al.](https://arxiv.org/abs/2211.07777), **FLUPS - a flexible and performant massively parallel Fourier transform library**, submitted, 2022
- [Caprace et al.](https://arxiv.org/abs/2006.09300), **FLUPS - A Fourier-based Library of Unbounded Poisson Solvers**, SIAM Journal on Scientific Computing, 2021


## Why should you use FLUPS?
- You can solve the Poisson on rectangular and uniform distributed grids;
- You can use either cell-centred or node-centred data layout;
- You can use any boundary conditions, including truly unbounded boundary and semi-unbounded conditions;
- You can solve many times the same Poisson problem at low cost using precomputed Green's function and communication patterns;
- You can use MPI to fasten the execution;
- You can use the profiler of `h3lpr` (see dependency) to optimize the execution speed;
- You can use any part of the library on its own, especially the pre-computed communications and the FFTs;
- You can apply filters or do any computation you want while in the Fourier space.

## Installation

FLUPS is a C++ library, with an API in C.
The compilation of FLUPS was tested with GCC (v9.4) and clang (v12.0).

### Dependencies
First, you need to install the dependencies, typically using the following configuration commands (for the mpich compilers)

- H3LPR in the `h3lpr_prefix` dir:
```shell
git clone git@github.com:vanreeslab/h3lpr.git
cd h3lpr 
ARCH_FILE=... make install -j 
```

- FFTW compatible implementation (e.g. `fftw3` > v3.3.8) in the `fftw_prefix` dir:
```shell
CC=mpicc CXX=mpic++ ./configure --prefix=fftw_prefix --enable-mpi --enable-openmp --disable-fortran --enable-shared
```

- For debugging purpose - HDF5 (> v1.10) in the `hdf5_prefix` dir:
```shell
CC=mpicc CXX=mpic++ FC=mpif90 ./configure --prefix=hdf5_prefix --enable-optimization=high --enable-build-mode=production
```


### Compilation
Then, you need to create a architecture/compiler dependent file in `make_arch` to define `CXX`, `CXXFLAGS`, `H3LPR_DIR`, `FFTW_DIR` and `HDF5_DIR`.
For example:
```makefile
#---------------------------------------------------------
# COMPILERS
#---------------------------------------------------------
# specify the compiler (intel in this case, may aslo be gcc)
CXX = mpicc

# set the flag (optimisation or not)
CXXFLAGS := -O3 -g -DNDEBUG -stdc++11
LDFLAGS := -fopenmp 

#---------------------------------------------------------
# DEPENDENCES DIRECTORIES
#---------------------------------------------------------
H3LPR_DIR := h3lpr_prefix
H3LPR_LIB := ${H3LPR_DIR}/lib
H3LPR_INC := ${H3LPR_DIR}/include

FFTW_DIR := fftw_prefix
FFTW_LIB := ${FFTW_DIR}/lib
FFTW_INC := ${FFTW_DIR}/include

# If needed
HDF5_DIR  := hdf5_prefix
HDF5_LIB := ${HDF5_DIR}/lib
HDF5_INC := ${HDF5_DIR}/include
```

By default, the Makefile is looking for `-lh3lpr`, `-lfftw3_openmp -lfftw3` and `-lhdf5`. You can overwrite this by changing the variable `H3LPR_LIBNAME`, `FFTW_LIBNAME` and `HDF5_LIBNAME` in your arch file.
For example:
```makefile
H3LPR_LIBNAME := -lh3lpr
FFTW_LIBNAME := -lfftw3_omp -lfftw3
HDF5_LIBNAME := -lhdf5_openmpi
```

Then you need to reference the created configuration file (using `ARCH_FILE`) and the prefix in you wish to install the library (using `PREFIX`).
You can either `export` the variables or reference them later while calling the Makefile.
If no prefix is given, `make install` uses the current working directory to install the library

Finally, go to the main folder and type the compilation command.
- Verify the compilation details before doing the installation
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

/!\ you must **install** the library as it copies some data required by the solver.
If you wish to keep everything local, simply do not give a prefix and the current directory will be selected.

**Performance notes**:  To increase the performance of the code, we highly recommend compiling it with _Link Time Optimisation_ (LTO). To do that, add the flag ` -flto` to your `CXXFLAGS` and `LDFLAGS` variables. In addition, you must ensure that your utility tool to create the library archive can build an archive file that libLTO can use at link time. Finally, if you have an architecture which supports LTO, overwrite the `AR` variable in your `make_arch`. 
```shell
AR := gcc-ar 
```
By default, the Makefile use the open-source utility tool `ar`.

### Available compilation flags
Here is an exhautstive list of the compilation flags that can be used to change the behavior of the code. To use `MY_FLAG`, simply add `-DMY_FLAG` to the variable `OPTS` in your `make_arch`.
- `HAVE_HDF5` : Enable the use of function to dump flups fields. When using this flag, you should detail your `HDF5` lib and include in your `make_arch`
- `COMM_NONBLOCK`: if specified, the code will use the non-blocking communication pattern instead of the all-to-all version.
- `PERF_VERBOSE`: requires an extensive I/O on the communication pattern used. For performance tuning and debugging purpose only.
- `NDEBUG`: use this flag to bypass various checks inside the library
- `PROF`: allow you to use the build-in profiler to have a detailed view of the timing in each part of the solve. Make sure you have created a folder `./prof` next to your executable.
- `REORDER_RANKS` (deprecated): try to reorder the MPI ranks based on the precomputed communication graph, using call to MPI_Dist_graph. We recommend the use of this feature when the number of processes > 128 and the nodes are allocated exclusive for your application, especially on fully unbounded domains.
- `HAVE_METIS` (deprecated): in combination with REORDER_RANKS, use METIS instead of MPI_Dist_graph to partition the call graph based on the allocated ressources. You must hence install metis for this functionality. This part of the code has never been demonstrated to show a real increase of performances and therefore is depracted. However we still conserve the code active with this flag.
- `COMM_DPREC`: will use the deprectated communication implementation (slower initalization time, kept for comparison purposes)
- `BALANCE_DPREC`: will use the deprecated distribution of unknowns on the ranks
- `MPI_40` : Use this flag to apply some fancy parameters to allow faster MPI calls if you have a MPI-4.0 compliant version
- `FFTW_FLAG` drives the flag used to init the fftw routines and can be set to ` FFTW_ESTIMATE`, ` FFTW_MEASURE`, ` FFTW_PATIENT`, or `FFTW_EXHAUSTIVE`.
- `MPI_NO_ALLOC` Use this flag to use the system allocation functions instead of the MPI ones when allocating data. 
- `MPI_BATCH_SEND=x` will have `x` non-blocking active send request, set to `INT_MAX` to send them all at once.
- `HAVE_WISDOM=\"path/to/filename\"` indicates that FFTW wisdom can be found at the given filename.


/!\ You may also change the memory alignement and the FFTW planner flag in the `flups.h` file.

### Documentation

The documentation is built using Doxygen.
To build the documentation, go to the `./doc` subfolder and type `doxygen`.
## How to use a solver?
### Detailed reference
The scientific background of the library is explained in Caprace et al., **FLUPS - A Fourier-based Library of Unbounded Poisson Solvers**, SIAM Journal on Scientific Computing, 2019 and in Balty et al., **FLUPS - a flexible and performant massively parallel Fourier transform library**, submitted 2022.

FLUPS solves two types of equations:
- laplacian(phi) = rhs, with phi and rhs either scalars or vectors 
- laplacian(phi) = rot(rhs), with phi and rhs vectors (also code Biot-Savart mode)

A detailed description of the API is provided (in the documentation)[doc/documentation.html] (@ref flups.h), as well as many implementation details.

### Memory layout
In this project we choose to handle the memory in a **Fortran** way of doing even if we are in C/C++.
So, the memory is aligned as a single row of size `n[0] * n[1] * n[2]`.
The fastest rotating index is set to be `n[0]` then `n[1]` and finally `n[2]`.

We have chosen this way of doing to reuse the 3D code in a 2D framework.
Indeed having the last dimension in the slower rotating index does not penalize the loops writting.

As an example, we here is how we access the memory for a scalar field:

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

Vector components are treated using a leading index of arrays (slowest rotating index), and thus corresponds to an additional outer loop.


### FLUPS in a nutshell
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

Then, you can define a new solver, its data-layout and its boundary condition
```cpp
// define the solver
FLUPS_BoundaryType* mybc[3][2];
for(int id=0; id<3; id++){
    for(int is=0; is<2; is++){
        mybc[id][is] = (FLUPS_BoundaryType*) flups_malloc(sizeof(int)*lda);
        for(int ida = 0; ida < lda; ida++) mybc[id][is][ida] = EVEN; 
    }
}
FLUPS_CenterType center_type[3] = {CELL_CENTER, CELL_CENTER, CELL_CENTER};
FLUPS_Solver *mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, prof);

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
for (int id = 0; id < 3; id++) {
    for (int is = 0; is < 2; is++) {
        flups_free(mybc[id][is]);
    }
}
```

### Code examples
Examples of usage of FLUPS in C programs are provided in the `./sample` subfolder.
This includes: 
* `validation`: the exe used for validation and scalability analysis (see our reference publication). This also constitutes an example of how to use FLUPS within a C++ client code, for the scalar Poisson equation.
* `solve_vtube`: another validation test case on a 2-D vortex tube. It may be used as an example on how to use FLUPS to solve the vector Poisson equation and the Biot-Savart mode.
* `solve_advanced_C`: an example showing how to embed flups in a C code, also showing how to use some advanced features (e.g. performing 3-D FFTs separately).


### Make the most of the parallel implementation

FLUPS features hybrid distributed (maintained)/shared(deprecated version) memory capabilities, enabling the library to adapt to a variety of software/hardware configurations. Also, two types of communications schemes are available: all-to-all and non-blocking. The user can select one option or the other at compilation time, through the `COMM_NONBLOCK` flag. Among the two non-blocking implementations, the user can choose to use _persistent_ communication or communication based on _MPI\_Datatype_.

The actual performance of the library (in terms of time-to-solution) depends a.o. on the number of unknowns per CPU, on the type of boundary conditions and on the architectures it runs on.  We here provide some guidelines for the user to determine the optimal setup (see reference publication for more details):
- We highly recommend the use of distributed memory when possible, even if FLUPS can run in a pure OpenMP mode.
- The all-to-all implementation should be considered as the default robust option. However, acceleration is possible using the non-blocking version, in particular when:
    - the number of unknowns per core is high (~128^3)
    - the total number of core is not too high (~< 10k)
- (deprecated) The mixed use of OpenMP and MPI is supported, and should only be considered in combination with the non-blocking implementation. However, the related performance is highly dependent on the computer architecture. 
- (deprecated) Should you use shared memory (`OMP_NUM_THREADS>1`), each thread must be handled by a distinct core (no hyper threading). Computer nodes providing non-uniform memory accesses 

We encourage the user seeking for optimal performance to run short dedicated tests on the targeted architecture. The `validation` executable, when compiled with the `PROF` option, can be used to time the execution. A basic comparison of performance on a typical-size problem should involve at least:
1. the all-to-all implementation without thread
2. the non-blocking implementation without thread
3. the isr implementation without thread


### Memory footprint
For the recommanded configuration of 128^3 unknowns per processor in full unbounded, we have measured the memory usage of FLUPS-v1.0 on a 2000 cores run:
- the all-to-all version uses ~530Mb (O.253kB/unknown)
- the non-blocking version uses ~560Mb (O.267kB/unknown)

<!--
(1500/(560/128^3))^(1/3)
For 1.5Go, max 168
14*12
21*8 
7*24-->

<!--
:warning: FLUPS was nerver tested above 1024^3 unknowns per core.
-->



## Implementation details and developers guide
### C++ use
We use the C++ language in a very limited way, on purpose.
The features used are the object oriented layout and some usefull features of the standard library.

### Conventions
- Put a `BEGIN_FUNC;` at the begining and a `END_FUNC;` at the end of each function
- Use `FLUPS_INFO` for verbosity (several levels available), `FLUPS_CHECK` for assertions and `FLUPS_ERROR` for error management
- Use `flups_malloc` and `flups_free` function to allocate/free memory
- how to name an action? `action_mySuperFunction` where `action` = `set`, `get`, `execute`, `switch`, `cmpt`
- how to name a function? `mySuperFunction`
- how to name an class? `MyClass`
- how to name an type? `MyType`

### Format Guide
We follow the Google formating rules, see https://google.github.io/styleguide/cppguide.html for more details

To configure the auto-formatter in VsCode, search in the settings for `C_Cpp.clang_format_fallbackStyle`.

Set then the value:
`{ BasedOnStyle: Google, ColumnLimit: 0, IndentWidth: 4, AlignConsecutiveAssignments: true, AlignConsecutiveDeclarations: true }`.

Inspired from https://clang.llvm.org/docs/ClangFormatStyleOptions.html (*Configurable Format Style Options* section)

### Debugging

FLUPS can be compiled with different levels of verbosity. The following compilation flags are accepted:
- `-DVERBOSE(=1)` provides basic output with essential information
- `-DVERBOSE=2` generates an output at the beginning and at the end of each function call. If the flag `-DPROF` is also defined, the execution of each function call is timed and displayed when exiting the function.
- `-DVERBOSE=3` or `-DVERBOSE=4` adds even more debugging information


### Testing 

The continuous integration of FLUPS is based on the tools provided by Gitlab. Different types of tests are performed depending on the situation:
- Any `push` event on any branches will trigger the _build test_. FLUPS is compiled with different compilation flags and coupled with various test cases (written in c++ or c). If there is a problem during the compilation, the test fails. 

- Any `merge request` triggers some _validation tests_. We test all the possible combination of boundary conditions (1000 possibilities), kernels (8 kernels) and data location (node-centred or cell-centred) using the [Google test library](https://github.com/google/googletest). Basically, we test the spatial convergence of all the kernels with all the combination of boundary conditions. The source code cand be found in the `test` directory while details and explantion of the test can be found [here](test/Readme.md). However, this extremly large amount of tests (16 000 in total)  require extensive computationnal resources. We hence rely on _daily_ testing routine, that uses the `sample/validation/` source code. <br/>
The _daily test_ is a smaller, in-house, test suite, that can be executed on a desktop machine. Diverse boundary conditions, domain size, resolution, procs repartition and kernels are tested and the results are compared to a dataset that has been generated with a validated version of the code. 


### Other resources and information

- [CI testing](doc/CI.md)
- [real to real FFT](doc/Mode_Correction.md)
