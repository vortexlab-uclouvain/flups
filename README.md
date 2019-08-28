# An unbounded-periodic Poisson solver

### How to use a solver?

To use the solver, you first need to create a topology
```cpp
int  axis      = 0;              // aligned along the first dimension
int  nglob[3]  = {64, 128, 64};  // global size of 64x64x64
int  nproc[3]  = {2, 1, 3};      // 6 procs; 2 x 1 x 3
bool isComplex = false;          // real data

Topology *topo = new Topology(axis, nglob, nproc, isComplex);

// define additional quantities
double L = {1.0, 2.0, 1.0};
double h = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};
```

Then, you can define a new solver and it's boundary condition
```cpp
// define the solver
const BoundaryType mybc[3][2] = {{UNB, UNB}, {EVEN, ODD}, {UNB, EVEN}};  // BC in X,Y,Z
FFTW_Solver *      mysolver   = new FFTW_Solver(topo, mybc, h, L);

// setup the solver
mysolver->set_GreenType(HEJ2);
mysolver->setup();
```

To solve a field `rhs` that has been defined on the topology, use
```cpp
mysolver->solve(topo, rhs, rhs, UP_SRHS);
```

Then, destroy the solver
```
delete (mysolver);
```

### Implementation details
#### C++ use
We use the C++ language in a very limited way, on purpose.
The features used are the object oriented layout and some usefull features of the standard library.

#### Conventions

- Put a ```BEGIN_FUNC``` at the begining of each function
- how to name an action? ```action_mySuperFunction``` where ```action``` = ```set```, ```get```, ```execute```, ```switch```, ```cmpt```
- how to name a function? ```mySuperFunction```
- how to name an class? ```MyClass```
- how to name an type? ```MyType```

#### Format Guide
We follow the Google formating rules, see https://google.github.io/styleguide/cppguide.html for more details

To configure the auto-formatter in VsCode, search in the settings for `C_Cpp.clang_format_fallbackStyle`.

Set then the value:
```{ BasedOnStyle: Google, ColumnLimit: 0, IndentWidth: 4, AlignConsecutiveAssignements: true, AlignConsecutiveDeclarations: true }```.

Inspired from https://clang.llvm.org/docs/ClangFormatStyleOptions.html (*Configurable Format Style Options* section)

#### Memory layout
In this project we choose to handle the memory in a **Fortran** way of doing iven if we are in C/C++.
So, the memory is aligned as a single row of size `n[0] * n[1] * n[2]`.
The fastest rotating index is set to be `n[0]` then `n[1]` and finally `n[2]`.

We have chosen this way of doing to reuse the 3D code in a 2D framework.
Indeed having the last dimension in the slower rotating index does not penalize the loops writting.

As an example, we here is how we access the memory

```cpp
double* data =(double*) malloc(n[0] * n[1] * n[2] * sizeof(double));

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

### Authors
By alphabetical order
- Denis-Gabriel Caprace (main author)
- Thomas Gillis (main author)


### Citation
Please cite the following paper