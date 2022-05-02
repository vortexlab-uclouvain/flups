# Testing in FLUPS

Tests performed on a cluster in Flups use the Google Test framework, which organizes individual tests into a series of test suites and returns a xml reports which can be post processed by the Gitlab platform. 

## Test description 
The tests verify the spatial convergence of the solver using all the available kernels with various combination of boundary conditions. For each test, we solve the Poisson equation 
<p align="center"><img src="https://latex.codecogs.com/gif.latex?\nabla^2 \phi = f, " title="Poisson equation" /></p>
 and the convergence of the solver is measured using the infinite norm of the error, defined as 
 <p align="center"><img src="https://latex.codecogs.com/gif.latex?  E_\infty = \sup_{x,y,z} \{|\phi(x,y,z) - \phi_{ref}(x,y,z)|\} " title="Error inf" /></p>
where <p align="center"><img src="https://latex.codecogs.com/gif.latex?  \phi_{ref}(x,y,z) = a(x)b(y)c(z) " title="Error inf" /></p> is an analytical solution composed as a product of 1D functions. Those functions, a, b, and c,  are defined depending on each set of boundary conditions. Given these analytical expression, the right and side f is computed as 
<p align="center"><img src="https://latex.codecogs.com/gif.latex?\nabla^2 f(x,y,z) = \frac{d^2}{dx^2}a(x)b(y)c(z) + a(x)\frac{d^2}{dy^2}b(y)c(z) + a(x)b(y)\frac{d^2}{dz^2}c(z), " title="Poisson equation" /></p>

The analytical solutions are detailed in the FLUPS paper. 

## Implementation notes 
The heart of the tests is composed by two classes: the `AnalyticalField` and the `ConvergenceTest`. 

  1. The `AnalyticalField` stores the analytical expression for the RHS and the solution based on a set of boundary conditions. In other words, given the bcs in a certain direction, it stores and provides the 1D analytical solution (a(x), b(y) and c(z) in the previous paragraph) and their second derivatives. 

  2. The `ConvergenceTest` is the main class of those tests. It inherits from the `testing::TestWithParam<int>` class of the Google library and implements all the function needed to perform a convergence test: the creation of the domain and of the solver, the computation of the RHS (based on an `AnalyticalField`), the solve itself, the computation of the error and of the convergence order. The templated parameters are: the location of the data (cell-centred (test number 0 to 5999) or node-centred (test number 6000 to 11 999), a specific kernel and the combination of the boundary conditions. The `ConvergenceTest` performs the tests using 6 different kernels. If one of the kernel does not have the right convergence order, the test fails. 


## Compilation 
Those tests depends on the Google test library. You should install it using :
```shell
 git clone https://github.com/google/googletest.git
 cd /googletest
 cmake . -DCMAKE_INSTALL_PREFIX=/soft/googletest
 make install -j
 rm -rf /googletest
 ```
Then, you should add the location of the library in your ARCH_FILE:
```shell 
## Specify the location of Googletest library (by default, assumed to be /usr/lib/)
GTEST_INC := /soft/googletest/include
GTEST_LIB := /soft/googletest/lib
```

You can launch the test using 
```shell
mpirun -n 8 ./flups_test_a2a 
```

You can select a certain test using 
```shell 
mpirun -n 8 ./flups_test_a2a --gtest_filter=FlupsValidation/ConvergenceTest.AllBoundaryConditions/0-11999
```

All the tests can be listed using 
```shell 
./flups_test_a2a --gtest_list_tests
```