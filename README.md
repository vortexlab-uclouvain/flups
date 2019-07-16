# An unbounded-periodic Poisson solver

## Example

## Implementation details
### C++ use
We use the C++ language in a very limited way, on purpose.
The features used are the object oriented layout and some usefull features of the standard library.

### Memory layout

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