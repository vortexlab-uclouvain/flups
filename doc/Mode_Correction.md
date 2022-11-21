# Correction of the modes 
Both versions of the solver (cell-centred or node-centred) use different FFT types to impose the correct boundary conditions. In some cases, points in memory (in general, the first or the last) are not taken into account by fftw as they are null, however they are required when performing the spectral convolution with the Green's function. Therefore, we postprocess the data after the 1D FFTs to ensure that the Green's function and the RHS are evaluated at the same wavenumber. The following tables detail those corrections.

Here is the list of all the post-processings: 
- `NULL_FIRST_POINT` resets the first point in memory to 0
- `NULL_FIRST_POINT` resets the last point in memory to 0
- `ENFORCE_PERIOD` corrects the last point and enforce the periodicity of the domain

All of this variables are set in the `FFTW_plan_dim` class (either node or cell).

## Cell-centred case 
In the cell-centred case, the RHS is computed in a cell-centred framework while the Green's functions are computed in a node-centred framework as it characterizes the influence of a source located in x = 0 on any other location in the domain. This is the main reasons to apply corrections.


| BCs            | Direction| RHS <br> FFT types |  RHS <br>  Correction | Green's function <br>  Correction |
| ----           | ---------|----     |-                  |-|                   
| PER - PER      | Forward  | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
|                | Backward | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
| ODD - ODD      | Forward  | DST-II - `RODFT10`   | `NULL_FIRST_POINT`| `POSTPRO_NONE`|
|                | Backward | DST-III - `RODFT01`  | `NULL_LAST_POINT`|`POSTPRO_NONE`|
| ODD - EVEN     | Forward  | DST-IV - `RODFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
|                | Backward | DST-IV - `RODFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
| EVEN - ODD     | Forward  | DCT-IV - `REDFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
|                | Backward | DCT-IV - `REDFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
| EVEN - EVEN    | Forward  | DCT-II - `REDFT10`   | `NULL_LAST_POINT`                          |`POSTPRO_NONE`|
|                | Backward | DCT-III - `REDFT01`  | `POSTPRO_NONE`                              | `POSTPRO_NONE`|




## Node-centred case
In the node-centred case, both the RHS and the Green's function are computed in a node-centred framework. However, FFTW discard some points hen dealing with even boundary, as they do not carry any information (the boundary is null). We therefore have to correct the output field of FFTW to have the RHS and GF evaluated at the same wavenumber. 

| BCs            | Direction| RHS <br> FFT types |  RHS <br>  Correction | Green's function <br>  Correction |
| ----           | ---------|----     |-                  |-|                   
| PER - PER      | Forward  | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
|                | Backward | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
| ODD - ODD      | Forward  | DST-I - `RODFT00`    | ` NULL_FIRST_POINT` + `NULL_LAST_POINT`| `POSTPRO_NONE`|
|                | Backward | DST-I - `RODFT00`    | ` NULL_FIRST_POINT` + `NULL_LAST_POINT`                         |`POSTPRO_NONE`|
| ODD - EVEN     | Forward  | DST-III - `RODFT01`  | `NULL_LAST_POINT`| modes shifted by 1/2|
|                | Backward | DST-II  - `RODFT10`  | `NULL_FIRST_POINT` | modes shifted by 1/2|
| EVEN - ODD     | Forward  | DCT-IV - `REDFT11`   | `NULL_LAST_POINT`                          | modes shifted by 1/2|
|                | Backward | DCT-IV - `REDFT11`   | `NULL_LAST_POINT`                              | modes shifted by 1/2|
| EVEN - EVEN    | Forward  | DCT-I - `REDFT00`    | `POSTPRO_NONE`                              |`POSTPRO_NONE`|
|                | Backward | DCT-I - `REDFT00`    | `POSTPRO_NONE`                              | `POSTPRO_NONE`|