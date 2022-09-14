# Correction of the modes 
In both versions of the solver (cell-centred or node-centred), the spectral representation of the right-hand side and the Green's function are not obtained through the same mechanism. The spectral representation of the RHS is computed thanks to Fourier transforms, while the spectral representation of the Green's function is either computed analyticaly or calculated thanks to Fourier transforms depending on the boundary conditions. 

To increase the solver performances and reduce the memory footprint, we rely on Discrete Sine Transform (DST) and Discrete Cosine Transform (DCT), when dealing with symmetric BCs. The transform type depends on the boundary conditions and the location of the data (cell or node). However, the number of points taken by FFTW and the wavenumber computed may vary from one transform to the other. Therefore, we need to apply some corrections to ensure that the Green's function and the RHS are evaluated at the same wavenumber when performing the spectral multiplication. The following tables detail those corrections.

Here is the list of all the corrections: 
- `ZEROMODE` resets the mode 0 to 0
- `FLIPFLOP` resets the flipflop mode to 0
- `SHIFT_RIGHT` shifts the modes by one index: i -> i+1, always done with `ZEROMODE`
- `SHIFT_LEFT`shifts the modes by one index: i -> i-1, always done with `FLIPFLOP`

All of this variables are set in the `FFTW_plan_dim` class (either node or cell).

## Cell-centred case 
In the cell-centred case, the RHS is computed in a cell-centred framework while the Green's functions are computed in a node-centred framework as it characterizes the influence of a source located in x = 0 on any other location in the domain. This is the main reasons to apply corrections.


| BCs            | Direction| RHS <br> FFT types |  RHS <br>  Correction | Green's function <br>  Correction |
| ----           | ---------|----     |-                  |-|                   
| PER - PER      | Forward  | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
|                | Backward | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
| ODD - ODD      | Forward  | DST-II - `RODFT10`   | `CORRECTION_SHIFTRIGHT` + `CORRECTION_ZEROMODE`| `POSTPRO_NONE`|
|                | Backward | DST-III - `RODFT01`  | `CORRECTION_SHIFTRIGHT` + `CORRECTION_ZEROMODE`|`POSTPRO_NONE`|
| ODD - EVEN     | Forward  | DST-IV - `RODFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
|                | Backward | DST-IV - `RODFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
| EVEN - ODD     | Forward  | DCT-IV - `REDFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
|                | Backward | DCT-IV - `REDFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
| EVEN - EVEN    | Forward  | DCT-II - `REDFT10`   | `CORRECTION_FLIPFLOP`                          |`POSTPRO_NONE`|
|                | Backward | DCT-III - `REDFT01`  | `POSTPRO_NONE`                              | `POSTPRO_NONE`|




## Node-centred case
In the node-centred case, both the RHS and the Green's function are computed in a node-centred framework. However, FFTW discard some points hen dealing with even boundary, as they do not carry any information (the boundary is null). We therefore have to correct the output field of FFTW to have the RHS and GF evaluated at the same wavenumber. 

| BCs            | Direction| RHS <br> FFT types |  RHS <br>  Correction | Green's function <br>  Correction |
| ----           | ---------|----     |-                  |-|                   
| PER - PER      | Forward  | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
|                | Backward | DFT     | `POSTPRO_NONE` |`POSTPRO_NONE`|
| ODD - ODD      | Forward  | DST-I - `RODFT00`    | `CORRECTION_ZEROMODE` + `CORRECTION_FLIPFLOP`| `POSTPRO_NONE`|
|                | Backward | DST-I - `RODFT00`    | `POSTPRO_NONE`                            |`POSTPRO_NONE`|
| ODD - EVEN     | Forward  | DST-III - `RODFT01`  | `CORRECTION_SHIFTLEFT` + `CORRECTION_FLIPFLOP`| modes shifted by 1/2|
|                | Backward | DST-II  - `RODFT10`  | `CORRECTION_SHIFTRIGHT` + `CORRECTION_ZEROMODE`| modes shifted by 1/2|
| EVEN - ODD     | Forward  | DCT-IV - `REDFT11`   | `CORRECTION_FLIPFLOP`                          | modes shifted by 1/2|
|                | Backward | DCT-IV - `REDFT11`   | `POSTPRO_NONE`                              | modes shifted by 1/2|
| EVEN - EVEN    | Forward  | DCT-I - `REDFT00`    | `POSTPRO_NONE`                              |`POSTPRO_NONE`|
|                | Backward | DCT-I - `REDFT00`    | `POSTPRO_NONE`                              | `POSTPRO_NONE`|