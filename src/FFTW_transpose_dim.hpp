/**
 * @file FFTW_transpose_dim.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef FFTW_TRANSPOSE_DIM_HPP
#define FFTW_TRANSPOSE_DIM_HPP

#include "fftw3.h"
#include "defines.hpp"
#include "mpi.h"

class FFTW_transpose_dim
{
    // we should use the fftw_mpi_plan_many_transpose function...
    // see FFTW p 63
    // see on the size using fftw_mpi_local_size_many_transposed

    int _size[3][3]; /**< @brief size[i] is the size of the data when aligned in the ith direction */
    size_t _pencil_start[3];
    size_t _pencil_end[3];

    MPI_Comm _comm[3];

    FFTW_transpose_dim(int global_size[3],const int nproc);

};

#endif
