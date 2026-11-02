/**
 * @file backend_gpu.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef BACKEND_GPU_HPP
#define BACKEND_GPU_HPP

typedef enum {
    FLUPS_DCT_I,
    FLUPS_DCT_II,
    FLUPS_DCT_III,
    FLUPS_DCT_IV,
    FLUPS_DST_I,
    FLUPS_DST_II,
    FLUPS_DST_III,
    FLUPS_DST_IV,
} flups_r2r_kind_t;

// Helper struct, see FFTW doc.
struct flups_iodim {
    int n;
    int is;
    int os;
};

struct flups_fft_arch {
    int dim;
    int howmany;
    int stride;
    int dist;
};

struct flups_shuffle_plan {
    int n0;
    int n1;
};

#endif  // BACKEND_GPU_HPP
