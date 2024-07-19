#ifndef BACKEND_FFTW_HPP
#define BACKEND_FFTW_HPP

#include <fftw3.h>

#include <string>

#include "defines.hpp"

using std::string;

// Datatypes
struct flups_fft_arch {
    int dim;
    int howmany;
    int stride;
    int dist;
};

using flups_fft_plan     = fftw_plan;
using flups_shuffle_plan = fftw_plan;
using flups_complex      = fftw_complex;
using flups_iodim        = fftw_iodim;

typedef enum {
    FLUPS_DCT_I   = FFTW_REDFT00,
    FLUPS_DCT_II  = FFTW_REDFT10,
    FLUPS_DCT_III = FFTW_REDFT01,
    FLUPS_DCT_IV  = FFTW_REDFT11,
    FLUPS_DST_I   = FFTW_RODFT00,
    FLUPS_DST_II  = FFTW_RODFT10,
    FLUPS_DST_III = FFTW_RODFT01,
    FLUPS_DST_IV  = FFTW_RODFT11,
} flups_r2r_kind_t;

typedef enum {
    memcpyD2H,
    memcpyH2D,
    memcpyD2D,
    memcpyH2H,
} memcpyKind_t;

typedef enum {
    memTypeHost,
    memTypeDevice,
    memTypeUnregistered,
} memType_t;

// Kernel execution
#define m_device
#define m_constant
#define m_host_device
#define m_backend_lambda [=]

template <typename O>
void for_loop(const O &op, const int n, string name) {
    FLUPS_INFO("Launching kernel: %s", name.c_str());
    for (int i = 0; i < n; i++) {
        op(i);
    }
    FLUPS_INFO("Kernel %s done", name.c_str());
}

#endif  // BACKEND_FFTW_HPP
