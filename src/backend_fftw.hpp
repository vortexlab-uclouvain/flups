#ifndef BACKEND_FFTW_HPP
#define BACKEND_FFTW_HPP

#include <fftw3.h>
#include <omp.h>

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

using flups_kind         = fftw_r2r_kind;
using flups_fft_plan     = fftw_plan;
using flups_shuffle_plan = fftw_plan;
using flups_complex      = fftw_complex;

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
