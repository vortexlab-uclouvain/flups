/**
 * @file backend_fftw.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef BACKEND_FFTW_HPP
#define BACKEND_FFTW_HPP

#include <fftw3.h>

#include <string>

#include "flups_interface.h"
#include "h3lpr/macros.hpp"

using std::string;

// FFTW planner flag driven by the NDEBUG flag
#ifndef FFTW_FLAG
#ifdef NDEBUG
#define FLUPS_FFTW_FLAG FFTW_PATIENT
#else
#define FLUPS_FFTW_FLAG FFTW_ESTIMATE
#endif
#else
#define FLUPS_FFTW_FLAG FFTW_FLAG
#endif

// Backend print & assert functions
#define BACKEND_INFO(format, ...)                          \
    ({                                                     \
        m_verb_def("FFTW backend", format, ##__VA_ARGS__); \
    })

#define BACKEND_CHECK(format, ...)                           \
    ({                                                       \
        m_assert_def("FFTW backend", format, ##__VA_ARGS__); \
    })

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

typedef int* __restrict __attribute__((aligned(FLUPS_ALIGNMENT)))           opt_int_ptr;
typedef double* __restrict __attribute__((aligned(FLUPS_ALIGNMENT)))        opt_double_ptr;
typedef flups_complex* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_complex_ptr;
typedef int                                                                 stream_t;
typedef int                                                                 event_t;

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

/**
 * @brief apply op(i) for i in [0, n)
 */
template <typename O>
void parallel_for(const O& op, const int n, string name, stream_t) {
    BACKEND_INFO("Launching kernel: %s", name.c_str());
    for (int i = 0; i < n; i++) {
        op(i);
    }
    BACKEND_INFO("Kernel %s done", name.c_str());
}

/**
 * @brief apply op(io, ii) for io in [0, n_outer) and ii in [0, n_inner), ii being the contiguous index
 *
 */
template <typename O>
void parallel_for_2d(const O& op, const int n_outer, const int n_inner, string name, stream_t) {
    BACKEND_INFO("Launching kernel: %s", name.c_str());
    for (int io = 0; io < n_outer; io++) {
        for (int ii = 0; ii < n_inner; ii++) {
            op(io, ii);
        }
    }
    BACKEND_INFO("Kernel %s done", name.c_str());
}

#endif  // BACKEND_FFTW_HPP
