#ifndef BACKEND_HIP_HPP
#define BACKEND_HIP_HPP

#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <rocblas/rocblas.h>

#include <string>

#include "backend_gpu.hpp"
#include "h3lpr/macros.hpp"

using std::string;

// Backend print & assert functions
#define BACKEND_INFO(format, ...)                         \
    ({                                                    \
        m_verb_def("HIP backend", format, ##__VA_ARGS__); \
    })

#define BACKEND_CHECK(format, ...)                          \
    ({                                                      \
        m_assert_def("HIP backend", format, ##__VA_ARGS__); \
    })

// Datatypes
using flups_complex = hipfftDoubleComplex;
typedef int*           opt_int_ptr;
typedef double*        opt_double_ptr;
typedef flups_complex* opt_complex_ptr;
typedef hipStream_t    stream_t;
typedef hipEvent_t     event_t;

struct flups_fft_plan {
    // memory pointer (used only for r2r plans)
    double* mem = NULL;
    // twiddle tables of the r2r pre/post-processing (NULL for the types that need none), see flups_plan_r2r
    double* twiddle = NULL;
    // FFT plan handle (opaque)
    hipfftHandle fft;
    // parameters for the r2r transform
    int              n;        // number of inputs/output real values
    int              dist;     // distance between two successive fft's
    int              stride;   // data stride
    int              howmany;  // number of fft's to perform
    flups_r2r_kind_t kind;     // kinf of the r2r transform
    bool             is_r2r = false;
    // parameters for the equivalent r2c transform
    int dim;      // size of the equivalent r2c fft
    int idist;    // input  distance between two successive fft's
    int odist;    // output distance between two successive fft's
    int memsize;  // size of the needed memory space
    // size in bytes of the FFT library work area needed by fft (auto-allocation is disabled,
    // the Solver provides one shared buffer to every plan through flups_plan_set_workarea)
    size_t worksize = 0;
};

// Memcopy kinds map to HIP's memcpy kinds
typedef enum {
    memcpyD2H = hipMemcpyDeviceToHost,
    memcpyH2D = hipMemcpyHostToDevice,
    memcpyD2D = hipMemcpyDeviceToDevice,
    memcpyH2H = hipMemcpyHostToHost,
} memcpyKind_t;

// Memory type enumeration (maps to HIP memory type constants)
typedef enum {
    memTypeHost         = hipMemoryTypeHost,
    memTypeDevice       = hipMemoryTypeDevice,
    memTypeUnregistered = hipMemoryTypeUnregistered,
} memType_t;

// Create a batched 1D plan without auto-allocated work area (see backend_hip.cpp)
void backend_make_plan(hipfftHandle* fft, size_t* worksize, const int dim, const int idist, const int odist, const hipfftType type, const int howmany);

// Check HIP / rocBLAS / hipFFT API errors
#define m_hiprt_check_errors(res)                                                                         \
    do {                                                                                                  \
        hipError_t flups_err_ = res;                                                                      \
        BACKEND_CHECK(hipSuccess == flups_err_, "Fatal HIP RT error: %s", hipGetErrorString(flups_err_)); \
    } while (0)

#define m_rocblas_check_errors(res)                                                                      \
    do {                                                                                                 \
        rocblas_status flups_err_ = res;                                                                 \
        BACKEND_CHECK(rocblas_status_success == flups_err_, "Fatal rocBLAS error: %d", (int)flups_err_); \
    } while (0)

#define m_hipfft_check_errors(res)                                                              \
    do {                                                                                        \
        hipfftResult flups_err_ = res;                                                          \
        BACKEND_CHECK(HIPFFT_SUCCESS == flups_err_, "Fatal hipFFT error: %d", (int)flups_err_); \
    } while (0)

// Kernel execution
#define BLOCK_SIZE 256
#define m_device __device__
#define m_constant __constant__
#define m_host_device __host__ __device__
#define m_backend_lambda [=] __device__

template <typename O>
__global__ void hip_kernel_(const O op, const int n) {
    // Grid-stride loop
    const int start  = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;
    for (int i = start; i < n; i += stride) {
        op(i);
    }
}

template <typename O>
void parallel_for(const O& op, const int n, string name, stream_t stream) {
    // nothing to do (e.g. a rank without any unknown): a 0-block launch is a HIP error
    if (n <= 0) {
        return;
    }
    BACKEND_INFO("Launching kernel: %s", name.c_str());
    const int grid = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    // Launch kernel using HIP runtime
    hipLaunchKernelGGL((hip_kernel_<O>),
                       dim3(grid),
                       dim3(BLOCK_SIZE),
                       0,
                       stream,
                       op,
                       n);
    // Note: hipLaunchKernelGGL is a macro that both launches and returns immediately (async w.r.t. host)
    BACKEND_INFO("Kernel %s done", name.c_str());
}

/**
 * @brief functor flattening a two-level op(io, ii) into the 1D index of parallel_for (ii is the fast index)
 */
template <typename O>
struct backend_flat2d_ {
    O             op;
    int           n_inner;
    m_device void operator()(const int idx) const { op(idx / n_inner, idx % n_inner); }
};

/**
 * @brief apply op(io, ii) for io in [0, n_outer) and ii in [0, n_inner): same launch as parallel_for over n_outer * n_inner
 */
template <typename O>
void parallel_for_2d(const O& op, const int n_outer, const int n_inner, string name, stream_t stream) {
    backend_flat2d_<O> flat{op, n_inner};
    parallel_for(flat, n_outer * n_inner, name, stream);
}

#endif  // BACKEND_HIP_HPP
