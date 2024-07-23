#ifndef BACKEND_CUDA_HPP
#define BACKEND_CUDA_HPP

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cufft.h>

#include <string>

#include "flups_interface.h"
#include "h3lpr/macros.hpp"

using std::string;

// Backend print & assert functions
#define BACKEND_INFO(format, ...)                          \
    ({                                                     \
        m_verb_def("CUDA backend", format, ##__VA_ARGS__); \
    })

#define BACKEND_CHECK(format, ...)                           \
    ({                                                       \
        m_assert_def("CUDA backend", format, ##__VA_ARGS__); \
    })

// Datatypes
using flups_complex = cufftDoubleComplex;
typedef int*           opt_int_ptr;
typedef double*        opt_double_ptr;
typedef flups_complex* opt_complex_ptr;

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

struct flups_fft_plan {
    // memory pointer (used only for r2r plans)
    double* mem = NULL;
    // cuFFT plan handle
    cufftHandle fft;
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
};

typedef enum {
    memcpyD2H = cudaMemcpyDeviceToHost,
    memcpyH2D = cudaMemcpyHostToDevice,
    memcpyD2D = cudaMemcpyDeviceToDevice,
    memcpyH2H = cudaMemcpyHostToHost,
} memcpyKind_t;

typedef enum {
    memTypeHost         = cudaMemoryTypeHost,
    memTypeDevice       = cudaMemoryTypeDevice,
    memTypeUnregistered = cudaMemoryTypeUnregistered,
} memType_t;

// Check CUDA API errors
#ifndef DNDEBUG
#define m_cudart_check_errors(res)                                                                 \
    do {                                                                                           \
        cudaError_t __err = res;                                                                   \
        BACKEND_CHECK(cudaSuccess == __err, "Fatal CUDA RT error: %s", cudaGetErrorString(__err)); \
    } while (0)

#define m_cublas_check_errors(res)                                                                             \
    do {                                                                                                       \
        cublasStatus_t __err = res;                                                                            \
        BACKEND_CHECK(CUBLAS_STATUS_SUCCESS == __err, "Fatal cuBLAS error: %s", cublasGetStatusString(__err)); \
    } while (0)

#define m_cufft_check_errors(res)                                              \
    do {                                                                       \
        cufftResult_t __err = res;                                             \
        BACKEND_CHECK(CUFFT_SUCCESS == __err, "Fatal cuFFT error: %d", __err); \
    } while (0)

#define m_check_errors(msg)                                                                           \
    do {                                                                                              \
        cudaError_t __err = cudaGetLastError();                                                       \
        BACKEND_CHECK(cudaSuccess == __err, "Fatal error: %s (%s).", msg, cudaGetErrorString(__err)); \
    } while (0)
#else
#define m_cudart_check_errors(msg) \
    { ((void)0); }
#define m_cublas_check_errors(msg) \
    { ((void)0); }
#define m_cufft_check_errors(msg) \
    { ((void)0); }
#define m_check_errors(msg) \
    { ((void)0); }
#endif

// Kernel execution
#define BLOCK_SIZE 256
#define m_device __device__
#define m_constant __constant__
#define m_host_device __host__ __device__
#define m_backend_lambda [=] __device__

template <typename O>
__global__ void for_loop_kernel(const O op, const int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        op(i);
    }
}

template <typename O>
void for_loop(const O& op, const int n, string name) {
    BACKEND_INFO("Launching kernel: %s", name.c_str());
    for_loop_kernel<O><<<(n + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(op, n);
    m_check_errors(name.c_str());
    // synchronize (kernel launches are asynchronous)
    m_cudart_check_errors(cudaDeviceSynchronize());
    BACKEND_INFO("Kernel %s done", name.c_str());
}

#endif  // BACKEND_CUDA_HPP
