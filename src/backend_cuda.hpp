#ifndef BACKEND_CUDA_HPP
#define BACKEND_CUDA_HPP

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <fftw3.h>

#include <string>

#include "defines.hpp"

using std::string;

// Datatypes
using flups_kind    = fftw_r2r_kind;
using flups_complex = cufftDoubleComplex;

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
    double *mem = NULL;
    // cuFFT plan handle
    cufftHandle fft;
    // parameters for the r2r transform
    int        n;        // number of inputs/output real values
    int        dist;     // distance between two successive fft's
    int        stride;   // data stride
    int        howmany;  // number of fft's to perform
    flups_kind kind;     // kinf of the r2r transform
    bool       is_r2r = false;
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
#define m_cudart_check_errors(res)                                                               \
    do {                                                                                         \
        cudaError_t __err = res;                                                                 \
        FLUPS_CHECK(cudaSuccess == __err, "Fatal CUDA RT error: %s", cudaGetErrorString(__err)); \
    } while (0)

#define m_cublas_check_errors(res)                                                                           \
    do {                                                                                                     \
        cublasStatus_t __err = res;                                                                          \
        FLUPS_CHECK(CUBLAS_STATUS_SUCCESS == __err, "Fatal cuBLAS error: %s", cublasGetStatusString(__err)); \
    } while (0)

#define m_cufft_check_errors(res)                                            \
    do {                                                                     \
        cufftResult_t __err = res;                                           \
        FLUPS_CHECK(CUFFT_SUCCESS == __err, "Fatal cuFFT error: %d", __err); \
    } while (0)

#define m_check_errors(msg)                                                                         \
    do {                                                                                            \
        cudaError_t __err = cudaGetLastError();                                                     \
        FLUPS_CHECK(cudaSuccess == __err, "Fatal error: %s (%s).", msg, cudaGetErrorString(__err)); \
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
void for_loop(const O &op, const int n, string name) {
    FLUPS_INFO("Launching kernel: %s", name.c_str());
    for_loop_kernel<O><<<(n + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(op, n);
    m_check_errors(name.c_str());
    // synchronize (kernel launches are asynchronous)
    m_cudart_check_errors(cudaDeviceSynchronize());
    FLUPS_INFO("Kernel %s done", name.c_str());
}

#endif  // BACKEND_CUDA_HPP
