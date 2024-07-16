#ifndef BACKEND_BACKEND_CUDA_HPP_
#define BACKEND_BACKEND_CUDA_HPP_

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <fftw3.h>

#include <string>

#include "defines.hpp"

using std::string;

// cuBLAS handle
extern cublasHandle_t cuHandle;

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

struct flups_plan {
    // memory pointer (used only for r2r plans)
    double *mem = NULL;
    // cuFFT plan handle
    cufftHandle cufft;
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

// Check CUDA RT API errors
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
// TODO: set default block size via compile param ?
#define BLOCK_SIZE 256
#define m_kernel [=] __device__

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

// Init & cleanup backend
inline void init_backend() {
    m_cublas_check_errors(cublasCreate(&cuHandle));

    // Set default GPU device
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int num_devices;
    m_cudart_check_errors(cudaGetDeviceCount(&num_devices));
    m_cudart_check_errors(cudaSetDevice(rank % num_devices));
}

inline void finalize_backend() {
    m_cublas_check_errors(cublasDestroy(cuHandle));
}

inline void flups_device_sync() {
    m_cudart_check_errors(cudaDeviceSynchronize());
}

// Allocation functions
inline void flups_memset(void *ptr, int val, size_t size) {
    m_cudart_check_errors(cudaMemset(ptr, val, size));
    flups_device_sync();
}

inline void *flups_balloc(size_t size) {
    void  *ptr;
    size_t m_calloc_size_        = size + FLUPS_ALIGNMENT - 1;
    size_t m_calloc_padded_size_ = m_calloc_size_ - m_calloc_size_ % FLUPS_ALIGNMENT;
    m_cudart_check_errors(cudaMalloc(&ptr, m_calloc_padded_size_));

    // Set the memory to zero
    flups_memset(ptr, 0, m_calloc_padded_size_);
    return ptr;
}

inline void flups_bfree(void *ptr) {
    m_cudart_check_errors(cudaFree(ptr));
}

inline flups_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_kind kind) {
    flups_plan plan;
    // parameters for the r2r transform
    plan.n       = arch.dim;
    plan.dist    = arch.dist;
    plan.stride  = arch.stride;
    plan.howmany = arch.howmany;
    plan.kind    = kind;
    plan.is_r2r  = true;
    // check whether we need to plan a r2c or c2r (by default it's a r2c)
    cufftType type = CUFFT_D2Z;

    switch (kind) {
        // DCT Type I
        case FFTW_REDFT00:
            plan.dim     = 2 * plan.n - 2;
            plan.idist   = 2 * plan.n;  // should be 2n-2 (but cufft doesn't like it)
            plan.odist   = plan.n;
            plan.memsize = (2 * plan.n) * plan.howmany;
            break;
        // DST Type I
        case FFTW_RODFT00:
            plan.dim     = 2 * plan.n + 2;
            plan.idist   = 2 * plan.n + 4;  // should be 2n+2 (but cufft doesn't like it)
            plan.odist   = plan.n + 2;
            plan.memsize = (2 * plan.n + 4) * plan.howmany;
            break;
        // DCT/DST Type II
        case FFTW_REDFT10:
        case FFTW_RODFT10:
            plan.dim     = 2 * plan.n;
            plan.idist   = 2 * plan.n + 2;  // should be 2n (but cufft doesn't like it)
            plan.odist   = plan.n + 1;
            plan.memsize = (2 * plan.n + 2) * plan.howmany;
            break;
        // DCT/DST Type III
        case FFTW_REDFT01:
        case FFTW_RODFT01:
            plan.dim     = 2 * plan.n;
            plan.idist   = plan.n + 1;
            plan.odist   = 2 * plan.n + 2;  // here too !
            plan.memsize = (2 * plan.n + 2) * plan.howmany;
            // this is a c2r transform
            type = CUFFT_Z2D;
            break;
        // DCT/DST Type IV
        case FFTW_REDFT11:
        case FFTW_RODFT11:
            plan.dim     = 8 * plan.n;
            plan.idist   = 8 * plan.n + 2;  //  should be 8n (but cufft doesn't like it)
            plan.odist   = 4 * plan.n + 1;
            plan.memsize = (8 * plan.n + 2) * plan.howmany;
            break;
        default:
            FLUPS_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    // Note: some plans (typically involved in the computation of the Green function
    //       for mixed BC's) require a non-unit stride of 2. Here, we handle that
    //       directly during the pre- and post-processing (thus the plan has unit stride).
    m_cufft_check_errors(cufftPlanMany(&plan.cufft, 1, &plan.dim,
                                       &plan.dim, 1, plan.idist,
                                       &plan.dim, 1, plan.odist,
                                       type, plan.howmany));
    return plan;
}

inline flups_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_plan plan;
    m_cufft_check_errors(cufftPlanMany(&plan.cufft, 1, &arch.dim,
                                       &arch.dim, 1, arch.dist,
                                       &arch.dim, 1, arch.dist / 2,
                                       CUFFT_D2Z, arch.howmany));
    return plan;
}

inline flups_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_plan plan;
    m_cufft_check_errors(cufftPlanMany(&plan.cufft, 1, &arch.dim,
                                       &arch.dim, 1, arch.dist / 2,
                                       &arch.dim, 1, arch.dist,
                                       CUFFT_Z2D, arch.howmany));
    return plan;
}

inline flups_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    flups_plan plan;
    m_cufft_check_errors(cufftPlanMany(&plan.cufft, 1, &arch.dim,
                                       &arch.dim, 1, arch.dist,
                                       &arch.dim, 1, arch.dist,
                                       CUFFT_Z2Z, arch.howmany));
    return plan;
}

inline flups_shuffle_plan flups_plan_shuffle(fftw_iodim dims[2], double *in, double *out, bool iscomplex) {
    flups_shuffle_plan plan;
    plan.n0 = dims[0].n;
    plan.n1 = dims[1].n;
    return plan;
}

// Execute functions
inline void flups_execute_c2c(flups_plan plan, flups_complex *in, flups_complex *out, int sign) {
    m_cufft_check_errors(cufftExecZ2Z(plan.cufft, in, out, sign));
    flups_device_sync();
}

inline void flups_execute_r2c(flups_plan plan, double *in, flups_complex *out) {
    m_cufft_check_errors(cufftExecD2Z(plan.cufft, in, out));
    flups_device_sync();
}

inline void flups_execute_c2r(flups_plan plan, flups_complex *in, double *out) {
    m_cufft_check_errors(cufftExecZ2D(plan.cufft, in, out));
    flups_device_sync();
}

inline void flups_execute_r2r(flups_plan plan, double *in, double *out) {
    int n       = plan.n;
    int dim     = plan.dim;
    int dist    = plan.dist;
    int idist   = plan.idist;
    int odist   = plan.odist;
    int stride  = plan.stride;
    int howmany = plan.howmany;

    // Working memory
    double *mem = plan.mem;

    switch (plan.kind) {
        // DCT Type I
        case FFTW_REDFT00: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + i] = in[j * dist + stride * i];

                if (i == 0 || i == n - 1) {
                    mem[j * idist + i] *= .5;
                }
            };
            for_loop(pre_op, n * howmany, "DCT Type I pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = 2.0 * mem[2 * (j * odist + i)];
            };
            for_loop(post_op, n * howmany, "DCT Type I post-processing");

            break;
        }
        // DCT Type II
        case FFTW_REDFT10: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + i] = in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DCT Type II pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                   = -i * M_PI / (2.0 * n);
                out[j * dist + stride * i] = 2.0 * (mem[2 * (j * odist + i)] * cos(t) - mem[2 * (j * odist + i) + 1] * sin(t));
            };
            for_loop(post_op, n * howmany, "DCT Type II post-processing");

            break;
        }
        // DCT Type III
        case FFTW_REDFT01: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                     = -i * M_PI / (2.0 * n);
                mem[2 * (j * idist + i)]     = in[j * dist + stride * i] * cos(t);
                mem[2 * (j * idist + i) + 1] = in[j * dist + stride * i] * sin(t);
            };
            for_loop(pre_op, n * howmany, "DCT Type III pre-processing");

            flups_execute_c2r(plan, (flups_complex *)mem, mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = mem[1 + j * odist + i];
            };
            for_loop(post_op, n * howmany, "DCT Type III post-processing");

            break;
        }
        // DCT Type IV
        case FFTW_REDFT11: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + (2 * i + 1)]           = in[j * dist + stride * i];
                mem[j * idist + (4 * n - (2 * i + 1))] = -in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DCT Type IV pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = mem[2 * (j * odist + (2 * i + 1))];
            };
            for_loop(post_op, n * howmany, "DCT Type IV post-processing");

            break;
        }
        // DST Type I
        case FFTW_RODFT00: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[1 + j * idist + i] = in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DST Type I pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                out[j * dist + stride * i] = -2.0 * mem[3 + 2 * (j * odist + i)];
            };
            for_loop(post_op, n * howmany, "DST Type I post-processing");

            break;
        }
        // DST Type II
        case FFTW_RODFT10: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                const int sign     = 2 * (i % 2 == 0) - 1;
                mem[j * idist + i] = sign * in[j * dist + stride * i];
            };
            for_loop(pre_op, n * howmany, "DST Type II pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                             = -i * M_PI / (2.0 * n);
                out[j * dist + (n - 1) - stride * i] = 2.0 * (mem[2 * (j * odist + i)] * cos(t) - mem[2 * (j * odist + i) + 1] * sin(t));
            };
            for_loop(post_op, n * howmany, "DST Type II post-processing");

            break;
        }
        // DST Type III
        case FFTW_RODFT01: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                double t                     = -i * M_PI / (2.0 * n);
                mem[2 * (j * idist + i)]     = in[j * dist + (n - 1) - stride * i] * cos(t);
                mem[2 * (j * idist + i) + 1] = in[j * dist + (n - 1) - stride * i] * sin(t);
            };
            for_loop(pre_op, n * howmany, "DST Type III pre-processing");

            flups_execute_c2r(plan, (flups_complex *)mem, mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                const int sign             = 2 * (i % 2 == 0) - 1;
                out[j * dist + stride * i] = sign * mem[1 + j * odist + i];
            };
            for_loop(post_op, n * howmany, "DST Type III post-processing");

            break;
        }
        // DST Type IV
        case FFTW_RODFT11: {
            // Pre-processing
            auto pre_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                mem[j * idist + (2 * i + 1)]           = in[j * dist + (n - 1) - stride * i];
                mem[j * idist + (4 * n - (2 * i + 1))] = -in[j * dist + (n - 1) - stride * i];
            };
            for_loop(pre_op, n * howmany, "DST Type IV pre-processing");

            flups_execute_r2c(plan, mem, (flups_complex *)mem);

            // Post-processing
            auto post_op = m_kernel(const int idx) {
                const int i = idx % n;
                const int j = idx / n;

                const int sign             = 2 * (i % 2 == 0) - 1;
                out[j * dist + stride * i] = sign * mem[2 * (j * odist + (2 * i + 1))];
            };
            for_loop(post_op, n * howmany, "DST Type IV post-processing");

            break;
        }
        default:
            FLUPS_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    // Reset memory
    flups_memset(mem, 0, plan.memsize * sizeof(double));
}

inline void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex) {
    if (!iscomplex) {
        double alpha = 1.0;
        m_cublas_check_errors(cublasDgeam(cuHandle,
                                          CUBLAS_OP_T, CUBLAS_OP_N,
                                          plan.n0, plan.n1,
                                          &alpha, in, plan.n1,
                                          NULL, NULL, plan.n0,
                                          out, plan.n0));
    } else {
        flups_complex alpha = {1.0, 0.0};
        m_cublas_check_errors(cublasZgeam(cuHandle,
                                          CUBLAS_OP_T, CUBLAS_OP_N,
                                          plan.n0, plan.n1,
                                          &alpha, (flups_complex *)in, plan.n1,
                                          NULL, NULL, plan.n0,
                                          (flups_complex *)out, plan.n0));
    }
    flups_device_sync();
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    m_cufft_check_errors(cufftDestroy(plan.cufft));
}

inline void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    // nothing needs to be done here
    ((void)0);
}

#endif
