#ifndef BACKEND_BACKEND_CUDA_HPP_
#define BACKEND_BACKEND_CUDA_HPP_

#include "omp.h"
#include "fftw3.h"
#include "cufft.h"
#include "defines.hpp"
#include "cublas_v2.h"
#include "cuda_runtime.h"

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
    // cuFFT plan handle
    cufftHandle cufft;
    // parameters for the r2r transform
    int n;              // number of inputs/output real values
    int dist;           // distance between two successive fft's
    int howmany;        // number of fft's to perform
    flups_kind kind;    // kinf of the r2r transform
    // parameters for the equivalent r2c transform
    int dim;            // size of the equivalent r2c fft
    int idist;          // input  distance between two successive fft's
    int odist;          // output distance between two successive fft's
    int nmem;           // size of the needed memory space
};


// Init & cleanup backend
inline void init_backend() {
    cublasCreate(&cuHandle);

    // Set default GPU device
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    omp_set_default_device(rank % omp_get_num_devices());
}

inline void finalize_backend() {
    cublasDestroy(cuHandle);
}

// Allocation functions
inline void flups_memset(void *ptr, int val, size_t size) {
    cudaMemset(ptr, val, size);
}

inline void* flups_balloc(size_t size) {
    size_t m_calloc_size_        = size + M_ALIGNMENT - 1;
    size_t m_calloc_padded_size_ = m_calloc_size_ - m_calloc_size_ % M_ALIGNMENT;
    void*  ptr                   = omp_target_alloc(m_calloc_padded_size_, omp_get_default_device());

    flups_memset(ptr, 0, m_calloc_padded_size_);
    return ptr;
}

inline void flups_bfree(void *ptr) {
    omp_target_free(ptr, omp_get_default_device());
}


inline flups_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_kind kind) {
    flups_plan plan;
    // parameters for the r2r transform
    plan.n = arch.dim;
    plan.dist = arch.dist;
    plan.howmany = arch.howmany;
    plan.kind = kind;
    // check whether we need to plan a r2c or c2r (by default it's a r2c)
    cufftType type = CUFFT_D2Z;

    switch (kind) {
        // DCT Type I
        case FFTW_REDFT00:
            plan.dim   = 2*plan.n - 2;
            plan.idist = 2*plan.n;      // should be 2n-2 (but cufft doesn't like it)
            plan.odist =   plan.n;
            plan.nmem  = (2*plan.n) * plan.howmany;
            break;
        // DCT Type II
        case FFTW_REDFT10:
            plan.dim   = 2*plan.n;
            plan.idist = 2*plan.n + 2;  // should be 2n (but cufft doesn't like it)
            plan.odist =   plan.n + 1;
            plan.nmem  = (2*plan.n + 2) * plan.howmany;
            break;
        // DCT Type III
        case FFTW_REDFT01:
            plan.dim   = 2*plan.n;
            plan.idist =   plan.n + 1;
            plan.odist = 2*plan.n;
            plan.nmem  = (2*plan.n + 2) * plan.howmany;
            // this is a c2r transform
            type = CUFFT_Z2D;
            break;
        // DCT Type IV
        case FFTW_REDFT11:
            plan.dim   = 8*plan.n;
            plan.idist = 8*plan.n + 2;  //  should be 8n (but cufft doesn't like it)
            plan.odist = 4*plan.n + 1;
            plan.nmem  = (8*plan.n + 2) * plan.howmany;
            break;
        // DST Type I
        case FFTW_RODFT00:
            plan.dim   = 2*plan.n + 2;
            plan.idist = 2*plan.n + 4;  // should be 2n+2 (but cufft doesn't like it)
            plan.odist =   plan.n + 3;
            plan.nmem  = (2*plan.n + 6) * plan.howmany;
            break;
    default:
        FLUPS_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    cufftPlanMany(&plan.cufft, 1, &plan.dim,
                    &plan.dim, 1, plan.idist,
                    &plan.dim, 1, plan.odist,
                    type, plan.howmany);
    return plan;
}

inline flups_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_plan plan;
    cufftPlanMany(&plan.cufft, 1, &arch.dim,
                  &arch.dim, 1, arch.dist,
                  &arch.dim, 1, arch.dist / 2,
                  CUFFT_D2Z, arch.howmany);
    return plan;
}

inline flups_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_plan plan;
    cufftPlanMany(&plan.cufft, 1, &arch.dim,
                  &arch.dim, 1, arch.dist / 2,
                  &arch.dim, 1, arch.dist,
                  CUFFT_Z2D, arch.howmany);
    return plan;
}

inline flups_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    flups_plan plan;
    cufftPlanMany(&plan.cufft, 1, &arch.dim,
                  &arch.dim, 1, arch.dist,
                  &arch.dim, 1, arch.dist,
                  CUFFT_Z2Z, arch.howmany);
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
    cufftExecZ2Z(plan.cufft, in, out, sign);
    cudaDeviceSynchronize();
}

inline void flups_execute_r2c(flups_plan plan, double *in, flups_complex *out) {
    cufftExecD2Z(plan.cufft, in, out);
    cudaDeviceSynchronize();
}

inline void flups_execute_c2r(flups_plan plan, flups_complex *in, double *out) {
    cufftExecZ2D(plan.cufft, in, out);
    cudaDeviceSynchronize();
}

inline void flups_execute_r2r(flups_plan plan, double *in, double *out) {
    int n = plan.n;
    int dim = plan.dim;
    int dist = plan.dist;
    int idist = plan.idist;
    int odist = plan.odist;
    int howmany = plan.howmany;
 
    // Allocate memory
    double *mem = (double*) flups_balloc(plan.nmem*sizeof(double));

    switch (plan.kind) {
        // DCT Type I
        case FFTW_REDFT00:
        {
            // Pre-processing
            #pragma omp target teams distribute is_device_ptr(in, mem)
            for (int j = 0; j < howmany; ++j) {
                #pragma omp parallel for
                for (int i = 0; i < n; ++i) {
                    mem[j * idist + i] = in[j * dist + i];
                }

                mem[j * idist]         *= .5;
                mem[j * idist + n - 1] *= .5;
            }

            flups_execute_r2c(plan, mem, (flups_complex*)mem);

            // Post-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(out, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    out[j * dist + i] = 2.0*mem[2*(j * odist + i)];
                }
            }
            break;
        }
        // DCT Type II
        case FFTW_REDFT10:
        {
            // Pre-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(in, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    mem[j * idist + i] = in[j * dist + i];
                }
            }

            flups_execute_r2c(plan, mem, (flups_complex*)mem);

            // Post-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(out, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    double t = -i*M_PI / (2.0*n);
                    out[j * dist + i] = 2.0 * (mem[2*(j*odist + i)]*std::cos(t) - mem[2*(j*odist + i) + 1]*std::sin(t)) ;
                }
            }
            break;
        }
        // DCT Type III
        case FFTW_REDFT01:
        {
            // Pre-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(in, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    double t = -i*M_PI / (2.0*n);
                    mem[2*(j * idist + i)]     = in[j * dist + i]*std::cos(t);
                    mem[2*(j * idist + i) + 1] = in[j * dist + i]*std::sin(t);
                }
            }

            flups_execute_c2r(plan, (flups_complex*)mem, mem);

            // Post-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(out, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    out[j * dist + i] = mem[1 + j * odist + i];
                }
            }
            break;
        }
        // DCT Type IV
        case FFTW_REDFT11:
        {
            // Pre-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(in, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    mem[j * idist + (2*i+1)]         =  in[j * dist + i];
                    mem[j * idist + (4*n - (2*i+1))] = -in[j * dist + i];
                }
            }

            flups_execute_r2c(plan, mem, (flups_complex*)mem);

            // Post-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(out, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    out[j * dist + i] = mem[2 * (j*odist + (2*i + 1))];
                }
            }
            break;
        }
        // DST Type I
        case FFTW_RODFT00:
        {
            // Pre-processing
            #pragma omp target teams distribute parallel for is_device_ptr(in, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    mem[1 + j * idist + i] = in[j * dist + i];
                }
            }

            flups_execute_r2c(plan, mem, (flups_complex*)mem);

            // Post-processing
            #pragma omp target teams distribute parallel for collapse(2) is_device_ptr(out, mem)
            for (int j = 0; j < howmany; ++j) {
                for (int i = 0; i < n; ++i) {
                    out[j * dist + i] = -2.0*mem[3 + 2*(j * odist + i)];
                }
            }
            break;
        }
    default:
        FLUPS_CHECK(false, "r2r transform of kind %d are not yet implemented with CUDA backend !", plan.kind);
    }

    // Free memory
    flups_bfree(mem);
}

inline void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex) {
    if (!iscomplex) {
        double alpha = 1.0;
        cublasDgeam(cuHandle,
                    CUBLAS_OP_T, CUBLAS_OP_N,
                    plan.n0, plan.n1,
                    &alpha, in, plan.n1,
                    NULL, NULL, plan.n0,
                    out, plan.n0);
    } else {
        flups_complex alpha = {1.0, 0.0};
        cublasZgeam(cuHandle,
                    CUBLAS_OP_T, CUBLAS_OP_N,
                    plan.n0, plan.n1,
                    &alpha, (flups_complex*)in, plan.n1,
                    NULL, NULL, plan.n0,
                    (flups_complex*)out, plan.n0);
    }
    cudaDeviceSynchronize();
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    cufftDestroy(plan.cufft);
}

inline void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    // nothing needs to be done here
    ((void)0);
}

#endif
