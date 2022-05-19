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

using flups_kind    = fftw_r2r_kind;      // currently unused
using flups_plan    = cufftHandle;
using flups_complex = cufftDoubleComplex;

// Init & cleanup backend
inline void init_backend() {
    cublasCreate(&cuHandle);
}

inline void finalize_backend() {
    cublasDestroy(cuHandle);
}

// Allocation functions
inline flups_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_kind kind) {
    FLUPS_CHECK(false, "r2r transform are not yet implemented with CUDA backend !");
    cufftHandle plan;
    // TODO
    return plan;
}

inline flups_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    cufftHandle plan;
    cufftPlanMany(&plan, 1, &arch.dim,
                  &arch.dim, 1, arch.dist,
                  &arch.dim, 1, arch.dist / 2,
                  CUFFT_D2Z, arch.howmany);
    return plan;
}

inline flups_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    cufftHandle plan;
    cufftPlanMany(&plan, 1, &arch.dim,
                  &arch.dim, 1, arch.dist / 2,
                  &arch.dim, 1, arch.dist,
                  CUFFT_Z2D, arch.howmany);
    return plan;
}

inline flups_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    cufftHandle plan;
    cufftPlanMany(&plan, 1, &arch.dim,
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
    #pragma omp target data use_device_ptr(in, out)
    cufftExecZ2Z(plan, in, out, sign);
    cudaDeviceSynchronize();
}

inline void flups_execute_r2c(flups_plan plan, double *in, flups_complex *out) {
    #pragma omp target data use_device_ptr(in, out)
    cufftExecD2Z(plan, in, out);
    cudaDeviceSynchronize();
}

inline void flups_execute_c2r(flups_plan plan, flups_complex *in, double *out) {
    #pragma omp target data use_device_ptr(in, out)
    cufftExecZ2D(plan, in, out);
    cudaDeviceSynchronize();
}

inline void flups_execute_r2r(flups_plan plan, double *in, double *out, int sign) {
    FLUPS_CHECK(false, "r2r transform are not yet implemented with CUDA backend !");
}

inline void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex) {
    if (!iscomplex) {
        double alpha = 1.0;
        #pragma omp target data use_device_ptr(in, out)
        cublasDgeam(cuHandle,
                    CUBLAS_OP_T, CUBLAS_OP_N,
                    plan.n0, plan.n1,
                    &alpha, in, plan.n1,
                    NULL  , NULL,    plan.n0,
                    out, plan.n1);
    } else {
        flups_complex alpha = {1.0, 0.0};
        #pragma omp target data use_device_ptr(in, out)
        cublasZgeam(cuHandle,
                    CUBLAS_OP_T, CUBLAS_OP_N,
                    plan.n0, plan.n1,
                    &alpha, (flups_complex*)in, plan.n1,
                    NULL  ,                    NULL, plan.n0,
                    (flups_complex*)out, plan.n0);
    }
    cudaDeviceSynchronize();
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    cufftDestroy(plan);
}

inline void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    // nothing needs to be done here
    ((void)0);
}

#endif
