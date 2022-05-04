#ifndef BACKEND_HPP
#define BACKEND_HPP


struct flups_fft_arch {
    int dim;
    int howmany;
    int dist;
};


#if defined(FLUPS_FFTW_BACKEND)
#include "fftw3.h"

// Datatypes
using flups_plan    = fftw_plan;
using flups_complex = fftw_complex;

// Allocation functions
inline flups_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    return fftw_plan_dft_r2c_1d(arch.dim, in, out, FFTW_FLAG);
}

inline flups_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    return fftw_plan_dft_c2r_1d(arch.dim, in, out, FFTW_FLAG);
}

inline flups_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    return fftw_plan_dft_1d(arch.dim, in, out, sign, FFTW_FLAG);
}

// Execute functions
inline void flups_execute_c2c(flups_plan plan, flups_complex *in, flups_complex *out, int sign) {
    fftw_execute_dft(plan, in, out);
}

inline void flups_execute_r2c(flups_plan plan, double *in, flups_complex *out) {
    fftw_execute_dft_r2c(plan, in, out);
}

inline void flups_execute_c2r(flups_plan plan, flups_complex *in, double *out) {
    fftw_execute_dft_c2r(plan, in, out);
}

inline void flups_execute_r2r(flups_plan plan, double *in, double *out, int sign) {
    fftw_execute_r2r(plan, in, out);
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    fftw_destroy_plan(plan);
}

#elif defined(FLUPS_FFTW_BATCH_BACKEND)
#include "fftw3.h"

// Datatypes
using flups_plan    = fftw_plan;
using flups_complex = fftw_complex;

// Allocation functions
inline flups_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    fftw_plan_with_nthreads(omp_get_max_threads());
    return fftw_plan_many_dft_r2c(1, &arch.dim, arch.howmany,
                              in, &arch.dim, 1, arch.dist,
                              out, &arch.dim, 1, arch.dist / 2,
                              FFTW_FLAG);
}

inline flups_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    fftw_plan_with_nthreads(omp_get_max_threads());
    return fftw_plan_many_dft_c2r(1, &arch.dim, arch.howmany,
                              in, &arch.dim, 1, arch.dist / 2,
                              out, &arch.dim, 1, arch.dist,
                              FFTW_FLAG);
}

inline flups_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    fftw_plan_with_nthreads(omp_get_max_threads());
    return fftw_plan_many_dft(1, &arch.dim, arch.howmany,
                              in, &arch.dim, 1, arch.dist,
                              out, &arch.dim, 1, arch.dist,
                              sign, FFTW_FLAG);
}

// Execute functions
inline void flups_execute_c2c(flups_plan plan, flups_complex *in, flups_complex *out, int sign) {
    fftw_execute_dft(plan, in, out);
}

inline void flups_execute_r2c(flups_plan plan, double *in, flups_complex *out) {
    fftw_execute_dft_r2c(plan, in, out);
}

inline void flups_execute_c2r(flups_plan plan, flups_complex *in, double *out) {
    fftw_execute_dft_c2r(plan, in, out);
}

inline void flups_execute_r2r(flups_plan plan, double *in, double *out, int sign) {
    fftw_execute_r2r(plan, in, out);
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    fftw_destroy_plan(plan);
}

#elif defined(FLUPS_CUDA_BACKEND)
#include "omp.h"
#include "cufft.h"
#include "cuda_runtime.h"

// Datatypes
using flups_plan    = cufftHandle;
using flups_complex = cufftDoubleComplex;

// Allocation functions
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
    // TODO
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    cufftDestroy(plan);
}

#else
#error "Unknown/no backend defined. Please specify a correct backend !"
#endif

#endif