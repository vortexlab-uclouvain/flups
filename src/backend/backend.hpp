#ifndef BACKEND_HPP
#define BACKEND_HPP

#define FLUPS_FFTW_BACKEND

#if defined(FLUPS_FFTW_BACKEND)
#include "fftw3.h"

// Datatypes
using flups_plan    = fftw_plan;
using flups_complex = fftw_complex;

// Allocation functions
inline flups_plan flups_plan_r2c_1d(int n, double *in, flups_complex *out) {
    return fftw_plan_dft_r2c_1d(n, in, out, FFTW_FLAG);
}

inline flups_plan flups_plan_c2r_1d(int n, flups_complex *in, double *out) {                                          \
    return fftw_plan_dft_c2r_1d(n, in, out, FFTW_FLAG);
}

inline flups_plan flups_plan_c2c_1d(int n, flups_complex *in, flups_complex *out, int sign) {
    return fftw_plan_dft_1d(n, in, out, sign, FFTW_FLAG);
}

// Execute functions
inline void flups_execute_c2c(flups_plan plan, flups_complex *in, flups_complex *out, int sign) {                            \
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
inline flups_plan flups_plan_r2c_1d(int n, double *in, flups_complex *out) {
   cufftHandle plan;
    cufftPlan1d(&plan, n, CUFFT_D2Z, 1);
    return plan; 
}

inline flups_plan flups_plan_c2r_1d(int n, flups_complex *in, double *out) {                                          \
    cufftHandle plan;
    cufftPlan1d(&plan, n, CUFFT_Z2D, 1);
    return plan;
}

inline flups_plan flups_plan_c2c_1d(int n, flups_complex *in, flups_complex *out, int sign) {
    cufftHandle plan;
    cufftPlan1d(&plan, n, CUFFT_Z2Z, 1);
    return plan;
}

// Execute functions
inline void flups_execute_c2c(flups_plan plan, flups_complex *in, flups_complex *out, int sign) {                            \
    cufftExecZ2Z(plan, in, out, sign);
}

inline void flups_execute_r2c(flups_plan plan, double *in, flups_complex *out) {
    cufftExecD2Z(plan, in, out);
}

inline void flups_execute_c2r(flups_plan plan, flups_complex *in, double *out) {
    cufftExecZ2D(plan, in, out);
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