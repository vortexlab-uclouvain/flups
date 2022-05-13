#ifndef BACKEND_HPP
#define BACKEND_HPP

#include "omp.h"
#include "fftw3.h"
#include "flups_interface.h"

// Datatypes
struct flups_fft_arch {
    int dim;
    int howmany;
    int stride;
    int dist;
};

using flups_kind    = fftw_r2r_kind;
using flups_plan    = fftw_plan;
using flups_complex = fftw_complex;

// Allocation functions
inline flups_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_kind kind) {
    fftw_plan_with_nthreads(omp_get_max_threads());
    return fftw_plan_many_r2r(1, &arch.dim, arch.howmany,
                              in, &arch.dim, arch.stride, arch.dist,
                              out, &arch.dim, arch.stride, arch.dist,
                              &kind, FFTW_FLAG);
}

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

#endif