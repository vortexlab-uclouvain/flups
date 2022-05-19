#ifndef BACKEND_BACKEND_FFTW_HPP_
#define BACKEND_BACKEND_FFTW_HPP_

#include "omp.h"
#include "fftw3.h"
#include "defines.hpp"

// Datatypes
struct flups_fft_arch {
    int dim;
    int howmany;
    int stride;
    int dist;
};

using flups_kind         = fftw_r2r_kind;
using flups_plan         = fftw_plan;
using flups_shuffle_plan = fftw_plan;
using flups_complex      = fftw_complex;

// Init & cleanup backend
inline void init_backend() {
    fftw_init_threads();
}

inline void finalize_backend() {
    fftw_cleanup_threads();
    fftw_cleanup();
}

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

inline flups_shuffle_plan flups_plan_shuffle(fftw_iodim dims[2], double *in, double *out, bool iscomplex) {
    if (!iscomplex) {
        return fftw_plan_guru_r2r(0, NULL, 2, dims, in, out, NULL, FFTW_FLAG);
    } else {
        return fftw_plan_guru_dft(0, NULL, 2, dims, (flups_complex*)in, (flups_complex*)out, FLUPS_FORWARD, FFTW_FLAG);
    }
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

inline void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex) {
    // We do this shuffle in-place
    if (!iscomplex) {
        fftw_execute_r2r(plan, in, out);
    } else {
        fftw_execute_dft(plan, (opt_complex_ptr)(in), (opt_complex_ptr)(out));

    }
}

// Cleanup functions
inline void flups_destroy_plan(flups_plan plan) {
    fftw_destroy_plan(plan);
}

inline void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    fftw_destroy_plan(plan);
}

#endif
