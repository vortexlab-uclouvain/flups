#ifndef BACKEND_HPP
#define BACKEND_HPP

#if defined(BACKEND_CUDA)
#define M_HAVE_GPU 1
#include "backend_cuda.hpp"
#elif defined(HAVE_HIP)
#define M_HAVE_GPU 1
#error "HIP backend not yet implemented"
#elif defined(BACKEND_FFTW)
#define M_HAVE_GPU 0
#include "backend_fftw.hpp"
#else
#error "Unknown or unspecified backend; please choose one of CUDA, HIP or FFTW"
#endif

// Bakcend init & cleanup functions
void backend_init();
void backend_finalize();

// Synchronization functions
void backend_sync();
void backend_sync_stream(stream_t stream);

// Allocation functions
void  backend_memset(void *ptr, int val, size_t size, stream_t stream);
void *backend_malloc(size_t size);
void  backend_free(void *ptr);

// Memory management functions
size_t    backend_alignment_of(double *ptr);
void      backend_memcpy(double *dst, double *src, size_t size, memcpyKind_t kind, stream_t stream);
memType_t backend_get_memtype(const void *ptr);

// Plan creation functions
flups_fft_plan     flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_r2r_kind_t kind);
flups_fft_plan     flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out);
flups_fft_plan     flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out);
flups_fft_plan     flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign);
flups_shuffle_plan flups_plan_shuffle(flups_iodim dims[2], double *in, double *out, bool iscomplex);

// Execute functions
void flups_execute_r2r(flups_fft_plan plan, double *in, double *out, stream_t stream);
void flups_execute_r2c(flups_fft_plan plan, double *in, flups_complex *out, stream_t stream);
void flups_execute_c2r(flups_fft_plan plan, flups_complex *in, double *out, stream_t stream);
void flups_execute_c2c(flups_fft_plan plan, flups_complex *in, flups_complex *out, int sign, stream_t stream);
void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex, stream_t stream);

// Cleanup functions
void flups_destroy_plan(flups_fft_plan plan);
void flups_destroy_shuffle_plan(flups_shuffle_plan plan);

#endif  // BACKEND_HPP