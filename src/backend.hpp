/**
 * @file backend.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef BACKEND_HPP
#define BACKEND_HPP

#if defined(BACKEND_CUDA)
#define M_HAVE_GPU 1
#include "backend_cuda.hpp"
#elif defined(BACKEND_HIP)
#define M_HAVE_GPU 1
#include "backend_hip.hpp"
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

// Stream events: record marks the position in the stream, sync blocks the HOST until every operation
// enqueued before the record has completed (a no-op on the FFTW backend)
event_t backend_event_create();
void    backend_event_destroy(event_t event);
void    backend_event_record(event_t event, stream_t stream);
void    backend_event_sync(event_t event);

// Allocation functions
void  backend_memset(void *ptr, int val, size_t size, stream_t stream);
void *backend_malloc(size_t size);
void  backend_free(void *ptr);

// Dense linear algebra
/**
 * @brief C = alpha * A * B^T (beta = 0), column-major, on the given stream
 *
 * A is (m x k) with leading dimension @p lda, B is (n x k) with leading dimension @p ldb and
 * C is (m x n) with leading dimension @p ldc, all column-major. The call is a no-op if any of
 * @p m, @p n or @p k is <= 0.
 *
 * @param m the number of rows of A and C
 * @param n the number of rows of B, i.e. the number of columns of C
 * @param k the number of columns of A and B, i.e. the contraction length
 * @param alpha the scaling of the product
 * @param A the (m x k) column-major matrix, `A[i + q * lda]`
 * @param lda the leading dimension of A, >= m
 * @param B the (n x k) column-major matrix, `B[j + q * ldb]`
 * @param ldb the leading dimension of B, >= n
 * @param C the (m x n) column-major result, `C[i + j * ldc]`, overwritten
 * @param ldc the leading dimension of C, >= m
 * @param stream the stream the operation is issued on
 */
void backend_dgemm_nt(const int m, const int n, const int k, const double alpha, const double *A, const int lda, const double *B, const int ldb, double *C, const int ldc, stream_t stream);

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

size_t flups_plan_worksize(const flups_fft_plan plan);
void   flups_plan_set_workarea(flups_fft_plan plan, void *work);

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
