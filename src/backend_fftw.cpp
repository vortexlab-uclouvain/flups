#if defined(BACKEND_FFTW)

#include "backend.hpp"
#include "h3lpr/ptr.hpp"

// Define the FFTW memory allocator & free (using H3LPR)
#define m_fftw_calloc(size)                                                        \
    ({                                                                             \
        H3LPR::m_ptr<H3LPR::H3LPR_ALLOC_POSIX, void *, FLUPS_ALIGNMENT> ptr(size); \
                                                                                   \
        ptr();                                                                     \
    })

/**
 * @brief free the memory allocated using m_calloc
 *
 */
#define m_fftw_free(data) \
    {                     \
        std::free(data);  \
    }

// Init & cleanup backend
void backend_init() {
#ifdef FLUPS_WISDOM_PATH
    FLUPS_WARNING("Importing wisdom from %s", FLUPS_WISDOM_PATH);
    fftw_import_wisdom_from_filename(FLUPS_WISDOM_PATH);
#endif
}

void backend_finalize() {
    fftw_cleanup();
}

// Sync functions -- actually no ops on the CPU
void backend_sync() { ((void)0); }
void backend_sync_stream(stream_t stream) { ((void)0); }

// Events -- everything is synchronous on the CPU, nothing to do
event_t backend_event_create() { return 0; }
void    backend_event_destroy(event_t event) { ((void)0); }
void    backend_event_record(event_t event, stream_t stream) { ((void)0); }
void    backend_event_sync(event_t event) { ((void)0); }

// Allocation functions
void backend_memset(void *ptr, int val, size_t size, stream_t) {
    std::memset(ptr, val, size);
}

void *backend_malloc(size_t size) {
    return m_fftw_calloc(size);
}

// Dense linear algebra: C = alpha * A * B^T, column-major
void backend_dgemm_nt(const int m, const int n, const int k, const double alpha, const double *A, const int lda, const double *B, const int ldb, double *C, const int ldc, stream_t) {
    if (m <= 0 || n <= 0 || k <= 0) {
        return;
    }
    // C = 0
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            C[i + (size_t)j * ldc] = 0.0;
        }
    }
    // C += alpha * A * B^T, one rank-1 update per contraction index
    for (int q = 0; q < k; ++q) {
        for (int j = 0; j < n; ++j) {
            const double b = alpha * B[j + (size_t)q * ldb];
            for (int i = 0; i < m; ++i) {
                C[i + (size_t)j * ldc] += A[i + (size_t)q * lda] * b;
            }
        }
    }
}

// Return alignment of pointer
size_t backend_alignment_of(double *ptr) {
    return (size_t)fftw_alignment_of(ptr);
}

// Memcopy function - from host to device and vice-versa
void backend_memcpy(double *dst, double *src, size_t size, memcpyKind_t kind, stream_t stream) {
    BACKEND_CHECK(memcpyH2H == kind, "Only host to host memcpy is suppoerted here");
    std::memcpy(dst, src, size);
}

// Helper function to find where the memory lives
memType_t backend_get_memtype(const void *ptr) {
    if (NULL == ptr) {
        return memTypeUnregistered;
    }
    return memTypeHost;
}

flups_fft_plan flups_plan_r2r(flups_fft_arch arch, double *in, double *out, flups_r2r_kind_t kind) {
    return fftw_plan_many_r2r(1, &arch.dim, arch.howmany,
                              in, &arch.dim, arch.stride, arch.dist,
                              out, &arch.dim, arch.stride, arch.dist,
                              (fftw_r2r_kind *)&kind, FLUPS_FFTW_FLAG);
}

flups_fft_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    return fftw_plan_many_dft_r2c(1, &arch.dim, arch.howmany,
                                  in, &arch.dim, 1, arch.dist,
                                  out, &arch.dim, 1, arch.dist / 2,
                                  FLUPS_FFTW_FLAG);
}

flups_fft_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    return fftw_plan_many_dft_c2r(1, &arch.dim, arch.howmany,
                                  in, &arch.dim, 1, arch.dist / 2,
                                  out, &arch.dim, 1, arch.dist,
                                  FLUPS_FFTW_FLAG);
}

flups_fft_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    return fftw_plan_many_dft(1, &arch.dim, arch.howmany,
                              in, &arch.dim, 1, arch.dist,
                              out, &arch.dim, 1, arch.dist,
                              sign, FLUPS_FFTW_FLAG);
}

flups_shuffle_plan flups_plan_shuffle(flups_iodim dims[2], double *in, double *out, bool iscomplex) {
    if (!iscomplex) {
        return fftw_plan_guru_r2r(0, NULL, 2, dims, in, out, NULL, FLUPS_FFTW_FLAG);
    } else {
        return fftw_plan_guru_dft(0, NULL, 2, dims, (flups_complex *)in, (flups_complex *)out, FLUPS_FORWARD, FLUPS_FFTW_FLAG);
    }
}

// Work area: FFTW manages its own buffers, nothing to share
size_t flups_plan_worksize(const flups_fft_plan plan) { return 0; }
void   flups_plan_set_workarea(flups_fft_plan plan, void *work) { ((void)0); }

// Execute functions
void flups_execute_c2c(flups_fft_plan plan, flups_complex *in, flups_complex *out, int sign, stream_t) {
    fftw_execute_dft(plan, in, out);
}

void flups_execute_r2c(flups_fft_plan plan, double *in, flups_complex *out, stream_t) {
    fftw_execute_dft_r2c(plan, in, out);
}

void flups_execute_c2r(flups_fft_plan plan, flups_complex *in, double *out, stream_t) {
    fftw_execute_dft_c2r(plan, in, out);
}

void flups_execute_r2r(flups_fft_plan plan, double *in, double *out, stream_t) {
    fftw_execute_r2r(plan, in, out);
}

void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex, stream_t) {
    // We do this shuffle in-place
    if (!iscomplex) {
        fftw_execute_r2r(plan, in, out);
    } else {
        fftw_execute_dft(plan, (opt_complex_ptr)(in), (opt_complex_ptr)(out));
    }
}

// Cleanup functions
void backend_free(void *ptr) {
    m_fftw_free(ptr);
}

void flups_destroy_plan(flups_fft_plan plan) {
    fftw_destroy_plan(plan);
}

void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    fftw_destroy_plan(plan);
}

#endif  // BACKEND_FFTW
