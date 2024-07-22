#if defined(BACKEND_FFTW)

#include "backend.hpp"

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

// Allocation functions
void backend_memset(void *ptr, int val, size_t size) {
    std::memset(ptr, val, size);
}

void *backend_malloc(size_t size) {
    return m_calloc(size);
}

// Return alignment of pointer
size_t backend_alignment_of(double *ptr) {
    return (size_t)fftw_alignment_of(ptr);
}

// Memcopy function - from host to device and vice-versa
void backend_memcpy(double *dst, double *src, size_t size, memcpyKind_t kind) {
    FLUPS_CHECK(memcpyH2H == kind, "Only host to host memcpy is suppoerted here");
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

// Execute functions
void flups_execute_c2c(flups_fft_plan plan, flups_complex *in, flups_complex *out, int sign) {
    fftw_execute_dft(plan, in, out);
}

void flups_execute_r2c(flups_fft_plan plan, double *in, flups_complex *out) {
    fftw_execute_dft_r2c(plan, in, out);
}

void flups_execute_c2r(flups_fft_plan plan, flups_complex *in, double *out) {
    fftw_execute_dft_c2r(plan, in, out);
}

void flups_execute_r2r(flups_fft_plan plan, double *in, double *out) {
    fftw_execute_r2r(plan, in, out);
}

void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex) {
    // We do this shuffle in-place
    if (!iscomplex) {
        fftw_execute_r2r(plan, in, out);
    } else {
        fftw_execute_dft(plan, (opt_complex_ptr)(in), (opt_complex_ptr)(out));
    }
}

// Cleanup functions
void backend_free(void *ptr) {
    m_free(ptr);
}

void flups_destroy_plan(flups_fft_plan plan) {
    fftw_destroy_plan(plan);
}

void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    fftw_destroy_plan(plan);
}

#endif  // BACKEND_FFTW