#if defined(HAVE_HIP)

#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>
#include <rocfft/rocfft.h>

#include <cstring>  // for NULL

#include "backend.hpp"
#include "defines.hpp"

// rocBLAS handle
rocblas_handle rocHandle;

// Init & cleanup backend
void backend_init() {
    // Initialize rocBLAS
    m_rocblas_check_errors(rocblas_create_handle(&rocHandle));
}

void backend_finalize() {
    m_rocblas_check_errors(rocblas_destroy_handle(rocHandle));
}

void backend_sync() {
    m_hiprt_check_errors(hipDeviceSynchronize());
}

void backend_sync_stream(stream_t stream) {
    m_hiprt_check_errors(hipStreamSynchronize(stream));
}

// Allocation functions
void backend_memset(void *ptr, int val, size_t size, stream_t stream) {
    m_hiprt_check_errors(hipMemsetAsync(ptr, val, size, stream));
}

void *backend_malloc(size_t size) {
    void  *ptr;
    size_t m_calloc_size_        = size + FLUPS_ALIGNMENT - 1;
    size_t m_calloc_padded_size_ = m_calloc_size_ - m_calloc_size_ % FLUPS_ALIGNMENT;
    m_hiprt_check_errors(hipMalloc(&ptr, m_calloc_padded_size_));

    return ptr;
}

void backend_free(void *ptr) {
    m_hiprt_check_errors(hipFree(ptr));
}

// Return alignment of pointer
size_t backend_alignment_of(double *ptr) {
    BACKEND_CHECK(memTypeDevice == backend_get_memtype(ptr), "Pointer is not on device");
    return ((size_t)ptr % FLUPS_ALIGNMENT);
}

// Memcopy function - from host to device and vice-versa
void backend_memcpy(double *dst, double *src, size_t size, memcpyKind_t kind, stream_t stream) {
    m_hiprt_check_errors(hipMemcpyAsync(dst, src, size, (hipMemcpyKind)kind, stream));
}

// Helper function to find where the memory lives
memType_t backend_get_memtype(const void *ptr) {
    hipPointerAttribute_t attr;
    // Try to get attributes; some HIP versions behave differently
    hipError_t herr = hipPointerGetAttributes(&attr, ptr);
    if (herr != hipSuccess) {
        return (NULL != ptr) ? memTypeHost : memTypeUnregistered;
    }
#if HIP_VERSION >= 40000000
    const memType_t memtype = (memType_t)attr.memoryType;
#else
    const memType_t memtype = (memType_t)attr.type;
#endif
    return (memTypeUnregistered == memtype && NULL != ptr) ? memTypeHost : memtype;
}

// Helper: create a rocFFT plan for 1D transforms with optional batch and strides
static rocfft_plan create_rocfft_plan_r2c(int n, int batch, int istride, int ostride, int idist, int odist) {
    rocfft_plan p      = nullptr;
    size_t      length = (size_t)n;
    // Use rocfft_plan_create with placement not-inplace, real->complex, double precision
    m_rocfft_check_errors(rocfft_plan_create(&p,
                                             rocfft_placement_inplace,
                                             rocfft_transform_type_real_forward,
                                             rocfft_precision_double,
                                             1,
                                             &length,
                                             (size_t)batch));
    return p;
}

static rocfft_plan create_rocfft_plan_c2r(int n, int batch, int istride, int ostride, int idist, int odist) {
    rocfft_plan p      = nullptr;
    size_t      length = (size_t)n;
    m_rocfft_check_errors(rocfft_plan_create(&p,
                                             rocfft_placement_inplace,
                                             rocfft_transform_type_real_inverse,
                                             rocfft_precision_double,
                                             1,
                                             &length,
                                             (size_t)batch));
    return p;
}

static rocfft_plan create_rocfft_plan_c2c(int n, int batch, int istride, int ostride, int idist, int odist, int sign) {
    rocfft_plan           p      = nullptr;
    size_t                length = (size_t)n;
    rocfft_transform_type tt     = (sign == FLUPS_FORWARD) ? rocfft_transform_type_complex_forward : rocfft_transform_type_complex_inverse;
    m_rocfft_check_errors(rocfft_plan_create(&p,
                                             rocfft_placement_inplace,
                                             tt,
                                             rocfft_precision_double,
                                             1,
                                             &length,
                                             (size_t)batch));
    return p;
}

// Plan creation functions
flups_fft_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_fft_plan plan;

    // Create rocFFT plan (1D batched)
    rocfft_plan p = create_rocfft_plan_r2c(arch.dim, arch.howmany, arch.stride, arch.dist / 2, arch.dist, arch.dist / 2);

    plan.fft     = p;
    plan.dim     = arch.dim;
    plan.howmany = arch.howmany;
    plan.dist    = arch.dist;
    plan.stride  = arch.stride;
    plan.is_r2r  = false;

    return plan;
}

flups_fft_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_fft_plan plan;

    rocfft_plan p = create_rocfft_plan_c2r(arch.dim, arch.howmany, arch.stride, arch.dist / 2, arch.dist / 2, arch.dist);

    plan.fft     = p;
    plan.dim     = arch.dim;
    plan.howmany = arch.howmany;
    plan.dist    = arch.dist;
    plan.stride  = arch.stride;
    plan.is_r2r  = false;

    return plan;
}

flups_fft_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    flups_fft_plan plan;

    rocfft_plan p = create_rocfft_plan_c2c(arch.dim, arch.howmany, arch.stride, arch.dist, arch.dist, arch.dist, sign);

    plan.fft     = p;
    plan.dim     = arch.dim;
    plan.howmany = arch.howmany;
    plan.dist    = arch.dist;
    plan.stride  = arch.stride;
    plan.is_r2r  = false;

    return plan;
}

flups_shuffle_plan flups_plan_shuffle(flups_iodim dims[2], double *in, double *out, bool iscomplex) {
    flups_shuffle_plan plan;
    plan.n0 = dims[0].n;
    plan.n1 = dims[1].n;
    return plan;
}

// Execute functions
void flups_execute_c2c(flups_fft_plan plan, flups_complex *in, flups_complex *out, int sign, stream_t stream) {
    BACKEND_CHECK(plan.fft != NULL, "flups_execute_c2c: plan not created");

    rocfft_execution_info info = nullptr;
    m_rocfft_check_errors(rocfft_execution_info_create(&info));
    m_rocfft_check_errors(rocfft_execution_info_set_stream(info, stream));

    void *inbuf  = (void *)in;
    void *outbuf = (void *)out;

    m_rocfft_check_errors(rocfft_execute(plan.fft, &inbuf, &outbuf, info));

    m_rocfft_check_errors(rocfft_execution_info_destroy(info));
}

void flups_execute_r2c(flups_fft_plan plan, double *in, flups_complex *out, stream_t stream) {
    BACKEND_CHECK(plan.fft != NULL, "flups_execute_r2c: plan not created");

    rocfft_execution_info info = nullptr;
    m_rocfft_check_errors(rocfft_execution_info_create(&info));
    m_rocfft_check_errors(rocfft_execution_info_set_stream(info, stream));

    void *inbuf  = (void *)in;
    void *outbuf = (void *)out;

    m_rocfft_check_errors(rocfft_execute(plan.fft, &inbuf, &outbuf, info));

    m_rocfft_check_errors(rocfft_execution_info_destroy(info));
}

void flups_execute_c2r(flups_fft_plan plan, flups_complex *in, double *out, stream_t stream) {
    BACKEND_CHECK(plan.fft != NULL, "flups_execute_c2r: plan not created");

    rocfft_execution_info info = nullptr;
    m_rocfft_check_errors(rocfft_execution_info_create(&info));
    m_rocfft_check_errors(rocfft_execution_info_set_stream(info, stream));

    void *inbuf  = (void *)in;
    void *outbuf = (void *)out;

    m_rocfft_check_errors(rocfft_execute(plan.fft, &inbuf, &outbuf, info));

    m_rocfft_check_errors(rocfft_execution_info_destroy(info));
}

void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex, stream_t stream) {
    m_rocblas_check_errors(rocblas_set_stream(rocHandle, stream));

    if (!iscomplex) {
        double alpha = 1.0;
        double beta  = 0.0;
        m_rocblas_check_errors(rocblas_dgeam(rocHandle,
                                             rocblas_operation_transpose, rocblas_operation_none,
                                             plan.n0, plan.n1,
                                             &alpha, in, plan.n1,
                                             &beta, nullptr, plan.n0,
                                             out, plan.n0));
    } else {
        rocblas_double_complex alpha{1.0, 0.0};
        rocblas_double_complex beta{0.0, 0.0};
        m_rocblas_check_errors(rocblas_zgeam(rocHandle,
                                             rocblas_operation_transpose, rocblas_operation_none,
                                             plan.n0, plan.n1,
                                             &alpha, (const rocblas_double_complex *)in, plan.n1,
                                             &beta, nullptr, plan.n0,
                                             (rocblas_double_complex *)out, plan.n0));
    }
}

// Cleanup functions
void flups_destroy_plan(flups_fft_plan plan) {
    if (plan.fft != NULL) {
        m_rocfft_check_errors(rocfft_plan_destroy(plan.fft));
    }
}

void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    // nothing needs to be done here
    ((void)0);
}

#endif  // HAVE_HIP
