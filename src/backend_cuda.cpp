#if defined(BACKEND_CUDA)

#include "backend.hpp"
#include "defines.hpp"

cublasHandle_t cuHandle;

// Init & cleanup backend
void backend_init() {
    // // Set default GPU device
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // int num_devices;
    // m_cudart_check_errors(cudaGetDeviceCount(&num_devices));
    // m_cudart_check_errors(cudaSetDevice(rank % num_devices));

    //  Initialize cuBLAS
    m_cublas_check_errors(cublasCreate(&cuHandle));
}

void backend_finalize() {
    m_cublas_check_errors(cublasDestroy(cuHandle));
}

void backend_sync() {
    m_cudart_check_errors(cudaDeviceSynchronize());
}

void backend_sync_stream(stream_t stream) {
    m_cudart_check_errors(cudaStreamSynchronize(stream));
}

// Events (no timing: cheaper, and we only use them for ordering)
event_t backend_event_create() {
    event_t event;
    m_cudart_check_errors(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
    return event;
}
void backend_event_destroy(event_t event) {
    m_cudart_check_errors(cudaEventDestroy(event));
}
void backend_event_record(event_t event, stream_t stream) {
    m_cudart_check_errors(cudaEventRecord(event, stream));
}
void backend_event_sync(event_t event) {
    m_cudart_check_errors(cudaEventSynchronize(event));
}

// Allocation functions
void backend_memset(void *ptr, int val, size_t size, stream_t stream) {
    m_cudart_check_errors(cudaMemsetAsync(ptr, val, size, stream));
}

void *backend_malloc(size_t size) {
    void  *ptr;
    size_t m_calloc_size_        = size + FLUPS_ALIGNMENT - 1;
    size_t m_calloc_padded_size_ = m_calloc_size_ - m_calloc_size_ % FLUPS_ALIGNMENT;
    m_cudart_check_errors(cudaMalloc(&ptr, m_calloc_padded_size_));

    return ptr;
}

void backend_free(void *ptr) {
    m_cudart_check_errors(cudaFree(ptr));
}

// Dense linear algebra: C = alpha * A * B^T, column-major
void backend_dgemm_nt(const int m, const int n, const int k, const double alpha, const double *A, const int lda, const double *B, const int ldb, double *C, const int ldc, stream_t stream) {
    if (m <= 0 || n <= 0 || k <= 0) {
        return;
    }
    const double beta = 0.0;
    m_cublas_check_errors(cublasSetStream(cuHandle, stream));
    m_cublas_check_errors(cublasDgemm(cuHandle,
                                      CUBLAS_OP_N, CUBLAS_OP_T,
                                      m, n, k,
                                      &alpha, A, lda,
                                      B, ldb,
                                      &beta, C, ldc));
}

// Return alignment of pointer
// TODO: what is the proper alignment on GPUs, and how to implement it ?
size_t backend_alignment_of(double *ptr) {
    BACKEND_CHECK(memTypeDevice == backend_get_memtype(ptr), "Pointer is not on device");
    return ((size_t)ptr % FLUPS_ALIGNMENT);
}

// Memcopy function - from host to device and vice-versa
void backend_memcpy(double *dst, double *src, size_t size, memcpyKind_t kind, stream_t stream) {
    m_cudart_check_errors(cudaMemcpyAsync(dst, src, size, (cudaMemcpyKind)kind, stream));
}

// Helper function to find where the memory lives
memType_t backend_get_memtype(const void *ptr) {
    cudaPointerAttributes attr;
    m_cudart_check_errors(cudaPointerGetAttributes(&attr, ptr));
    const memType_t memtype = (memType_t)attr.type;
    // Host pointer returns cudaMemoryTypeUnregistered --> turn it into host if ptr is not NULL
    return (memTypeUnregistered == memtype && NULL != ptr) ? memTypeHost : memtype;
}

/**
 * @brief create a batched 1D cuFFT plan WITHOUT its own work area
 *
 * The work area is provided later by the Solver (one buffer shared by every plan, they all
 * execute serialized on the same stream), see flups_plan_set_workarea. The size needed by this
 * plan is returned in worksize.
 */
void backend_make_plan(cufftHandle *fft, size_t *worksize, const int dim, const int idist, const int odist, const cufftType type, const int howmany) {
    long long n = dim, id = idist, od = odist;
    m_cufft_check_errors(cufftCreate(fft));
    m_cufft_check_errors(cufftSetAutoAllocation(*fft, 0));
    m_cufft_check_errors(cufftMakePlanMany64(*fft, 1, &n,
                                             &n, 1, id,
                                             &n, 1, od,
                                             type, howmany, worksize));
}

size_t flups_plan_worksize(const flups_fft_plan plan) {
    return plan.worksize;
}

void flups_plan_set_workarea(flups_fft_plan plan, void *work) {
    m_cufft_check_errors(cufftSetWorkArea(plan.fft, work));
}

flups_fft_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_fft_plan plan;
    backend_make_plan(&plan.fft, &plan.worksize, arch.dim, arch.dist, arch.dist / 2, CUFFT_D2Z, arch.howmany);
    return plan;
}

flups_fft_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_fft_plan plan;
    backend_make_plan(&plan.fft, &plan.worksize, arch.dim, arch.dist / 2, arch.dist, CUFFT_Z2D, arch.howmany);
    return plan;
}

flups_fft_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    flups_fft_plan plan;
    backend_make_plan(&plan.fft, &plan.worksize, arch.dim, arch.dist, arch.dist, CUFFT_Z2Z, arch.howmany);
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
    m_cufft_check_errors(cufftSetStream(plan.fft, stream));
    m_cufft_check_errors(cufftExecZ2Z(plan.fft, in, out, sign));
}

void flups_execute_r2c(flups_fft_plan plan, double *in, flups_complex *out, stream_t stream) {
    m_cufft_check_errors(cufftSetStream(plan.fft, stream));
    m_cufft_check_errors(cufftExecD2Z(plan.fft, in, out));
}

void flups_execute_c2r(flups_fft_plan plan, flups_complex *in, double *out, stream_t stream) {
    m_cufft_check_errors(cufftSetStream(plan.fft, stream));
    m_cufft_check_errors(cufftExecZ2D(plan.fft, in, out));
}

void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex, stream_t stream) {
    m_cublas_check_errors(cublasSetStream(cuHandle, stream));

    // out = alpha * in^T + beta * B with beta = 0: B is not read (NULL), but the scalar must be a valid pointer
    if (!iscomplex) {
        const double alpha = 1.0;
        const double beta  = 0.0;
        m_cublas_check_errors(cublasDgeam(cuHandle,
                                          CUBLAS_OP_T, CUBLAS_OP_N,
                                          plan.n0, plan.n1,
                                          &alpha, in, plan.n1,
                                          &beta, NULL, plan.n0,
                                          out, plan.n0));
    } else {
        const flups_complex alpha = {1.0, 0.0};
        const flups_complex beta  = {0.0, 0.0};
        m_cublas_check_errors(cublasZgeam(cuHandle,
                                          CUBLAS_OP_T, CUBLAS_OP_N,
                                          plan.n0, plan.n1,
                                          &alpha, (flups_complex *)in, plan.n1,
                                          &beta, NULL, plan.n0,
                                          (flups_complex *)out, plan.n0));
    }
}

// Cleanup functions
void flups_destroy_plan(flups_fft_plan plan) {
    m_cufft_check_errors(cufftDestroy(plan.fft));
    if (plan.twiddle != NULL) backend_free(plan.twiddle);
}

void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    // nothing needs to be done here
    ((void)0);
}

#endif  // BACKEND_CUDA
