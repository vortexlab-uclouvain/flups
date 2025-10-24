#if defined(BACKEND_CUDA)

#include "backend.hpp"
#include "defines.hpp"

cublasHandle_t cuHandle;
cudaStream_t   comp_stream = 0;

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

// Allocation functions
void backend_memset(void *ptr, int val, size_t size) {
    m_cudart_check_errors(cudaMemset(ptr, val, size));
    backend_sync();
}

void *backend_malloc(size_t size) {
    backend_sync_stream(comp_stream);
    void  *ptr;
    size_t m_calloc_size_        = size + FLUPS_ALIGNMENT - 1;
    size_t m_calloc_padded_size_ = m_calloc_size_ - m_calloc_size_ % FLUPS_ALIGNMENT;
    m_cudart_check_errors(cudaMalloc(&ptr, m_calloc_padded_size_));

    // Set the memory to zero
    backend_memset(ptr, 0, m_calloc_padded_size_);
    return ptr;
}

void backend_free(void *ptr) {
    m_cudart_check_errors(cudaFree(ptr));
}

// Return alignment of pointer
// TODO: what is the proper alignment on GPUs, and how to implement it ?
size_t backend_alignment_of(double *ptr) {
    BACKEND_CHECK(memTypeDevice == backend_get_memtype(ptr), "Pointer is not on device");
    return ((size_t)ptr % FLUPS_ALIGNMENT);
}

// Memcopy function - from host to device and vice-versa
void backend_memcpy(double *dst, double *src, size_t size, memcpyKind_t kind) {
    m_cudart_check_errors(cudaMemcpyAsync(dst, src, size, (cudaMemcpyKind)kind, comp_stream));
}

// Helper function to find where the memory lives
memType_t backend_get_memtype(const void *ptr) {
    cudaPointerAttributes attr;
    m_cudart_check_errors(cudaPointerGetAttributes(&attr, ptr));
    const memType_t memtype = (memType_t)attr.type;
    // Host pointer returns cudaMemoryTypeUnregistered --> turn it into host if ptr is not NULL
    return (memTypeUnregistered == memtype && NULL != ptr) ? memTypeHost : memtype;
}

flups_fft_plan flups_plan_r2c(flups_fft_arch arch, double *in, flups_complex *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_fft_plan plan;
    m_cufft_check_errors(cufftPlanMany(&plan.fft, 1, &arch.dim,
                                       &arch.dim, 1, arch.dist,
                                       &arch.dim, 1, arch.dist / 2,
                                       CUFFT_D2Z, arch.howmany));
    m_cufft_check_errors(cufftSetStream(plan.fft, comp_stream));
    return plan;
}

flups_fft_plan flups_plan_c2r(flups_fft_arch arch, flups_complex *in, double *out) {
    FLUPS_CHECK(arch.dist % 2 == 0, "the fft stride (or rather 'dist') should be even");
    flups_fft_plan plan;
    m_cufft_check_errors(cufftPlanMany(&plan.fft, 1, &arch.dim,
                                       &arch.dim, 1, arch.dist / 2,
                                       &arch.dim, 1, arch.dist,
                                       CUFFT_Z2D, arch.howmany));
    m_cufft_check_errors(cufftSetStream(plan.fft, comp_stream));
    return plan;
}

flups_fft_plan flups_plan_c2c(flups_fft_arch arch, flups_complex *in, flups_complex *out, int sign) {
    flups_fft_plan plan;
    m_cufft_check_errors(cufftPlanMany(&plan.fft, 1, &arch.dim,
                                       &arch.dim, 1, arch.dist,
                                       &arch.dim, 1, arch.dist,
                                       CUFFT_Z2Z, arch.howmany));
    m_cufft_check_errors(cufftSetStream(plan.fft, comp_stream));
    return plan;
}

flups_shuffle_plan flups_plan_shuffle(flups_iodim dims[2], double *in, double *out, bool iscomplex) {
    flups_shuffle_plan plan;
    plan.n0 = dims[0].n;
    plan.n1 = dims[1].n;
    return plan;
}

// Execute functions
void flups_execute_c2c(flups_fft_plan plan, flups_complex *in, flups_complex *out, int sign) {
    m_cufft_check_errors(cufftExecZ2Z(plan.fft, in, out, sign));
}

void flups_execute_r2c(flups_fft_plan plan, double *in, flups_complex *out) {
    m_cufft_check_errors(cufftExecD2Z(plan.fft, in, out));
}

void flups_execute_c2r(flups_fft_plan plan, flups_complex *in, double *out) {
    m_cufft_check_errors(cufftExecZ2D(plan.fft, in, out));
}

void flups_execute_shuffle(flups_shuffle_plan plan, double *in, double *out, bool iscomplex) {
    // Set stream
    m_cublas_check_errors(cublasSetStream(cuHandle, comp_stream));

    if (!iscomplex) {
        double alpha = 1.0;
        m_cublas_check_errors(cublasDgeam(cuHandle,
                                          CUBLAS_OP_T, CUBLAS_OP_N,
                                          plan.n0, plan.n1,
                                          &alpha, in, plan.n1,
                                          NULL, NULL, plan.n0,
                                          out, plan.n0));
    } else {
        flups_complex alpha = {1.0, 0.0};
        m_cublas_check_errors(cublasZgeam(cuHandle,
                                          CUBLAS_OP_T, CUBLAS_OP_N,
                                          plan.n0, plan.n1,
                                          &alpha, (flups_complex *)in, plan.n1,
                                          NULL, NULL, plan.n0,
                                          (flups_complex *)out, plan.n0));
    }
}

// Cleanup functions
void flups_destroy_plan(flups_fft_plan plan) {
    m_cufft_check_errors(cufftDestroy(plan.fft));
}

void flups_destroy_shuffle_plan(flups_shuffle_plan plan) {
    // nothing needs to be done here
    ((void)0);
}

#endif  // BACKEND_CUDA
