/**
 * @file test_cufftmp.cpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#include <cufftMp.h>
#include <mpi.h>

#include <complex>
#include <vector>

int main() {
    // Initialize MPI, pick a device
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int ndevices;
    cudaGetDeviceCount(&ndevices);
    cudaSetDevice(rank % ndevices);

    // Allocate CPU memory
    size_t                           N = 32;
    std::vector<std::complex<float>> cpu_data((N / size) * N * N, {314, 0});

    // Create plan, and make a multi-process plan with communicator attachment
    cufftHandle plan = 0;
    size_t      workspace;
    cufftCreate(&plan);
    cufftMpMakePlan3d(plan, N, N, N, CUFFT_C2C, &comm, CUFFT_COMM_MPI, &workspace);

    // Allocate memory, copy CPU data to GPU
    cudaLibXtDesc *desc;
    cufftXtMalloc(plan, &desc, CUFFT_XT_FORMAT_INPLACE);
    cufftXtMemcpy(plan, desc, cpu_data.data(), CUFFT_COPY_HOST_TO_DEVICE);

    // Run C2C FFT Forward
    cufftXtExecDescriptor(plan, desc, desc, CUFFT_FORWARD);

    // Copy back to CPU
    std::fill(cpu_data.begin(), cpu_data.end(), std::complex<float>(0, 0));
    cufftXtMemcpy(plan, cpu_data.data(), desc, CUFFT_COPY_DEVICE_TO_HOST);

    // Data in cpu_data is now distributed along the Y dimension, of size N * (N / size) * N
    // Test by comparing the very first entry on rank 0
    if (rank == 0) {
        if (cpu_data[0].real() == 314 * N * N * N) {
            printf("PASSED\n");
        } else {
            printf("FAILED\n");
        }
    }

    // Cleanup
    cufftXtFree(desc);
    cufftDestroy(plan);
    MPI_Finalize();
}
