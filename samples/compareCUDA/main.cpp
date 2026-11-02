/**
 * @file main.cpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "flups.h"
#include "h3lpr/parser.hpp"
#include "h3lpr/profiler.hpp"

#ifdef HAVE_HEFFTE
#include "heffte.h"
#endif

#ifdef HAVE_CUFFTMP
#include "cufftMp.h"
#endif

int main(int argc, char *argv[]) {
    //-------------------------------------------------------------------------
    // Initialize MPI
    //-------------------------------------------------------------------------
    MPI_Init(&argc, &argv);

    int      rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    //--------------------------------------------------------------------------
    // Get info from the command line
    //--------------------------------------------------------------------------
    H3LPR::Parser parser(argc, (const char **)argv);
    const auto    arg_nglob = parser.GetValues<int, 3>("--nglob", "the global resolution, will be used for both ACCFFT and FLUPS", {64, 64, 64});
    const auto    arg_nproc = parser.GetValues<int, 3>("--nproc", "the proc distribution, for FLUPS only", {1, 1, 1});
    const auto    arg_dom   = parser.GetValues<double, 3>("--dom", "the size of the domain, must be compatible with nglob", {1.0, 1.0, 1.0});
    const int     n_iter    = parser.GetValue<int>("--niter", "the number of iterations to perform", 20);
    const int     n_warm    = parser.GetValue<int>("--warm", "the number of iterations to perform when warming up", 1);
    const bool    profile   = parser.GetFlag("--profile", "forward the profiler to flups");
    parser.Finalize();

    //--------------------------------------------------------------------------
    // Definition of the problem
    //--------------------------------------------------------------------------
    const int    nglob[3] = {arg_nglob[0], arg_nglob[1], arg_nglob[2]};
    const int    nproc[3] = {arg_nproc[0], arg_nproc[1], arg_nproc[2]};
    const double L[3]     = {arg_dom[0], arg_dom[1], arg_dom[2]};
    assert(1 == nproc[0]);

    // get the grid spacing
    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};

    // get the PER PER PER BC everywhere
    const FLUPS_CenterType center_type[3] = {CELL_CENTER, CELL_CENTER, CELL_CENTER};
    FLUPS_BoundaryType    *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            mybc[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
            mybc[id][is][0] = PER;
        }
    }

    // local topology; defined by FLUPS and used by the other solvers
    int topo_start[3];
    int topo_end[3];
    int spectral_start[3];
    int spectral_end[3];
    int topo_nloc[3];
    int spectral_nloc[3];

    //..........................................................................
    // Display
    flups_info(argc, argv);

    //--------------------------------------------------------------------------
    std::string     prof_name = "beatme_nglob" + std::to_string(nglob[0]) + "_" + std::to_string(nglob[1]) + "_" + std::to_string(nglob[2]) + "_nrank" + std::to_string(comm_size);
    H3LPR::Profiler prof(prof_name);

    //--------------------------------------------------------------------------
    if (rank == 0) {
        printf("--------------------------------------------------------------\n");
        printf("                   FLUPS arguments\n");
        printf("Flups Topology initialised with:\n");
        printf("    --nglob= %d %d %d \n", nglob[0], nglob[1], nglob[2]);
        printf("    --nproc= %d %d %d \n", nproc[0], nproc[1], nproc[2]);
        printf("    --L    = %f %f %f \n", L[0], L[1], L[2]);
        printf("    --h    = %f %f %f \n", h[0], h[1], h[2]);
        printf("--------------------------------------------------------------\n");
    }

    //--------------------------------------------------------------------------
    /** - Benchmark FLUPS */
    //--------------------------------------------------------------------------
    {
        //--------------------------------------------------------------------------
        /** - Initialize FLUPS */
        //--------------------------------------------------------------------------
        if (rank == 0) printf("Initialization of FLUPS\n");

        // create a real topology
        FLUPS_Profiler *flups_prof = (profile) ? (FLUPS_Profiler *)&prof : nullptr;
        FLUPS_Topology *topoTmp    = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);
        FLUPS_Solver   *mysolver   = flups_init_timed(topoTmp, mybc, h, L, NOD, center_type, flups_prof);

        // set the CHAT2 green type (even if it's not used)
        flups_set_greenType(mysolver, CHAT_2);
        flups_setup(mysolver, true);
        double *solFLU = flups_get_innerBuffer(mysolver);

        // to fill the data we use the inner topo
        FLUPS_Topology *topoIn = flups_get_innerTopo_physical(mysolver);
        // instruct the solver to skip the first ST
        flups_skip_firstSwitchtopo(mysolver);

        // Get some info from the topology (everyone is going to use this)
        flups_topo_get_istartGlob(topoIn, topo_start);
        topo_nloc[0] = flups_topo_get_nloc(topoIn, 0);
        topo_nloc[1] = flups_topo_get_nloc(topoIn, 1);
        topo_nloc[2] = flups_topo_get_nloc(topoIn, 2);
        topo_end[0]  = topo_start[0] + topo_nloc[0];
        topo_end[1]  = topo_start[1] + topo_nloc[1];
        topo_end[2]  = topo_start[2] + topo_nloc[2];

        // Get some info on the spectral topology
        FLUPS_Topology *topoSpec = flups_get_innerTopo_spectral(mysolver);
        flups_topo_get_istartGlob(topoSpec, spectral_start);
        spectral_nloc[0] = flups_topo_get_nloc(topoSpec, 0);
        spectral_nloc[1] = flups_topo_get_nloc(topoSpec, 1);
        spectral_nloc[2] = flups_topo_get_nloc(topoSpec, 2);
        spectral_end[0]  = spectral_start[0] + spectral_nloc[0];
        spectral_end[1]  = spectral_start[1] + spectral_nloc[1];
        spectral_end[2]  = spectral_start[2] + spectral_nloc[2];

        //--------------------------------------------------------------------------
        /** - warm up */
        //--------------------------------------------------------------------------
        for (int iter = 0; iter < n_warm; ++iter) {
            flups_do_FFT(mysolver, solFLU, FLUPS_FORWARD);
            MPI_Barrier(comm);
            flups_do_FFT(mysolver, solFLU, FLUPS_BACKWARD);
        }

        //--------------------------------------------------------------------------
        /** - Proceed to the solve in place */
        //--------------------------------------------------------------------------
        m_profStart(&prof, "beatme");
        for (int iter = 0; iter < n_iter; iter++) {
            MPI_Barrier(comm);
            m_profStart(&prof, "flups");
            m_profStart(&prof, "flups - forward");
            flups_do_FFT(mysolver, solFLU, FLUPS_FORWARD);
            cudaDeviceSynchronize();
            m_profStop(&prof, "flups - forward");
            m_profStop(&prof, "flups");

            MPI_Barrier(comm);

            m_profStart(&prof, "flups");
            m_profStart(&prof, "flups - backward");
            flups_do_FFT(mysolver, solFLU, FLUPS_BACKWARD);
            cudaDeviceSynchronize();
            m_profStop(&prof, "flups - backward");
            m_profStop(&prof, "flups");
        }
        m_profStop(&prof, "beatme");

        //--------------------------------------------------------------------------
        /** - cleanup */
        //-------------------------------------------------------------------------
        flups_cleanup_solver(mysolver);
        flups_topo_free(topoTmp);
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                flups_free(mybc[id][is]);
            }
        }
        flups_cleanup_backend();
    }

    //--------------------------------------------------------------------------
    /** - Benchmark cuFFTMp */
    //--------------------------------------------------------------------------
#define CUFFT_CHECK(call)                                                                               \
    do {                                                                                                \
        cufftResult result = (call);                                                                    \
        if (result != CUFFT_SUCCESS) {                                                                  \
            int _rank;                                                                                  \
            MPI_Comm_rank(MPI_COMM_WORLD, &_rank);                                                      \
            fprintf(stderr, "[rank %d] cuFFT error: %d at %s:%d\n", _rank, result, __FILE__, __LINE__); \
            MPI_Abort(MPI_COMM_WORLD, 1);                                                               \
        }                                                                                               \
    } while (0)

    {
        //--------------------------------------------------------------------------
        /** - Initialize cuFFTMp */
        //--------------------------------------------------------------------------
        if (rank == 0) printf("Initialization of cuFFTMp\n");

        int nx             = nglob[0];
        int ny             = nglob[1];
        int nz             = nglob[2];
        int n[3]           = {nx, ny, nz};
        int nz_real        = nz;
        int nz_complex     = (nz / 2 + 1);
        int nz_real_padded = 2 * nz_complex;

        // helper stuff
        const int nranks1d = std::sqrt(comm_size);
        assert(nranks1d * nranks1d == comm_size);
        const auto displacement = [](int length, int rank, int size) {
            int ranks_cutoff = length % size;
            return (rank < ranks_cutoff ? rank * (length / size + 1) : ranks_cutoff * (length / size + 1) + (rank - ranks_cutoff) * (length / size));
        };

        // get the disctribution (see cufftMp sample: https://github.com/NVIDIA/CUDALibrarySamples/blob/main/cuFFTMp/samples/r2c_c2r_pencils_cufftMpMakePlan/cufftmp_r2c_c2r_pencils_cufftMpMakePlan.cu)
        long long input_start[3], input_end[3], output_start[3], output_end[3], input_stride[3], output_stride[3];
        for (int i = 0; i < nranks1d; i++) {
            for (int j = 0; j < nranks1d; j++) {
                // skip ranks that don't belong to this pencil
                if (rank != i * nranks1d + j) continue;

                // input
                input_start[0] = displacement(nx, i, nranks1d);
                input_start[1] = displacement(ny, j, nranks1d);
                input_start[2] = 0;

                input_end[0] = displacement(nx, i + 1, nranks1d);
                input_end[1] = displacement(ny, j + 1, nranks1d);
                input_end[2] = nz_real;

                input_stride[0] = nz_real_padded * (input_end[1] - input_start[1]);
                input_stride[1] = nz_real_padded;
                input_stride[2] = 1;

                // output
                output_start[0] = displacement(nx, i, nranks1d);
                output_start[1] = 0;
                output_start[2] = displacement(nz_complex, j, nranks1d);

                output_end[0] = displacement(nx, i + 1, nranks1d);
                output_end[1] = ny;
                output_end[2] = displacement(nz_complex, j + 1, nranks1d);

                output_stride[0] = (output_end[2] - output_start[2]) * (output_end[1] - output_start[1]);
                output_stride[1] = (output_end[2] - output_start[2]);
                output_stride[2] = 1;
            }
        }

        // Create plan
        cufftHandle plan_fwd, plan_bwd;
        size_t      workspace_fwd, workspace_bwd;

        CUFFT_CHECK(cufftCreate(&plan_fwd));
        CUFFT_CHECK(cufftCreate(&plan_bwd));

        CUFFT_CHECK(cufftMpMakePlanDecomposition(
            plan_fwd, 3, n,
            input_start, input_end, input_stride,
            output_start, output_end, output_stride,
            CUFFT_R2C, &comm, CUFFT_COMM_MPI, &workspace_fwd));
        CUFFT_CHECK(cufftMpMakePlanDecomposition(
            plan_bwd, 3, n,
            input_start, input_end, input_stride,
            output_start, output_end, output_stride,
            CUFFT_C2R, &comm, CUFFT_COMM_MPI, &workspace_bwd));

        // Allocate multi-GPU descriptor
        cudaLibXtDesc *desc;
        CUFFT_CHECK(cufftXtMalloc(plan_fwd, &desc, CUFFT_XT_FORMAT_DISTRIBUTED_INPUT));

        //--------------------------------------------------------------------------
        /** - warm up */
        //--------------------------------------------------------------------------
        for (int iter = 0; iter < n_warm; ++iter) {
            CUFFT_CHECK(cufftXtExecDescriptor(plan_fwd, desc, desc, CUFFT_FORWARD));
            CUFFT_CHECK(cufftXtExecDescriptor(plan_bwd, desc, desc, CUFFT_INVERSE));
        }

        //--------------------------------------------------------------------------
        /** - Proceed to the solve in place */
        //--------------------------------------------------------------------------
        m_profStart(&prof, "beatme");
        for (int iter = 0; iter < n_iter; iter++) {
            MPI_Barrier(comm);
            m_profStart(&prof, "cufftMp");
            m_profStart(&prof, "cufftMp - forward");
            CUFFT_CHECK(cufftXtExecDescriptor(plan_fwd, desc, desc, CUFFT_FORWARD));
            m_profStop(&prof, "cufftMp - forward");
            m_profStop(&prof, "cufftMp");

            MPI_Barrier(comm);

            m_profStart(&prof, "cufftMp");
            m_profStart(&prof, "cufftMp - backward");
            CUFFT_CHECK(cufftXtExecDescriptor(plan_bwd, desc, desc, CUFFT_INVERSE));
            m_profStop(&prof, "cufftMp - backward");
            m_profStop(&prof, "cufftMp");
        }
        m_profStop(&prof, "beatme");

        //--------------------------------------------------------------------------
        /** - cleanup */
        //-------------------------------------------------------------------------
        CUFFT_CHECK(cufftXtFree(desc));
        CUFFT_CHECK(cufftDestroy(plan_fwd));
        CUFFT_CHECK(cufftDestroy(plan_bwd));
    }
#undef CUFFT_CHECK

#ifdef HAVE_HEFFTE
    //--------------------------------------------------------------------------
    /** - Benchmark HEFFTE */
    //--------------------------------------------------------------------------
    {
        //--------------------------------------------------------------------------
        /** - Initialize HeFFTE */
        //--------------------------------------------------------------------------
        if (rank == 0) printf("Initialization of HeFFTE\n");

        // create the box using FLUPS' topology
        using heffte_backend        = heffte::backend::cufft;
        heffte::box3d<> const inbox = {
            {topo_start[0], topo_start[1], topo_start[2]},
            {topo_end[0] - 3, topo_end[1] - 1, topo_end[2] - 1}};
        heffte::box3d<> const outbox = {
            {spectral_start[0], spectral_start[1], spectral_start[2]},
            {spectral_end[0] - 1, spectral_end[1] - 1, spectral_end[2] - 1}};

        // set some options
        heffte::plan_options options = heffte::default_options<heffte_backend>();
        options.use_pencils          = true;
        options.use_gpu_aware        = true;
        options.use_reorder          = false;
        options.algorithm            = heffte::reshape_algorithm::p2p_plined;
        heffte::fft3d_r2c<heffte_backend> heffte_fft(inbox, outbox, 0, comm, options);

        // create the data vectors
        std::vector<double> heffte_input(heffte_fft.size_inbox());
        std::iota(heffte_input.begin(), heffte_input.end(), 0);

        heffte::gpu::vector<double>               heffte_gpu_input = heffte::gpu::transfer().load(heffte_input);
        heffte::gpu::vector<std::complex<double>> heffte_gpu_output(heffte_fft.size_outbox());

        //--------------------------------------------------------------------------
        /** - warm up */
        //--------------------------------------------------------------------------
        for (int iter = 0; iter < n_warm; ++iter) {
            heffte_fft.forward(heffte_gpu_input.data(), heffte_gpu_output.data(), heffte::scale::full);
            MPI_Barrier(comm);
            heffte_fft.backward(heffte_gpu_output.data(), heffte_gpu_input.data());
        }

        //--------------------------------------------------------------------------
        /** - Proceed to the solve in place */
        //--------------------------------------------------------------------------
        m_profStart(&prof, "beatme");
        for (int iter = 0; iter < n_iter; iter++) {
            m_profStart(&prof, "heffte");
            m_profStart(&prof, "heffte - forward");
            heffte_fft.forward(heffte_gpu_input.data(), heffte_gpu_output.data(), heffte::scale::full);
            m_profStop(&prof, "heffte - forward");
            m_profStop(&prof, "heffte");

            MPI_Barrier(comm);

            m_profStart(&prof, "heffte");
            m_profStart(&prof, "heffte - backward");
            heffte_fft.backward(heffte_gpu_output.data(), heffte_gpu_input.data());
            m_profStop(&prof, "heffte - backward");
            m_profStop(&prof, "heffte");
        }
        m_profStop(&prof, "beatme");

        //--------------------------------------------------------------------------
        /** - Check result */
        //--------------------------------------------------------------------------
        std::vector<double> heffte_output = heffte::gpu::transfer::unload(heffte_gpu_input);

        double err_loc = 0.0;
        for (int i = 0; i < heffte_input.size(); i++) {
            err_loc = std::max(err_loc, std::abs(heffte_output[i] - heffte_input[i]));
        }

        double err_glob = 0.0;
        MPI_Reduce(&err_loc, &err_glob, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        if (0 == rank) printf("HeFFTE error = %e\n", err_glob);
    }
#endif

    prof.Disp();
    MPI_Finalize();
}
