#include <cmath>
#include <iostream>

#include "flups.h"
#include "h3lpr/parser.hpp"
#include "h3lpr/profiler.hpp"
#include "heffte.h"

int main(int argc, char *argv[]) {
    //-------------------------------------------------------------------------
    // Initialize MPI
    //-------------------------------------------------------------------------
    int      rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;

    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);

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
        const int topo_nmem[3] = {flups_topo_get_nmem(topoIn, 0), flups_topo_get_nmem(topoIn, 1), flups_topo_get_nmem(topoIn, 2)};
        const int topo_nloc[3] = {flups_topo_get_nloc(topoIn, 0), flups_topo_get_nloc(topoIn, 1), flups_topo_get_nloc(topoIn, 2)};
        topo_end[0]            = topo_start[0] + topo_nloc[0];
        topo_end[1]            = topo_start[1] + topo_nloc[1];
        topo_end[2]            = topo_start[2] + topo_nloc[2];

        // Get some info on the spectral topology
        FLUPS_Topology *topoSpec = flups_get_innerTopo_spectral(mysolver);
        flups_topo_get_istartGlob(topoSpec, spectral_start);
        const int spectral_nloc[3] = {flups_topo_get_nloc(topoSpec, 0), flups_topo_get_nloc(topoSpec, 1), flups_topo_get_nloc(topoSpec, 2)};
        spectral_end[0]            = spectral_start[0] + spectral_nloc[0];
        spectral_end[1]            = spectral_start[1] + spectral_nloc[1];
        spectral_end[2]            = spectral_start[2] + spectral_nloc[2];

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
            m_profStop(&prof, "flups - forward");
            m_profStop(&prof, "flups");

            MPI_Barrier(comm);

            m_profStart(&prof, "flups");
            m_profStart(&prof, "flups - backward");
            flups_do_FFT(mysolver, solFLU, FLUPS_BACKWARD);
            m_profStop(&prof, "flups - backward");
            m_profStop(&prof, "flups");
        }
        m_profStop(&prof, "beatme");

        //--------------------------------------------------------------------------
        /** - cleanup */
        //-------------------------------------------------------------------------
        flups_cleanup(mysolver);
        flups_topo_free(topoTmp);
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                flups_free(mybc[id][is]);
            }
        }
        flups_cleanup_backend();
    }

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

    prof.Disp();
    MPI_Finalize();
}
