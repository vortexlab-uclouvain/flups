#include <cmath>
#include <iostream>

#include "Solver.hpp"
#include "accfft.h"
#include "h3lpr/profiler.hpp"
#include "h3lpr/parser.hpp"
#include "flups.h"

using H3LPR::Profiler;
using H3LPR::Parser;

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
    const auto    arg_nglob = parser.GetValues<int, 3>("--nglob", "the global resolution, will be used for both ACCFFT and FLUPS", {64,64,64});
    const auto    arg_nproc = parser.GetValues<int, 3>("--nproc", "the proc distribution, for FLUPS only", {1, 1, 1});
    const auto    arg_dom   = parser.GetValues<double, 3>("--dom", "the size of the domain, must be compatible with nglob", {1.0, 1.0, 1.0});
    const int     n_iter    = parser.GetValue<int>("--niter", "the number of iterations to perform", 20);
    const int     n_warm    = parser.GetValue<int>("--warm", "the number of iterations to perform when warming up", 1);
    parser.Finalize();

    //--------------------------------------------------------------------------
    // Definition of the problem
    //--------------------------------------------------------------------------
    const int nglob[3] = {arg_nglob[0], arg_nglob[1], arg_nglob[2]};
    const int nproc[3] = {arg_nproc[0], arg_nproc[1], arg_nproc[2]};
    const double L[3]     = {arg_dom[0], arg_dom[1], arg_dom[2]};
    //const double L[3] = {1.0, 1.0, 1.0};

    // get the grid spacing
    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};
    FLUPS_CHECK(h[0] == h[1] && h[1] == h[2], "The grid spacing must be the same");

    // get the PER PER PER BC everywhere
    const FLUPS_CenterType center_type[3] = {CELL_CENTER, CELL_CENTER, CELL_CENTER};
    FLUPS_BoundaryType    *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            mybc[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
            mybc[id][is][0] = PER;
        }
    }

    if (rank == 0) {
        printf("form the command line: %d %d %d unknowns on %d %d %d proc\n", nglob[0], nglob[1], nglob[2], nproc[0], nproc[1], nproc[2]);
    }
    if (nproc[0] != 1 && rank == 0) {
        printf("--------------------------------------------------------------\n");
        printf("                      WARNING\n");
        printf("ACCFFT assumes a distribution already aligned in the first dimension they will proceed\n");
        printf(" To be fair with FLUPS, you should impose that nproc[0] = 1\n");
        printf("--------------------------------------------------------------\n");
    }

    //--------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //--------------------------------------------------------------------------
    FLUPS_INFO("Initialization of FLUPS");

    // create a real topology
    FLUPS_Topology topoIn(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);
    FLUPS_Solver  *mysolver = new FLUPS_Solver(&topoIn, mybc, h, L, NOD, center_type, nullptr);

    // set the CHAT2 green type (even if it's not used)
    mysolver->set_GreenType(CHAT_2);
    double *solFLU = mysolver->setup(true);

    //..........................................................................
    // set some straightforward data
    int topo_nmem[3] = {topoIn.nmem(0), topoIn.nmem(1), topoIn.nmem(2)};

    // set a simple expression
    double val = 0.0;
    for (int i2 = 0; i2 < topoIn.nloc(2); i2++) {
        for (int i1 = 0; i1 < topoIn.nloc(1); i1++) {
            for (int i0 = 0; i0 < topoIn.nloc(0); i0++) {
                double x   = 2.0 * M_PI / nglob[0] * (i0 + topoIn.cmpt_start_id(0));
                double y   = 2.0 * M_PI / nglob[1] * (i1 + topoIn.cmpt_start_id(1));
                double z   = 2.0 * M_PI / nglob[2] * (i2 + topoIn.cmpt_start_id(2));
                size_t id  = localIndex(0, i0, i1, i2, 0, topo_nmem, 1, 0);
                solFLU[id] = sin(x) + sin(y) + sin(z);
            }
        }
    }
    //--------------------------------------------------------------------------
    /** - Initialize ACCFFT */
    //--------------------------------------------------------------------------
    // from https://github.com/amirgholami/accfft/blob/master/steps/step2/step2.cpp
    // ACCFFT starts with a pencil in Z distribution, always!

    /* Create Cartesian Communicator */
    int      c_dims[2];
    MPI_Comm c_comm;
    accfft_create_comm(MPI_COMM_WORLD, c_dims, &c_comm);

    // let ACCFFT decide on the topology choice, pencil in Z, as always
    int    isize[3], osize[3], istart[3], ostart[3];

    int n_acc[3] = {nglob[0],nglob[1],nglob[2]};
    size_t alloc_max = accfft_local_size_dft_r2c(n_acc, isize, istart, osize, ostart, c_comm);

    double *data_acc = (double *)accfft_alloc(alloc_max);

    // init the ACCFFT with 1 thread
    accfft_init(1);
    // get the plan
    accfft_plan *plan = accfft_plan_dft_3d_r2c(n_acc, data_acc, data_acc, c_comm, ACCFFT_MEASURE);

    // setup the data
    int n2_ = (nglob[2] / 2 + 1) * 2;
    for (int i = 0; i < isize[0]; i++) {
        for (int j = 0; j < isize[1]; j++) {
            for (int k = 0; k < isize[2]; k++) {
                double x      = 2.0 * M_PI / nglob[0] * (i + istart[0]);
                double y      = 2.0 * M_PI / nglob[1] * (j + istart[1]);
                double z      = 2.0 * M_PI / nglob[2] * k;
                size_t ptr    = i * isize[1] * n2_ + j * n2_ + k;
                data_acc[ptr] = sin(x) + sin(y) + sin(z);
            }
        }
    }

    //--------------------------------------------------------------------------
    /** - warm up everybody */
    //--------------------------------------------------------------------------
    for (int iter = 0; iter < n_warm; ++iter) {
        // warm up accfft
        accfft_execute_r2c(plan, data_acc, (Complex *)data_acc);
        accfft_execute_c2r(plan, (Complex *)data_acc, data_acc);

        // warm up flups
        mysolver->do_FFT(solFLU, FLUPS_FORWARD);
        mysolver->do_FFT(solFLU, FLUPS_BACKWARD);
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /** - Proceed to the solve in place */
    //--------------------------------------------------------------------------
    std::string prof_name = "beatme_nglob" + std::to_string(nglob[0]) + "_nrank" + std::to_string(comm_size);
    Profiler    prof(prof_name);

    m_profStart(&prof, "beatme");
    for (int iter = 0; iter < n_iter; iter++) {
        if (rank == 0) {
            printf("doing solve number %d/%d\n", iter, n_iter);
            printf("------- ACCFFT \n");
        }
        fflush(stdout);

        //......................................................................
        // ACCFFT
        //......................................................................
        MPI_Barrier(comm);
        m_profStart(&prof, "ACCFFT");
        m_profStart(&prof, "ACCFFT - forward");
        accfft_execute_r2c(plan, data_acc, (Complex *)data_acc);
        m_profStop(&prof, "ACCFFT - forward");

        MPI_Barrier(comm);
        m_profStart(&prof, "ACCFFT - backward");
        accfft_execute_c2r(plan, (Complex *)data_acc, data_acc);
        m_profStop(&prof, "ACCFFT - backward");
        m_profStop(&prof, "ACCFFT");

        //......................................................................
        // flups
        //......................................................................
        if (rank == 0) {
            printf("------- flups \n");
        }
        fflush(stdout);
        MPI_Barrier(comm);
        m_profStart(&prof, "flups");
        m_profStart(&prof, "flups - forward");
        mysolver->do_FFT(solFLU, FLUPS_FORWARD);
        m_profStop(&prof, "flups - forward");

        MPI_Barrier(comm);
        m_profStart(&prof, "flups - backward");
        mysolver->do_FFT(solFLU, FLUPS_BACKWARD);
        m_profStop(&prof, "flups - backward");
        m_profStop(&prof, "flups");
    }
    m_profStop(&prof, "beatme");

    //--------------------------------------------------------------------------
    /** - get timings */
    //-------------------------------------------------------------------------
    prof.Disp();

    //--------------------------------------------------------------------------
    /** - cleanup */
    //-------------------------------------------------------------------------
    // ACCFFT
    accfft_free(data_acc);
    accfft_destroy_plan(plan);
    accfft_cleanup();

    // --- FLUPS -------

    // force the call to destructor of the solver to cleanup the comm patterns
    delete (mysolver);

    // free the bcs
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            // mybc[id][is] = (FLUPS_BoundaryType*) flups_malloc(sizeof(int)*1);
            flups_free(mybc[id][is]);
        }
    }

    MPI_Finalize();
}

// end-of-file
