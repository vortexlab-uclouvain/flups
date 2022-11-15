#include "gtest/gtest.h"

#include <limits>

#include "mpi.h"

#include "h3lpr/macros.hpp"
#include "h3lpr/profiler.hpp"
#include "flups.h"

#include "analytical_field.hpp"

#define m_log(format, ...) m_log_def("test", format, ##__VA_ARGS__)
#define ZERO_TOL 1000.0 * std::numeric_limits<double>::epsilon() 


TEST(PostProcessing, TestDst){
    //-------------------------------------------------------------------------
    // Init MPI stuff
    int      rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    //-------------------------------------------------------------------------
    // Definition of the problem
    //-------------------------------------------------------------------------
    const int nglob[3] = {33, 33, 33};
    const int nproc[3] = {2, 2, 2};
    const double L[3] = {1., 1., 1.};
    const double h[3] = {L[0] / (nglob[0] - 1), L[1] / (nglob[1] - 1), L[2] / (nglob[2] - 1)};  
    int nsolve = 1;

    FLUPS_CenterType    center_type[3] = {NODE_CENTER, NODE_CENTER, NODE_CENTER};
    FLUPS_BoundaryType *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            mybc[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
        }
    }
    
    mybc[0][0][0] = ODD; 
    mybc[0][1][0] = ODD; 
    
    mybc[1][0][0] = EVEN;
    mybc[1][1][0] = ODD;
    
    mybc[2][0][0] = ODD;
    mybc[2][1][0] = EVEN;
    
    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------
    // create a real topology
    FLUPS_Topology *topo      = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);

    // Solver creation and init
    FLUPS_Solver *mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, nullptr);
    flups_set_greenType(mysolver,CHAT_2);
    flups_setup(mysolver, true);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *sol = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    memset(sol, 0, sizeof(double) * flups_topo_get_memsize(topo));
    memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);

    const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;

    AnalyticalField analytics[3];
    for (int dir = 0; dir < 3; ++dir) {
        analytics[dir].SetParam(*mybc[dir][0], *mybc[dir][1]);
    }
    analytics[0].SetFreq(1.);
    analytics[1].SetFreq(2.);
    analytics[2].SetFreq(4.);

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                sol[id]         = 1.0;
            }
        }
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id   = flups_locID(0, i0, i1, i2, 0, ax0, nmem, 1);
                const double x[3] = {(istart[ax0] + i0) * h[ax0],
                                     (istart[ax1] + i1) * h[ax1],
                                     (istart[ax2] + i2) * h[ax2]};
                for (int dir0 = 0; dir0 < 3; dir0++) {
                    const int dir1 = (dir0 + 1) % 3;
                    const int dir2 = (dir0 + 2) % 3;
                    sol[id] *= analytics[dir0].Sol(x[dir0], L[dir0]);
                    rhs[id] += analytics[dir0].Rhs(x[dir0], L[dir0]) * analytics[dir1].Sol(x[dir1], L[dir1]) * analytics[dir2].Sol(x[dir2], L[dir2]);
                }
            }
        }
    }

    //-------------------------------------------------------------------------
    /** Add some spurious values in the data supposed to be corrected        */
    //-------------------------------------------------------------------------
    // Check if my rank contains the first or the last data
    bool have_first_point[3];
    bool have_last_point[3];
    for (int id = 0; id < 3; id++) {
        have_first_point[id] = (istart[id] == 0);
        have_last_point[id]  = (istart[id] + flups_topo_get_nloc(topo, (flups_topo_get_axis(topo) + id) % 3) == nglob[id]);
    }

    //.........................................................................
    // Add spurious data in the x-direction
    // This is an Odd-Odd case. Both first and last data are discarded by fftw and corrected by flups
    if(have_first_point[0]){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                int          first_idx_x  = 0;
                double       first_pos[3] = {(istart[ax0] + first_idx_x) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t first_loc    = flups_locID(0, first_idx_x, i1, i2, 0, 0, nmem, 1);
                rhs[first_loc]            = 1.0;
            }
        }
    }

    if(have_last_point[0]){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                int          last_idx_x   = flups_topo_get_nloc(topo, ax0) - 1;
                double       last_pos[3]  = {(istart[ax0] + last_idx_x) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t last_loc     = flups_locID(0, last_idx_x, i1, i2, 0, 0, nmem, 1);
                rhs[last_loc]             = 1.0;
            }
        }
    }

    //.........................................................................
    // Add spurious data in the y-direction
    // This is an Even-Odd case. Only the last point is corrected by flups
    if(have_last_point[1]){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                int          last_idx_y   = flups_topo_get_nloc(topo, ax1) - 1;
                double       last_pos[3]  = {(istart[ax0] + i0) * h[ax0], (istart[ax1] + last_idx_y) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t last_loc     = flups_locID(0, i0, last_idx_y, i2, 0, 0, nmem, 1);
                rhs[last_loc]             = 10.0;
            }
        }
    }

    //.........................................................................
    // Add spurious data in the z-direction
    // This is an Odd-even case. Only the last point is corrected by flups
    if(have_first_point[2]){
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                int          first_idx_z  = 0;
                double       first_pos[3] = {(istart[ax0] + i0) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + first_idx_z) * h[ax2]};
                const size_t first_loc    = flups_locID(0, i0, i1, first_idx_z, 0, 0, nmem, 1);
                rhs[first_loc]            = 100.0;
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - solve the equation */
    //-------------------------------------------------------------------------
    for (int is = 0; is < nsolve; is++) {
        flups_solve(mysolver, field, rhs, STD);
    }

    //-------------------------------------------------------------------------
    /** - Compute the error */
    //-------------------------------------------------------------------------
    double err = std::numeric_limits<double>::lowest();
    if((0 == flups_topo_get_nloc(topo,ax2) - have_last_point[2]) || (0 == flups_topo_get_nloc(topo,ax1) - have_last_point[1]) || (0 == flups_topo_get_nloc(topo,ax1) - have_last_point[0])){
        m_log("You don't have enough points to fill in all the procs - my error (rank %d) is set to 0 by default \n ", rank);
        err = 0.0;
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                double       pos[3] = {(istart[ax0] + i0) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t id     = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                // Gaussian
                err = std::max(err, abs(field[id] - sol[id]));
            }
        }
    }

    // --------------------------------------------------------------------------
    // Global error 
    double err_glob = 0.0;
    MPI_Allreduce(&err, &err_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    m_log("--------------- Global Error ---------------");
    m_log("The global error == %e", err_glob);

    ASSERT_NEAR(err_glob, 0.0, ZERO_TOL);

    //-------------------------------------------------------------------------
    /** - CLEAN */
    //-------------------------------------------------------------------------
    flups_free(rhs);
    flups_free(sol);
    flups_free(field);
    flups_topo_free(topo);
    flups_cleanup(mysolver);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }
    //-------------------------------------------------------------------------
}


TEST(PostProcessing, TestDft){
    //-------------------------------------------------------------------------
    // Init MPI stuff
    int      rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);

    //-------------------------------------------------------------------------
    // Definition of the problem
    //-------------------------------------------------------------------------
    const int nglob[3] = {33, 33, 33};
    const int nproc[3] = {2, 2, 2};
    const double L[3] = {1., 1., 1.};
    const double h[3] = {L[0] / (nglob[0] - 1), L[1] / (nglob[1] - 1), L[2] / (nglob[2] - 1)};  
    int nsolve = 1;

    FLUPS_CenterType    center_type[3] = {NODE_CENTER, NODE_CENTER, NODE_CENTER};
    FLUPS_BoundaryType *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            mybc[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
            mybc[id][is][0] = PER;
        }
    }

    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------
    // create a real topology
    FLUPS_Topology *topo      = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);

    // Solver creation and init
    FLUPS_Solver *mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, nullptr);
    flups_set_greenType(mysolver,CHAT_2);
    flups_setup(mysolver, true);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *sol = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    memset(sol, 0, sizeof(double) * flups_topo_get_memsize(topo));
    memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);

    const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;

    AnalyticalField analytics[3];
    for (int dir = 0; dir < 3; ++dir) {
        analytics[dir].SetParam(*mybc[dir][0], *mybc[dir][1]);
    }
    analytics[0].SetFreq(1.);
    analytics[1].SetFreq(2.);
    analytics[2].SetFreq(4.);

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                sol[id]         = 1.0;
            }
        }
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id   = flups_locID(0, i0, i1, i2, 0, ax0, nmem, 1);
                const double x[3] = {(istart[ax0] + i0) * h[ax0],
                                     (istart[ax1] + i1) * h[ax1],
                                     (istart[ax2] + i2) * h[ax2]};
                for (int dir0 = 0; dir0 < 3; dir0++) {
                    const int dir1 = (dir0 + 1) % 3;
                    const int dir2 = (dir0 + 2) % 3;
                    sol[id] *= analytics[dir0].Sol(x[dir0], L[dir0]);
                    rhs[id] += analytics[dir0].Rhs(x[dir0], L[dir0]) * analytics[dir1].Sol(x[dir1], L[dir1]) * analytics[dir2].Sol(x[dir2], L[dir2]);
                }
            }
        }
    }

    //-------------------------------------------------------------------------
    /** Add some spurious values in the data supposed to be corrected        */
    //-------------------------------------------------------------------------
    // Check if my rank contains the first or the last data
    bool have_last_point[3];
    for (int id = 0; id < 3; id++) {
        have_last_point[id]  = (istart[id] + flups_topo_get_nloc(topo, (flups_topo_get_axis(topo) + id) % 3) == nglob[id]);
    }

    //.........................................................................
    // Add spurious data in the x-direction
    if(have_last_point[0]){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                int          last_idx_x   = flups_topo_get_nloc(topo, ax0) - 1;
                double       last_pos[3]  = {(istart[ax0] + last_idx_x) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t last_loc     = flups_locID(0, last_idx_x, i1, i2, 0, 0, nmem, 1);
                rhs[last_loc]             = 100.0;
            }
        }
    }

    //.........................................................................
    // Add spurious data in the y-direction
    if(have_last_point[1]){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                int          last_idx_y   = flups_topo_get_nloc(topo, ax1) - 1;
                double       last_pos[3]  = {(istart[ax0] + i0) * h[ax0], (istart[ax1] + last_idx_y) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t last_loc     = flups_locID(0, i0, last_idx_y, i2, 0, 0, nmem, 1);
                rhs[last_loc]             = 100.0;
            }
        }
    }

    //.........................................................................
    // Add spurious data in the z-direction
    if(have_last_point[2]){
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                int          last_idx_z   = flups_topo_get_nloc(topo, ax2) - 1;
                double       last_pos[3]  = {(istart[ax0] + i0) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + last_idx_z) * h[ax2]};
                const size_t last_loc     = flups_locID(0, i0, i1, last_idx_z, 0, 0, nmem, 1);
                rhs[last_loc]             = 100.0;
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - solve the equation */
    //-------------------------------------------------------------------------
    for (int is = 0; is < nsolve; is++) {
        flups_solve(mysolver, field, rhs, STD);
    }

    //-------------------------------------------------------------------------
    /** - Compute the error */
    //-------------------------------------------------------------------------
    double err = std::numeric_limits<double>::lowest();
    if((0 == flups_topo_get_nloc(topo,ax2) - have_last_point[2]) || (0 == flups_topo_get_nloc(topo,ax1) - have_last_point[1]) || (0 == flups_topo_get_nloc(topo,ax1) - have_last_point[0])){
        m_log("You don't have enough points to fill in all the procs - my error (rank %d) is set to 0 by default \n ", rank);
        err = 0.0;
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                double       pos[3] = {(istart[ax0] + i0) * h[ax0], (istart[ax1] + i1) * h[ax1], (istart[ax2] + i2) * h[ax2]};
                const size_t id     = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                // Gaussian
                err = std::max(err, abs(field[id] - sol[id]));
            }
        }
    }

    // --------------------------------------------------------------------------
    // Global error 
    double err_glob = 0.0;
    MPI_Allreduce(&err, &err_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    m_log("--------------- Global Error ---------------");
    m_log("The global error == %e", err_glob);

    ASSERT_NEAR(err_glob, 0.0, ZERO_TOL);

    //-------------------------------------------------------------------------
    /** - CLEAN */
    //-------------------------------------------------------------------------
    flups_free(rhs);
    flups_free(sol);
    flups_free(field);
    flups_topo_free(topo);
    flups_cleanup(mysolver);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }
    //-------------------------------------------------------------------------
}