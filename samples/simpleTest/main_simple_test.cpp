/**
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"
#include "h3lpr/profiler.hpp"
#include "h3lpr/macros.hpp"
#include "flups.h"

#define m_log(format, ...) m_log_def("test", format, ##__VA_ARGS__)

void print_res(double *A, const FLUPS_Topology* topo, std::string filename);

int main(int argc, char *argv[])  {    
    
    MPI_Init(&argc, &argv);
    int rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
    m_log("Starting...");

    //-------------------------------------------------------------------------
    // Definition of the problem
    //-------------------------------------------------------------------------
    bool is_node = true;
    const int nglob[3] = {9, 9, 9};
    const int nproc[3] = {1, 1, comm_size};
    const double L[3] = {1., 1., 1.};
    const double h[3] = {L[0] / (nglob[0] - is_node), L[1] / (nglob[1] - is_node), L[2] / (nglob[2] - is_node)};  
    int nsolve = 1;

    FLUPS_CenterType center_type[3]; 
    FLUPS_BoundaryType *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        center_type[id] = is_node ? NODE_CENTER : CELL_CENTER;
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

    // Profiler creation 
    FLUPS_Profiler* prof = flups_profiler_new();

    // Solver creation and init
    FLUPS_Solver *mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, prof);
    flups_set_greenType(mysolver,CHAT_2);
    flups_setup(mysolver, true);

    // recompute the communicator and the rank
    comm = flups_topo_get_comm(topo);
    MPI_Comm_rank(comm,&rank);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);

    const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;

    const double k = 2.0 * M_PI ;
    auto fun = [k](double pos[3]) -> double { return cos(k * pos[0]); };
    auto d2fundx = [k](double pos[3]) -> double { return (-(k*k)*cos(k * pos[0])); };

    const double shift = is_node? 0.0 : 0.5;
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2) ; i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1) ; i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                double pos[3] ={(istart[ax0] + i0 + shift) * h[ax0], (istart[ax1] + i1 + shift) * h[ax1], (istart[ax2] + i2 + shift) * h[ax2]};
                const size_t id    = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                rhs[id] = d2fundx(pos);
                // rhs[id] = (istart[ax0] + i0) + nmem[0]*((istart[ax1] + i1) + nglob[1]*(i2 + istart[ax2])); //d2fundx(pos);
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - solve the equation */
    //-------------------------------------------------------------------------
    for (int is = 0; is < nsolve; is++) {
        m_log("Iteration %d/%d", is, nsolve);
        flups_solve(mysolver, field, rhs, STD);
    }


    //-------------------------------------------------------------------------
    /** - Compute the error */
    //-------------------------------------------------------------------------
    double err = std::numeric_limits<double>::lowest();
    bool is_last_topo[3]; 
    for(int id = 0; id < 3; id++){
        is_last_topo[id] = (istart[id] + flups_topo_get_nloc(topo,(flups_topo_get_axis(topo) + id)%3) == nglob[id]);
    }

    if((0 == flups_topo_get_nloc(topo,ax2) - is_last_topo[2]) || (0 == flups_topo_get_nloc(topo,ax1) - is_last_topo[1]) || (0 == flups_topo_get_nloc(topo,ax1) - is_last_topo[0])){
        printf("[test - %d] You don't have enough points to fill in all the procs - my error (rank %d) is set to 0 by default \n ", rank, rank);
        err = 0.0;
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2) - is_last_topo[2] ; i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1) - is_last_topo[1] ; i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0) - is_last_topo[0]; i0++) {
                double pos[3] ={(istart[ax0] + i0 + shift) * h[ax0], (istart[ax1] + i1 + shift) * h[ax1], (istart[ax2] + i2 + shift) * h[ax2]};
                const size_t id    = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                // Gaussian
                err = std::max(err, abs(field[id] - fun(pos)));
            }
        }
    }

    //m_log("Printing results");
    //printf("[test %d] Error computation - end %d %d %d ", rank,flups_topo_get_nloc(topo,ax0) - is_last_topo[0], flups_topo_get_nloc(topo,ax1) - is_last_topo[1], flups_topo_get_nloc(topo,ax2) - is_last_topo[2] );
    //std::string arg_name = (argc == 2) ? argv[1] : "default";
    //print_res(field, topo, "data/sol_" + arg_name) ;

    // --------------------------------------------------------------------------
    // Local error 
    printf("[test - %d] my error is %e \n", rank, err);

    // --------------------------------------------------------------------------
    // Global error 
    double err_glob = 0.0;
    MPI_Allreduce(&err, &err_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    m_log("--------------- Global Error ---------------");
    m_log("The global error == %e", err_glob);
    
    flups_profiler_disp(prof);

    //-------------------------------------------------------------------------
    /** - CLEAN */
    //-------------------------------------------------------------------------
    flups_free(rhs);
    flups_free(field);
    flups_topo_free(topo);
    flups_cleanup(mysolver);
    flups_profiler_free(prof);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }

    MPI_Finalize();
}

void print_res(double *A, const FLUPS_Topology* topo, std::string filename) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE *fptr; 
    std::string name = filename + "_" + std::to_string(rank) + ".txt";
    fptr = fopen(name.c_str(), "w");

    if(fptr == NULL){ printf("Error!");   exit(1); }

    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    int nmem[3] = {flups_topo_get_nmem(topo,0), flups_topo_get_nmem(topo,1), flups_topo_get_nmem(topo,2)};
    
    // int gstart[3];
    // flups_topo_get_istartGlob(topo,gstart);

    if (flups_topo_get_isComplex(topo)){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                //local indexes start
                const size_t id = flups_locID(ax0, 0, i1, i2, 0, ax0, nmem,2);
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    // printf("(%d %d %d) %lg +i1 * %lg\n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id + i0 * 2], A[id + i0 * 2 + 1]);
                    fprintf(fptr, "\t (%e, i1 * %e)", A[id + i0 * 2], A[id + i0 * 2 + 1]);
                }
                fprintf(fptr, "\n");
            }
            fprintf(fptr, "\n \n");
        }
    
    } else {
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                //local indexes start
                const size_t id = flups_locID(ax0, 0, i1, i2, 0, ax0, nmem,1);
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    // printf("(%d %d %d) %lg \n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id + i0]);
                    fprintf(fptr, "\t %e", A[id + i0]);
                }
                fprintf(fptr, "\n");
            }
            fprintf(fptr, "\n \n");
        }
    }
    fclose(fptr);
}
