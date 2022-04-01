#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"
#include "h3lpr/profiler.hpp"
#include "flups.h"

void print_res(double *A, const FLUPS_Topology* topo);

int main(int argc, char *argv[])  {    
    printf("Starting...\n");
    MPI_Init(&argc, &argv);
    int rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    //-------------------------------------------------------------------------
    //Definition of the problem
    //-------------------------------------------------------------------------
    const int     nglob[3] = {9, 9, 9};
    const int     nproc[3] = {1, 1, 1};
    const double  L[3]     = {1., 1., 1.};

    const double h[3] = {L[0] / (nglob[0]-1), L[1] / (nglob[1]-1), L[2] / (nglob[2]-1)};

    const FLUPS_CenterType center_type[3] = {NODE_CENTER, NODE_CENTER, NODE_CENTER};
    // FLUPS_BoundaryType* mybc[3][2] = {{EVEN, EVEN}, {EVEN, EVEN}, {EVEN, EVEN}};

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


    // solver creation and init
    FLUPS_Solver *mysolver = flups_init(topo, mybc, h, L, NOD, center_type);
    flups_set_greenType(mysolver,CHAT_2);
    // flups_set_greenType(mysolver, HEJ_2);
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

    const double cos_cons = 2.0 * M_PI ;
    
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2) ; i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1) ; i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                double       x     = (istart[ax0] + i0) * h[ax0];
                double       y     = (istart[ax1] + i1) * h[ax1];
                double       z     = (istart[ax2] + i2) * h[ax2];
                const size_t id    = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                // second derivative of a cosine
                // rhs[id] = - (cos_cons) * (cos_cons) * cos(cos_cons * x);
                // rhs[id] = - ( cos_cons * cos_cons ) * sin(cos_cons * x); // * sin(cos_cons * y) * sin(cos_cons * z);
                rhs[id] = -(cos_cons * cos_cons) * sin(cos_cons * x);
                // printf("% \t", rhs[id]);
            }
            // printf("\n");
        }
        // printf("\n \n");
    }

    //-------------------------------------------------------------------------
    /** - solve the equations */
    //-------------------------------------------------------------------------
    for (int is = 0; is < 1; is++) {
        flups_solve(mysolver, field, rhs, STD);
    }

    printf("The errrooooor isss \n ");
    printf("ax0 = %d nmeme = %d %d %d-- end = %d %d %d \n", ax0, nmem[0], nmem[1], nmem[2], flups_topo_get_nloc(topo,ax0), flups_topo_get_nloc(topo,ax1), flups_topo_get_nloc(topo,ax2));
    double err = std::numeric_limits<double>::lowest();
    for (int i2 = 1; i2 < flups_topo_get_nloc(topo,ax2) -1   ; i2++) {
        for (int i1 = 1; i1 < flups_topo_get_nloc(topo,ax1) -1; i1++) {
            for (int i0 = 1; i0 < flups_topo_get_nloc(topo,ax0)-1; i0++) {
                double       x     = (istart[ax0] + i0) * h[ax0];
                // double       y     = (istart[1] + i1) * h[1];
                // double       z     = (istart[2] + i2) * h[2];
                const size_t id    = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                
                // Gaussian
                // err = std::max(err, abs(field[id] - cos(cos_cons* x)));
                err = std::max(err, abs(field[id] - sin(cos_cons * x)));
                // err = std::max(err, abs(field[id] - id));
                printf("%e \t",abs(field[id] - sin(cos_cons* x)) );
                // printf("%e \t", field[id] - cos(cos_cons* x))
                // rhs[id] = cos(2.0 * M_PI * 3.0 * x); // sign * c_1o4pi * oosigma3 * sqrt(2.0 / M_PI) * exp(-rho2 * 0.5);
            }
            printf("\n");
        }
        printf("\n \n");
    }

    printf("Error == %e \n", err);
    //SOL
    

    //RHS
    

    // //-------------------------------------------------------------------------
    // /** - CLEAN */
    // //-------------------------------------------------------------------------

    flups_free(rhs);
    flups_free(field);

    flups_topo_free(topo);
    flups_cleanup(mysolver);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }

    MPI_Finalize();
}

void print_res(double *A, const FLUPS_Topology* topo) {
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    int nmem[3];
    
    for(int i=0;i<3;i++){
        nmem[i] = flups_topo_get_nmem(topo,i);
    }

    int gstart[3];
    flups_topo_get_istartGlob(topo,gstart);

    if (flups_topo_get_isComplex(topo)){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                //local indexes start
                const size_t id = flups_locID(ax0, 0, i1, i2, 0, ax0, nmem,2);
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    printf("(%d %d %d) %lg +i1* %lg\n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id + i0 * 2], A[id + i0 * 2 + 1]);
                }
            }
        }
    } else {
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                //local indexes start
                const size_t id = flups_locID(ax0, 0, i1, i2, 0, ax0, nmem,1);
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    printf("(%d %d %d) %lg \n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id + i0]);
                }
            }
        }
    }
}
