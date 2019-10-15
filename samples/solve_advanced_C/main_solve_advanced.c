#include <math.h>

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "flups.h"
#include "mpi.h"

// void print_res(double *A, const Topology* topo);


int main(int argc, char *argv[]) {    
    printf("Starting...\n");
    //-------------------------------------------------------------------------
    // Initialize MPI
    //-------------------------------------------------------------------------
    int rank, comm_size;

    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        printf("Invalid number of procs\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
   
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //-----------mak--------------------------------------------------------------
    //Definition of the problem
    //-------------------------------------------------------------------------
    const int     nglob[3] = {64, 64, 64};
    const int     nproc[3] = {1, 1, 2};
    const double  L[3]     = {1., 1., 1.};;

    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};

    const FLUPS_BoundaryType mybc[3][2] = {PER, PER,
                                            PER, PER,
                                            PER, PER};

    if(comm_size!=nproc[0]*nproc[1]*nproc[2]){
        printf("Invalid number of procs\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }


    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------

    // create a real topology

    const FLUPS_Topology *topoIn      = flups_topo_new(0, nglob, nproc, false, NULL, FLUPS_ALIGNMENT);
    const int             nprocOut[3] = {1, 2, 1};

    
    // solver creation and init
    FLUPS_Solver *mysolver = flups_init(topoIn, mybc, h, L);
    flups_set_greenType(mysolver,CHAT_2);
    double* solFLU = flups_setup(mysolver);

    // retrieveing internal info from the solver
    const FLUPS_Topology *topoSpec = flups_get_topo_spectral(mysolver);

    int nmemIn[3];
    int istartIn[3], istartSpec[3];
    int nmemIN[3],nmemSpec[3];

    int memsizeIN = flups_topo_get_memsize(topoIn);
    for (int i = 0; i < 3; i++) {
        nmemIn[i]   = flups_topo_get_nmem(topoIn, i);
        nmemSpec[i] = flups_topo_get_nmem(topoSpec, i);
    }
    flups_topo_get_istartGlob(topoIn,istartIn);
    flups_topo_get_istartGlob(topoSpec,istartSpec);
    
    // //-------------------------------------------------------------------------
    // /** - allocate rhs and solution */
    // //-------------------------------------------------------------------------
    
    // printf("[FLUPS] topo IN glob : %d %d %d \n",topoIn->nglob(0),topoIn->nglob(1),topoIn->nglob(2));
    // printf("[FLUPS] topo IN loc : %d*%d*%d = %d (check: %d %d %d)\n",topoIn->nmem(0),topoIn->nmem(1),topoIn->nmem(2),topoIn->memsize(),topoIn->nloc(0),topoIn->nloc(1),topoIn->nloc(2));
    // printf("[FLUPS] topo OUT glob : %d %d %d \n",topoSpec->nglob(0),topoSpec->nglob(1),topoSpec->nglob(2));
    // printf("[FLUPS] topo OUT loc  : nmem: %d*%d*%d nf:%d (nloc: %d %d %d)  \n",topoSpec->nmem(0),topoSpec->nmem(1),topoSpec->nmem(2),topoSpec->nf(),topoSpec->nloc(0),topoSpec->nloc(1),topoSpec->nloc(2));

    

    printf("I am going to allocate rhs: %d \n",memsizeIN);
    fflush(stdout);
    
 
    double *rhsFLU   = (double *)flups_malloc(sizeof(double) * memsizeIN);
    memset(rhsFLU, 0, sizeof(double ) * memsizeIN);
    // memset(solFLU, 0, sizeof(double ) * FLUmemsizeOUT); 
    
    // printf("istart: %d %d %d =? %d %d %d(P3D)\n",istartIn[0],istartIn[1],istartIn[2],glob_startIN[0],glob_startIN[1],glob_startIN[2]);

    double f = 1; //frequency of the wave
    const double c_2pi = 2.0*M_PI;

    for (int i2 = 0; i2 < flups_topo_get_nloc(topoIn,2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topoIn,1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topoIn,0); i0++) {
                double       x    = (istartIn[0] + i0 + 0.5) * h[0];
                double       y    = (istartIn[1] + i1 + 0.5) * h[1];
                double       z    = (istartIn[2] + i2 + 0.5) * h[2];

                FLUPS_SIZE id;
                id    = flups_locID(0, i0, i1, i2, 0, nmemIn, 1);
                rhsFLU[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
            }
        }
    }
    memcpy(solFLU, rhsFLU, sizeof(double ) * memsizeIN);

    // //-------------------------------------------------------------------------
    // /** - Proceed to the FWD 3D FFT */
    // //-------------------------------------------------------------------------
    
    flups_do_FFT(mysolver,solFLU,FLUPS_FORWARD);

    // //-------------------------------------------------------------------------
    // /** - Multiplication in spectral space */
    // //-------------------------------------------------------------------------
    flups_do_mult(mysolver,solFLU,SRHS);

    // //-------------------------------------------------------------------------
    // /** - Gaussian filtering of the solution */
    // //-------------------------------------------------------------------------

    for (int i2 = 0; i2 < flups_topo_get_nloc(topoIn,2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topoIn,1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topoIn,0); i0++) {
                double       x    = (istartIn[0] + i0 + 0.5) * h[0];
                double       y    = (istartIn[1] + i1 + 0.5) * h[1];
                double       z    = (istartIn[2] + i2 + 0.5) * h[2];

                FLUPS_SIZE id;
                id    = flups_locID(0, i0, i1, i2, 0, nmemIn, 1);
                rhsFLU[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
            }
        }
    }

    // //-------------------------------------------------------------------------
    // /** - Proceed to the BACKWARD 3D FFT */
    // //-------------------------------------------------------------------------
    flups_do_FFT(mysolver,solFLU,FLUPS_BACKWARD);

// #define PRINT_RES 
#ifdef PRINT_RES
    /* normalize*/
    print_res(solFLU,topoIn);
#endif

#ifdef H5_DUMP
    
#endif

    // //-------------------------------------------------------------------------
    // /** - CLEAN */
    // //-------------------------------------------------------------------------
    
    if(rank==0)
        printf("Done ! Now let's clean...\n");
    fflush(stdout);
    
    flups_free(rhsFLU);

    // topoIn->~Topology();
    // topoSpec->~Topology();
    // mysolver->~Solver();

    MPI_Finalize();
}

// void print_res(double *A, const Topology* topo) {
//     const int ax0     = topo->axis();
//     const int ax1     = (ax0 + 1) % 3;
//     const int ax2     = (ax0 + 2) % 3;
//     int nmem[3];
    
//     for(int i=0;i<3;i++){
//         nmem[i] = topo->nmem(i);
//     }

//     int gstart[3];
//     topo->get_istart_glob(gstart);

//     if (topo->isComplex()){
//         for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
//             for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
//                 //local indexes start
//                 const FLUPS_SIZE id = localIndex(ax0, 0, i1, i2, ax0, nmem,2);
//                 for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
//                     printf("(%d %d %d) %lg +i1* %lg\n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id + i0 * 2], A[id + i0 * 2 + 1]);
//                 }
//             }
//         }
//     } else {
//         for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
//             for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
//                 //local indexes start
//                 const FLUPS_SIZE id = localIndex(ax0, 0, i1, i2, ax0, nmem,1);
//                 for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
//                     printf("(%d %d %d) %lg \n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id + i0]);
//                 }
//             }
//         }
//     }
// }
