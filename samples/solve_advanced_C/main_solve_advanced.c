#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"
#include "flups.h"

void print_res(double *A, const FLUPS_Topology* topo);

int main(int argc, char *argv[]) {    
    printf("Starting...\n");
    //-------------------------------------------------------------------------
    // Initialize MPI
    //-------------------------------------------------------------------------
    int rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;

    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        printf("Invalid number of procs\n");
        MPI_Abort(comm,1);
    }
   
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

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
        MPI_Abort(comm,1);
    }


    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------

    // create a real topology

    FLUPS_Topology *topoIn      = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT,comm);
    const int             nprocOut[3] = {1, 2, 1};

    
    // solver creation and init
    FLUPS_Solver *mysolver = flups_init(topoIn, mybc, h, L, 1);
    flups_set_greenType(mysolver,CHAT_2);
    double* solFLU = flups_setup(mysolver, true);

    // recompute the communicator and the rank
    comm = flups_topo_get_comm(topoIn);
    MPI_Comm_rank(comm,&rank);

    // retrieveing internal info from the solver
    const FLUPS_Topology *topoSpec = flups_get_innerTopo_spectral(mysolver);

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
    
    printf("[FLUPS] topo IN glob : %d %d %d \n",flups_topo_get_nglob(topoIn,0),flups_topo_get_nglob(topoIn,1),flups_topo_get_nglob(topoIn,2));
    printf("[FLUPS] topo IN loc : %d*%d*%d = %d (check: %d %d %d)\n",flups_topo_get_nmem(topoIn,0),flups_topo_get_nmem(topoIn,1),flups_topo_get_nmem(topoIn,2),flups_topo_get_memsize(topoIn),flups_topo_get_nglob(topoIn,0),flups_topo_get_nglob(topoIn,1),flups_topo_get_nglob(topoIn,2));
    printf("[FLUPS] topo OUT glob : %d %d %d \n",flups_topo_get_nglob(topoSpec,0),flups_topo_get_nglob(topoSpec,1),flups_topo_get_nglob(topoSpec,2));
    printf("[FLUPS] topo OUT loc  : nmem: %d*%d*%d complex?:%d (nloc: %d %d %d)  \n",flups_topo_get_nmem(topoSpec,0),flups_topo_get_nmem(topoSpec,1),flups_topo_get_nmem(topoSpec,2),flups_topo_get_isComplex(topoSpec),flups_topo_get_nloc(topoSpec,0),flups_topo_get_nloc(topoSpec,1),flups_topo_get_nloc(topoSpec,2));

    printf("I am going to allocate rhs: %d \n",memsizeIN);
    fflush(stdout);
    
 
    double *rhsFLU   = (double *)flups_malloc(sizeof(double) * memsizeIN);
    memset(rhsFLU, 0, sizeof(double ) * memsizeIN);
    // memset(solFLU, 0, sizeof(double ) * FLUmemsizeOUT); 
    
    double f = 1; //frequency of the wave
    const double c_2pi = 2.0*3.141592653589793;

    for (int i2 = 0; i2 < flups_topo_get_nloc(topoIn,2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topoIn,1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topoIn,0); i0++) {
                double       x    = (istartIn[0] + i0 + 0.5) * h[0];
                double       y    = (istartIn[1] + i1 + 0.5) * h[1];
                double       z    = (istartIn[2] + i2 + 0.5) * h[2];

                size_t id;
                id    = flups_locID(0, i0, i1, i2, 0, 0, nmemIn, 1);
                rhsFLU[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
            }
        }
    }
    memcpy(solFLU, rhsFLU, sizeof(double ) * memsizeIN);

    // //-------------------------------------------------------------------------
    // /** - Skip the copy of data
    // //-------------------------------------------------------------------------
    // By using the advanced features of the API, we skip the copy which is done
    // in flups_solve: data are copied from the location the user allocated himself,
    // to an inner memory reserved by flups, and ensuring proper alignment.
    // In this example, we have directly filled that "inner memory" (solFlu) with 
    // the initial condition
    
    // flups_do_copy(mysolver,topoIn,solFLU,FLUPS_FORWARD);

    // //-------------------------------------------------------------------------
    // /** - Proceed to the FWD 3D FFT */
    // //-------------------------------------------------------------------------
    
    flups_do_FFT(mysolver,solFLU,FLUPS_FORWARD);

    // //-------------------------------------------------------------------------
    // /** - Multiplication in spectral space */
    // //-------------------------------------------------------------------------
    flups_do_mult(mysolver,solFLU,RHS);

    // //-------------------------------------------------------------------------
    // /** - Gaussian filtering of the solution */
    // //-------------------------------------------------------------------------

    double kfact[3], koffset[3], symstart[3];
    flups_get_spectralInfo(mysolver,kfact,koffset,symstart);

    printf("%lf %lf %lf  %lf %lf %lf  %lf %lf %lf  \n",kfact[0],kfact[1],kfact[2],koffset[0],koffset[1],koffset[2],symstart[0],symstart[1],symstart[2]);

    const int ax0     = flups_topo_get_axis(topoSpec);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nf      = 2; //topoSpec is full spectral, there are two doubles per element  
    
    // int nmemSpec[3];
    // for(int i=0;i<3;i++){
    //     nmemSpec[i] = flups_topo_get_nmem(topoSpec,i);
    // }
        
    for (int i2 = 0; i2 < flups_topo_get_nloc(topoSpec,ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topoSpec,ax1); i1++) {
            //local indexes start
            const size_t id = flups_locID(ax0, 0, i1, i2, 0, ax0, nmemSpec,nf);
            for (int i0 = 0; i0 < flups_topo_get_nloc(topoSpec,ax0); i0++) {
                int is[3];
                flups_symID(ax0, i0, i1, i2, istartSpec, symstart, 0, is);

                // (symmetrized) wave number
                const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];

                const double ksqr   = k0 * k0 + k1 * k1 + k2 * k2;

                solFLU[id + i0 * nf] *= exp(-ksqr); //REAL part
                solFLU[id + i0 * nf + 1] *= exp(-ksqr); //COMPLEX part
            }
        }
    }

    flups_hdf5_dump(topoSpec,"spectrum",solFLU);

    // //-------------------------------------------------------------------------
    // /** - Proceed to the BACKWARD 3D FFT */
    // //-------------------------------------------------------------------------
    flups_do_FFT(mysolver,solFLU,FLUPS_BACKWARD);

    // //-------------------------------------------------------------------------
    // /** - Skip the copy back of data
    // //-------------------------------------------------------------------------

    // flups_do_copy(mysolver,topoIn,solFLU,FLUPS_BACKWARD);

    // //-------------------------------------------------------------------------
    // /** - Export results */
    // //-------------------------------------------------------------------------

// #define PRINT_RES
#ifdef PRINT_RES
    /* normalize*/
    print_res(solFLU,topoIn);
#endif

    flups_hdf5_dump(topoIn,"result",solFLU);

    // //-------------------------------------------------------------------------
    // /** - CLEAN */
    // //-------------------------------------------------------------------------
    
    if(rank==0)
        printf("Done ! Now let's clean...\n");
    fflush(stdout);
    
    flups_free(rhsFLU);

    flups_topo_free(topoIn);
    flups_topo_free(topoSpec);
    flups_cleanup(mysolver);
    
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
