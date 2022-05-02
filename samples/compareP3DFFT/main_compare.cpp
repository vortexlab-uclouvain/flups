#include <cmath>
#include <iostream>

#include "Solver.hpp"
#include "p3dfft.h"

void print_all_P3D(double *A,long int nar);
void print_res(double *A, const Topology* topo);

int main(int argc, char *argv[]) {
    
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
        FLUPS_ERROR("The MPI-provided thread behavior does not match");
    }
   
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    //-------------------------------------------------------------------------
    //Definition of the problem
    //-------------------------------------------------------------------------
    const int     nglob[3] = {16, 16, 16};
    const int     nproc[3] = {1, 1, 2}; //nproc[0] has to be 1 //CAUTION FOR THIS: nproc[1]<=nproc[2] !!!
    const double  L[3]     = {1., 1., 1.};

    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};
    
    const FLUPS_CenterType center_type[3] = {CELL_CENTER, CELL_CENTER, CELL_CENTER};

    FLUPS_BoundaryType* mybc[3][2];
    for(int id=0; id<3; id++){
        for(int is=0; is<2; is++){
            mybc[id][is] = (FLUPS_BoundaryType*) flups_malloc(sizeof(int)*1);
            mybc[id][is][0] = PER;
        }
    }

    unsigned char op_f[]="fft", op_b[]="tff";

    const int n_iter = 100;



    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------

    // create a real topology
    int FLUnmemIn[3],FLUnmemOUT[3];
    FLUPS_Topology *topoIn = new FLUPS_Topology(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);
    const int  nprocOut[3] = {1, 2, 1};
    const int  nglobOut[3] = {17, 32, 64};
    
    std::string FLUPSprof = "compare_FLUPS_res" + std::to_string((int)(nglob[0]/L[0])) + "_nrank" + std::to_string(comm_size)+"_nthread" + std::to_string(omp_get_max_threads());
    Profiler* FLUprof = new Profiler(FLUPSprof);

    FLUPS_Solver *mysolver = new FLUPS_Solver(topoIn, mybc, h, L, NOD, center_type, FLUprof);

    mysolver->set_GreenType(CHAT_2);
    double *solFLU = mysolver->setup(true);
    // update the comm and the rank
    comm = flups_topo_get_comm(topoIn);
    MPI_Comm_rank(comm, &rank);

    const Topology *topoSpec    = mysolver->get_innerTopo_spectral();

    for (int i = 0; i < 3; i++) {
        FLUnmemIn[i]  = topoIn->nmem(i);
        FLUnmemOUT[i] = topoSpec->nmem(i);
    }
    
    int istartGloOut[3], istartGlo[3];
    topoIn->get_istart_glob(istartGlo);
    topoSpec->get_istart_glob(istartGloOut);

    int FLUmemsizeIN  = topoIn->memsize();
    int FLUmemsizeOUT = topoSpec->memsize();

    //-------------------------------------------------------------------------
    /** - Initialize P3DFFT */
    //-------------------------------------------------------------------------

    int conf;
    int dims[2] = {nproc[1],nproc[2]};
    int P3Dmemsize[3];
    int istart[3],isize[3],iend[3];
    int fstart[3],fsize[3],fend[3];
    int tstart[3],tsize[3],tend[3];

    if (rank == 0)
        printf("Using processor grid %d x %d with pencils in x\n", dims[0], dims[1]);

    std::string P3DFFTprof = "compare_P3DFFT_res" + std::to_string((int)(nglob[0]/L[0])) + "_nrank" + std::to_string(comm_size)+"_nthread" + std::to_string(omp_get_max_threads());
    Profiler* P3Dprof = new Profiler(P3DFFTprof);

    P3Dprof->create("init", "root");
    P3Dprof->create("FFTandSwitch", "root");

    /* Initialize P3DFFT */
    P3Dprof->start("init");
    Cp3dfft_setup(dims, nglob[0], nglob[1], nglob[2], MPI_Comm_c2f(comm), nglob[0], nglob[1], nglob[2], 1, P3Dmemsize);
    P3Dprof->stop("init");
    
    /* Get dimensions for input array - real numbers, X-pencil shape.
      Note that we are following the Fortran ordering, i.e.
      the dimension  with stride-1 is X. */
    //really? well x seems to be the last dimension...
    conf = 1;
    Cp3dfft_get_dims(istart, iend, isize, conf);
    
    /* Get dimensions for output array - complex numbers, Z-pencil shape.
      Stride-1 dimension could be X or Z, depending on how the library
      was compiled (stride1 option) */
    conf = 2;
    Cp3dfft_get_dims(fstart, fend, fsize, conf);
    
    /* Get what you should allocate in place, in double, per proc. Large enough for padding */
    conf = 3;
    Cp3dfft_get_dims(tstart, tend, tsize, conf);



    // //-------------------------------------------------------------------------
    // /** - allocate rhs and solution */
    // //-------------------------------------------------------------------------
    printf("[FLUPS] topo IN glob : %d %d %d \n",topoIn->nglob(0),topoIn->nglob(1),topoIn->nglob(2));
    printf("[FLUPS] topo IN loc : %d*%d*%d = %d (check: %d %d %d)\n",topoIn->nmem(0),topoIn->nmem(1),topoIn->nmem(2),topoIn->memsize(),topoIn->nloc(0),topoIn->nloc(1),topoIn->nloc(2));
    printf("[FLUPS] topo OUT glob : %d %d %d \n",topoSpec->nglob(0),topoSpec->nglob(1),topoSpec->nglob(2));
    printf("[FLUPS] topo OUT loc  : nmem: %d*%d*%d nf:%d (nloc: %d %d %d)  \n",topoSpec->nmem(0),topoSpec->nmem(1),topoSpec->nmem(2),topoSpec->nf(),topoSpec->nloc(0),topoSpec->nloc(1),topoSpec->nloc(2));

    printf("[P3DFFT] topo 0 loc : %d %d %d  - (is: %d %d %d ie: %d %d %d size: %d %d %d) \n",P3Dmemsize[0],P3Dmemsize[1],P3Dmemsize[2],istart[0],istart[1],istart[2], iend[0],iend[1],iend[2], isize[0],isize[1],isize[2]);
    printf("[P3DFFT] topo 2 loc : %d %d %d  - (is: %d %d %d ie: %d %d %d) \n",fsize[0],fsize[1],fsize[2], fstart[0],fstart[1],fstart[2], fend[0],fend[1],fend[2]);
    printf("[P3DFFT] for alloc: %d %d %d  =? %d %d %d \n",P3Dmemsize[0],P3Dmemsize[1],P3Dmemsize[2], tsize[0],tsize[1],tsize[2]);
    // => P3Dmemsize = tsize


    // int nm = isize[0]*isize[1]*isize[2];
    // printf("isize: %d %d %d - fsize: %d %d %d\n", isize[0],isize[1],isize[2], fsize[0],fsize[1],fsize[2]);
    // if(nm < fsize[0]*fsize[1]*fsize[2]*2){ //switch to complex?
    //     nm = fsize[0]*fsize[1]*fsize[2]*2;
    // }
    int nm = tsize[0]*tsize[1]*tsize[2]; //this is in double
    printf("I am going to allocate FLUPS: %d (out %d R) , P3D: %d (out %d C) \n",FLUmemsizeIN,FLUmemsizeOUT,nm, tsize[0]*tsize[1]*tsize[2]);
    
 
    double *rhsFLU   = (double *)fftw_malloc(sizeof(double) * FLUmemsizeIN);
    // solFLU   allocated by setup
    double *rhsP3D   = (double *)fftw_malloc(sizeof(double) * nm);
    double *solP3D   = (double *)fftw_malloc(sizeof(double) * nm);

    std::memset(rhsFLU, 0, sizeof(double ) * FLUmemsizeIN);
    std::memset(solFLU, 0, sizeof(double ) * FLUmemsizeOUT); 
    std::memset(rhsP3D, 0, sizeof(double ) * nm);
    std::memset(solP3D, 0, sizeof(double ) * nm);
    
    printf("istart: %d %d %d =? %d %d %d(P3D)\n",istartGlo[0],istartGlo[1],istartGlo[2],istart[0],istart[1],istart[2]);

    double f = 1;

    for (int i2 = 0; i2 < topoIn->nloc(2); i2++) {
        for (int i1 = 0; i1 < topoIn->nloc(1); i1++) {
            for (int i0 = 0; i0 < topoIn->nloc(0); i0++) {
                double       x    = (istartGlo[0] + i0 ) * h[0]; //+ 0.5
                double       y    = (istartGlo[1] + i1 ) * h[1]; //+ 0.5
                double       z    = (istartGlo[2] + i2 ) * h[2]; //+ 0.5

                size_t id;
                id         = localIndex(0, i0, i1, i2, 0, FLUnmemIn, 1, 0);
                rhsFLU[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
                
                //p3d does not care of the size you allocate, juste fill the first isize elements
                id =  localIndex(0, i0, i1, i2, 0, isize, 1, 0);
                rhsP3D[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
            }
        }
    }

    // //-------------------------------------------------------------------------
    // /** - Proceed to the solve */
    // //-------------------------------------------------------------------------
    
    int Ntot = fsize[0]*fsize[1];
    Ntot *= fsize[2]*2; //the number of real corresponding to the complex numbers
    
    double factor = 1.0/(nglob[0]* nglob[1]*nglob[2]);

    for (int iter=0;iter<n_iter;iter++){
        
        MPI_Barrier(comm);
        if(rank==0){
            printf("Iter %i\n",iter);
        }
        fflush(stdout);

        // ------------------P3D---------------:

        P3Dprof->start("FFTandSwitch");
        Cp3dfft_ftran_r2c(rhsP3D,solP3D,op_f);
        P3Dprof->stop("FFTandSwitch");

//#define PRINT_RES
#ifdef PRINT_RES
        /* normallize */
        for(int id = 0; id<Ntot; id++){
            solP3D[id] *= factor;
        }   
        print_all_P3D(solP3D,Ntot);
#endif

        // ------------------FLUPS---------------:
        MPI_Barrier(comm);

        std::memcpy(solFLU, rhsFLU, sizeof(double ) * FLUmemsizeIN);
        mysolver->do_FFT(solFLU,FLUPS_FORWARD);

#ifdef PRINT_RES
        /* normalize*/
        for(int id = 0; id<FLUmemsizeOUT; id++){
            solFLU[id] *= factor;
        }
        print_res(solFLU,topoSpec);
#endif
        //Note: if we do several FFT in a raw, the results are wrong with FLUPS for iter>=1. This is because
        // FLUPS would require to do the backward transform before doing a new forward transform (or a reset
        // function, which we choose not to implement).

    }

    // //-------------------------------------------------------------------------
    // /** - get timings */
    // //-------------------------------------------------------------------------
    
    // --- FLUPS -------
    FLUprof->disp("solve");
    delete(FLUprof);

    // --- P3DFFT -------
    double timers[12], gtmean[12], gtmax[12], gtmin[12];

    get_timers(timers);
    MPI_Reduce(&timers,&gtmean,12,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Reduce(&timers,&gtmax ,12,MPI_DOUBLE,MPI_MAX,0,comm);
    MPI_Reduce(&timers,&gtmin ,12,MPI_DOUBLE,MPI_MIN,0,comm);

    for (int i=0;i < 12;i++) {
        gtmean[i] = gtmean[i]/ ((double) comm_size);
    }

    double timeRef      = P3Dprof->get_timeAcc("FFTandSwitch");
    double timeInit     = P3Dprof->get_timeAcc("init");
    string P3DNames[12] = {"All2All(v) 1", "All2All(v) 2", "?", "?", "fft 1", "reorder 1", "fft 2", "reorder+fft 3", "?", "?", "?", "?"};
    int    order[12]    = {5, 1, 6, 7, 2, 8, 3, 4, 9, 10, 11, 12};


    P3Dprof->disp("FFTandSwitch");
    if(rank == 0) {
        // printf("===================================================================================================================================================\n");
        // printf("          TIMER P3DFFT %s\n",P3DFFTprof.c_str());
        // printf("%25s|  %-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\n","-NAME-    ", "-% global-", "-% local-", "-Total time-", "-Self time-", "-time/call-", "-Min time-", "-Max time-","-Mean cnt-","-(MB/s)-");

        string filename = "prof/" + P3DFFTprof + "_time.csv";
        FILE* file      = fopen(filename.c_str(), "w+");
    
        printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", "init", 100*timeInit/timeRef, 100*timeInit/(timeRef+timeInit), timeInit, 0., timeInit/n_iter, 0., 0.,(double) n_iter,0.);
        fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f;%09.2f\n", "init", 100*timeInit/timeRef, 100*timeInit/(timeRef+timeInit), timeInit, 0., timeInit/n_iter, 0., 0.,(double) n_iter,0);
        printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", "FFTandSwitch", 100., 100*timeRef/(timeRef+timeInit), timeRef, 0., timeRef/n_iter, 0., 0., (double) n_iter,0.);
        fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f;%09.2f\n", "FFTandSwitch", 100., 100*timeRef/(timeRef+timeInit), timeRef, 0., timeRef/n_iter, 0., 0., (double) n_iter,0.);
        for(int i=0; i<6; i++){
            int j = order[i]-1;
            printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", ("-- "+P3DNames[j]).c_str(), 100.*gtmean[j]/timeRef, 100*gtmean[j]/timeRef, gtmean[j], 0., gtmean[j]/n_iter, gtmin[j]/n_iter, gtmax[j]/n_iter, (double) n_iter,0.);
            fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f;%09.2f\n", ("-- "+P3DNames[j]).c_str(), 100.*gtmean[j]/timeRef, 100*gtmean[j]/timeRef, gtmean[j], 0., gtmean[j]/n_iter, gtmin[j]/n_iter, gtmax[j]/n_iter, (double) n_iter,0.);
        }
        fclose(file);

        // -- backup of all timers --
        filename = "prof/" + P3DFFTprof + "_backup.csv";
        file            = fopen(filename.c_str(), "w+");
        for(int i=0;i < 12;i++) {
            fprintf(file,"timer[%d] (avg/max/min): %lE %lE %lE\n",i+1,gtmean[i],gtmax[i],gtmin[i]);
        }
        fclose(file);
    }   
    delete(P3Dprof);

    if(rank==0)
        printf("Done ! Now let's clean...\n");
    fflush(stdout);

    fftw_free(solP3D);
    fftw_free(rhsP3D);
    fftw_free(rhsFLU);

    topoIn->~Topology();
    topoSpec->~Topology();
    
    mysolver->~Solver();
    
    for(int id=0; id<3; id++){
        for(int is=0; is<2; is++){
            mybc[id][is] = (FLUPS_BoundaryType*) flups_malloc(sizeof(int)*1);
            flups_free(mybc[id][is]);
        }
    }
    Cp3dfft_clean();

    MPI_Finalize();
}



void print_all_P3D(double *A,long int nar)
{
  int x,y,z,conf,Fstart[3],Fsize[3],Fend[3];
  long int i;

  conf = 2;
  Cp3dfft_get_dims(Fstart,Fend,Fsize,conf);
  /*
  Fsize[0] *= 2;
  Fstart[0] = (Fstart[0]-1)*2;
  */
  for(i=0;i < nar;i+=2)
    if(fabs(A[i]) + fabs(A[i+1]) > 1.25e-4) {
      z = i/(2*Fsize[0]*Fsize[1]);
      y = i/(2*Fsize[0]) - z*Fsize[1];
      x = i/2-z*Fsize[0]*Fsize[1] - y*Fsize[0];
      printf("P3D(%d,%d,%d) %.16lg %.16lg\n",x+Fstart[0],y+Fstart[1],z+Fstart[2],A[i],A[i+1]);
    }
}

void print_res(double *A, const Topology* topo) {
    const int ax0     = topo->axis();
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nf = 2;//topo->nf()
    int nmem[3];
    
    for(int i=0;i<3;i++){
        nmem[i] = topo->nmem(i);
    }

    int gstart[3];
    topo->get_istart_glob(gstart);

    if (topo->isComplex()){
        for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
                //local indexes start
                const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, 2, 0);
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                    
                    if (std::fabs(A[id+i0*2]) + std::fabs(A[id+i0*2 + 1]) > 1.25e-4) {
                        printf("FLU(%d %d %d) %lg +i1* %lg\n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id+i0*2], A[id+i0*2 + 1]);
                    }
                }
            }
        }
    } else {
        for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
                //local indexes start
                const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, 1, 0);
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                    
                    if (std::fabs(A[id+i0]) > 1.25e-4) {
                        printf("FLU(%d %d %d) %lg \n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id+i0]);
                    }
                }
            }
        }
    }
}