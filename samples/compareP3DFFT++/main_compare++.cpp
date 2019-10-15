#include <cmath>
#include <iostream>

#include "Solver.hpp"
#include "p3dfft.h"

void print_res(p3dfft::complex_double *A, int *sdims, int *gstart);
void print_res(double *A, const Topology* topo);

int main(int argc, char *argv[]) {
    //-------------------------------------------------------------------------
    // Default values
    //-------------------------------------------------------------------------
    int n_iter = 10;
    int res[3] = {16,16,16};
    int nproc2D[2] = {1,2};

    //-------------------------------------------------------------------------
    // Parse arguments
    //-------------------------------------------------------------------------
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            printf(" --nprocs, -np Nj Nk :          the number of MPI processes for a pencil decomposition in x (!! Nj<=Np !!)\n");
            printf(" --resolution, -res Rx Ry Rz :  Rx,Ry,Rz is the total number of points in each direction \n");
            printf(" --niter, -ni Ni :              Ni is the number of times we call the same 3D FFT (for statistics on the profiler) \n");
            return 0;
        } else if ((arg == "-np") || (arg == "--nprocs")) {
            for (int j = 0; j<2;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    nproc2D[j] = atoi(argv[i+j+1]); 
                    if(nproc2D[j]<1){
                        fprintf(stderr, "nprocs must be >0\n");
                        return 1;
                    }
                } else { //Missing argument
                    fprintf(stderr, "missing argument in --nprocs\n");
                    return 1;
                }  
            }
            i+=2;
        }  else if ((arg == "-res")|| (arg== "--resolution") ) {
            for (int j = 0; j<3;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    res[j] = atoi(argv[i+j+1]); 
                    if(res[j]<=0.0){
                        fprintf(stderr, "res must be >0\n");
                        return 1;
                    }
                } else { //Missing argument
                    fprintf(stderr, "missing argument in -res\n");
                    return 1;
                }  
            }
            i+=3;
        }  else if ((arg == "-ni")|| (arg== "--niter") ) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                n_iter = atoi(argv[i+1]); 
                if(n_iter<1){
                    fprintf(stderr, "niter must be >0\n");
                    return 1;
                }
            } else { //Missing argument
                fprintf(stderr, "missing --niter\n");
                return 1;
            }  
            i++;
        // } else if ((arg == "-bc")|| (arg== "--boundary-conditions") ) {
        //     for (int j = 0; j<6;j++){
        //         if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
        //             bcdef[j/2][j%2] = (BoundaryType) atoi(argv[i+j+1]); 
        //         } else { //Missing argument
        //             fprintf(stderr, "missing argument in --boundary-conditions\n");
        //             return 1;
        //         }  
        //     }
        //     i+=6;
        }
    }
    
    //-------------------------------------------------------------------------
    // Initialize MPI
    //-------------------------------------------------------------------------
    int rank, comm_size;

    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        FLUPS_ERROR("The MPI-provided thread behavior does not match", LOCATION);
    }
   
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //-------------------------------------------------------------------------
    //Definition of the problem
    //-------------------------------------------------------------------------
    const int     nglob[3] = {res[0], res[1], res[2]};
    const int     nproc[3] = {1, nproc2D[0], nproc2D[1]}; //nproc[0] has to be 1 //CAUTION FOR THIS: nproc[1]<=nproc[2] !!!
    const double  L[3]     = {1., 1., 1.};;

    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};

    const BoundaryType mybc[3][2] = {PER, PER,
                                    PER, PER,
                                    PER, PER};

    if(comm_size!=nproc[0]*nproc[1]*nproc[2])
        FLUPS_ERROR("Invalid number of procs",LOCATION);


    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------

    // create a real topology
    int FLUnmemIn[3],FLUnmemOUT[3];
    const Topology *topoIn    = new Topology(0, nglob, nproc, false, NULL,FLUPS_ALIGNMENT);
    const int  nprocOut[3] = {1, 2, 1};
    const int  nglobOut[3] = {17, 32, 64};
    
    // prepare profiling
    std::string FLUPSprof = "compare_FLUPS_res" + std::to_string((int)(nglob[0]/L[0])) + "_nrank" + std::to_string(comm_size)+"_nthread" + std::to_string(omp_get_max_threads());
    Profiler* FLUprof = new Profiler(FLUPSprof);

    // solver creation and init
    Solver *mysolver = new Solver(topoIn, mybc, h, L, FLUprof);
    mysolver->set_GreenType(CHAT_2);
    double* solFLU = mysolver->setup(); //already allocated

    // retrieveing internal info from the solver
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
    int  proc_order[3] = {0, 1, 2};
    
    if (rank == 0)
        printf("Using processor grid %d x %d with pencils in x\n", dims[0], dims[1]);

    std::string P3DFFTprof = "compare_P3DFFT_res" + std::to_string((int)(nglob[0]/L[0])) + "_nrank" + std::to_string(comm_size)+"_nthread" + std::to_string(omp_get_max_threads());
    Profiler* P3Dprof = new Profiler(P3DFFTprof);

    P3Dprof->create("init", "root");
    P3Dprof->create("FFTandSwitch", "root");

    /* Initialize P3DFFT */
    P3Dprof->start("init");
    if(rank==0)
        printf("Initilizing P3D...\n");    fflush(stdout);
    p3dfft::setup();
    if(rank==0)
        printf("...set up.\n");    fflush(stdout);

    // Define the transform types for these two transforms
    int type_ids1[3] = {p3dfft::R2CFFT_D,p3dfft::CFFT_FORWARD_D ,p3dfft::CFFT_FORWARD_D };
    int type_ids2[3] = {p3dfft::C2RFFT_D,p3dfft::CFFT_BACKWARD_D,p3dfft::CFFT_BACKWARD_D};
    p3dfft::trans_type3D type_rcc(type_ids1); 
    p3dfft::trans_type3D type_ccr(type_ids2);

    // input grid
    int  gdimsIN[3]      = {topoIn->nglob(0), topoIn->nglob(1), topoIn->nglob(2)};
    int  mem_orderIN[3]  = {0, 1, 2};
    int  nprocIN[3]      = {nproc[0], nproc[1], nproc[2]};
    p3dfft::grid gridIN(gdimsIN,-1,nprocIN,proc_order,mem_orderIN,MPI_COMM_WORLD); 

    if(rank==0)
        printf("...input grid.\n");    fflush(stdout);

    // ouput grid
    int gdimsOUT[3];
    for(int i=0; i < 3;i++) 
        gdimsOUT[i] = gdimsIN[i];
    gdimsOUT[0] = gdimsOUT[0]/2+1;
    int mem_orderOUT[3]  = {2, 1, 0}; //blindly mimicking samples, anything else produces a segfault anyway...
    int nprocOUT[3]      = {dims[0],dims[1],1};
    p3dfft::grid gridOUT(gdimsOUT,0,nprocOUT,proc_order,mem_orderOUT,MPI_COMM_WORLD); 

    if(rank==0)
        printf("...output grid.\n");    fflush(stdout);

    // Set up 3D transforms, including stages and plans, for forward trans.
    p3dfft::transform3D<double,p3dfft::complex_double> trans_f(gridIN,gridOUT,&type_rcc);
    // Set up 3D transforms, including stages and plans, for backward trans.
    p3dfft::transform3D<p3dfft::complex_double,double> trans_b(gridOUT,gridIN,&type_ccr);
    if(rank==0)
        printf("...transforms.\n");    fflush(stdout);

    P3Dprof->stop("init");
    if(rank==0)
        printf("Done with P3D init.\n");    fflush(stdout);

    p3dfft::timers.init();


    // Find local dimensions in storage order, and also the starting position of the local array in the global array
  
    //input grid - SIZE IN DOUBLEs
    int P3DnlocIN[3],glob_startIN[3];
    for(int i=0;i<3;i++) {
        P3DnlocIN[mem_orderIN[i]] = gridIN.ldims[i];
        glob_startIN[mem_orderIN[i]] = gridIN.glob_start[i];
    }
    int P3DmemsizeIN = P3DnlocIN[0]*P3DnlocIN[1]*P3DnlocIN[2];

    //output grid - SIZE IN COMPLEXes
    int P3DnlocOUT[3],glob_startOUT[3];
    for(int i=0;i<3;i++) {
        P3DnlocOUT[mem_orderOUT[i]] = gridOUT.ldims[i];
        glob_startOUT[mem_orderOUT[i]] = gridOUT.glob_start[i];
    }
    int P3DmemsizeOUT = P3DnlocOUT[0]*P3DnlocOUT[1]*P3DnlocOUT[2];

    

    // //-------------------------------------------------------------------------
    // /** - allocate rhs and solution */
    // //-------------------------------------------------------------------------
    
    printf("[FLUPS] topo IN glob : %d %d %d \n",topoIn->nglob(0),topoIn->nglob(1),topoIn->nglob(2));
    printf("[FLUPS] topo IN loc : %d*%d*%d = %d (check: %d %d %d)\n",topoIn->nmem(0),topoIn->nmem(1),topoIn->nmem(2),topoIn->memsize(),topoIn->nloc(0),topoIn->nloc(1),topoIn->nloc(2));
    printf("[FLUPS] topo OUT glob : %d %d %d \n",topoSpec->nglob(0),topoSpec->nglob(1),topoSpec->nglob(2));
    printf("[FLUPS] topo OUT loc  : nmem: %d*%d*%d nf:%d (nloc: %d %d %d)  \n",topoSpec->nmem(0),topoSpec->nmem(1),topoSpec->nmem(2),topoSpec->nf(),topoSpec->nloc(0),topoSpec->nloc(1),topoSpec->nloc(2));

    printf("[P3DFFT++] topo IN glob  : %d %d %d  \n",gdimsIN[0],gdimsIN[1],gdimsIN[2]);
    printf("[P3DFFT++] topo IN loc   : %d %d %d (is: %d %d %d) \n",P3DnlocIN[0],P3DnlocIN[1],P3DnlocIN[2],glob_startIN[0],glob_startIN[1],glob_startIN[2]);
    printf("[P3DFFT++] topo OUT glob : %d %d %d  \n",gdimsOUT[0],gdimsOUT[1],gdimsOUT[2]);
    printf("[P3DFFT++] topo OUT loc  : %d %d %d (is: %d %d %d) \n",P3DnlocOUT[0],P3DnlocOUT[1],P3DnlocOUT[2],glob_startOUT[0],glob_startOUT[1],glob_startOUT[2]);



    printf("I am going to allocate FLUPS: %d (inside FLUPS: %d) , P3D: %d (out %d C) \n",FLUmemsizeIN,FLUmemsizeOUT,P3DmemsizeIN,P3DmemsizeOUT);
    
 
    double *rhsFLU   = (double *)fftw_malloc(sizeof(double) * FLUmemsizeIN);
    //solFLU   allocated by sflups setup with the larger size
    double *rhsP3D   = (double *)fftw_malloc(sizeof(double) * P3DmemsizeIN);
    p3dfft::complex_double *solP3D   = (p3dfft::complex_double *)fftw_malloc(sizeof(p3dfft::complex_double) * P3DmemsizeOUT);

    std::memset(rhsFLU, 0, sizeof(double ) * FLUmemsizeIN);
    std::memset(solFLU, 0, sizeof(double ) * FLUmemsizeOUT); 
    std::memset(rhsP3D, 0, sizeof(double ) * P3DmemsizeIN);
    std::memset(solP3D, 0, sizeof(p3dfft::complex_double) * P3DmemsizeOUT);
    
    // printf("istart: %d %d %d =? %d %d %d(P3D)\n",istartGlo[0],istartGlo[1],istartGlo[2],glob_startIN[0],glob_startIN[1],glob_startIN[2]);

    double f = 1; //frequency of the wave

    for (int i2 = 0; i2 < topoIn->nloc(2); i2++) {
        for (int i1 = 0; i1 < topoIn->nloc(1); i1++) {
            for (int i0 = 0; i0 < topoIn->nloc(0); i0++) {
                double       x    = (istartGlo[0] + i0 ) * h[0]; //+ 0.5
                double       y    = (istartGlo[1] + i1 ) * h[1]; //+ 0.5
                double       z    = (istartGlo[2] + i2 ) * h[2]; //+ 0.5

                size_t id;
                id    = localIndex(0, i0, i1, i2, 0, FLUnmemIn, 1);
                rhsFLU[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
                
                //p3d does not care of the size you allocate, juste fill the first isize elements
                id =  localIndex(0, i0, i1, i2, 0, P3DnlocIN, 1);
                rhsP3D[id] = sin((c_2pi / L[0] * f) * x) * sin((c_2pi / L[1] * f) * y) * sin((c_2pi / L[2] * f) * z);
            }
        }
    }

    // //-------------------------------------------------------------------------
    // /** - Proceed to the solve */
    // //-------------------------------------------------------------------------
    
       
    double factor = 1.0/(nglob[0]* nglob[1]*nglob[2]);
    //Warmup
    // trans_f.exec(rhsP3D,solP3D,false);

    for (int iter=0;iter<n_iter;iter++){
        
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==0){
            printf("Iter %i\n",iter);
        }
        fflush(stdout);

        // ------------------P3D---------------:

        P3Dprof->start("FFTandSwitch");
        trans_f.exec(rhsP3D,solP3D,false);  // Execute forward real-to-complex FFT
        P3Dprof->stop("FFTandSwitch");

// #define PRINT_RES
#ifdef PRINT_RES
        /* normalize */
        for(int id = 0; id<P3DmemsizeOUT; id++){
            solP3D[id] *= factor;
        }
        print_res(solP3D,P3DnlocOUT,glob_startOUT);
#endif

        // ------------------FLUPS---------------:
        MPI_Barrier(MPI_COMM_WORLD);

        //reinit sol to prepare for inplace 3D FFT
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
#ifdef P3DMODIF
    double times[8][3];
    p3dfft::timers.get(times,MPI_COMM_WORLD);
    string P3DNames[8] = {"Reorder_trans","Reorder_out","Reorder_in","Trans_exec","Packsend","Packsend_trans","Unpackrecv","Alltoall"};
#endif

    double timeRef      = P3Dprof->get_timeAcc("FFTandSwitch");
    double timeInit     = P3Dprof->get_timeAcc("init");

    P3Dprof->disp("FFTandSwitch");
#ifdef P3DMODIF    
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
        for(int i=0; i<8; i++){
            int j = i;
            printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", ("-- "+P3DNames[j]).c_str(), 100.*times[j][0]/timeRef, 100*times[j][0]/timeRef, times[j][0], 0., times[j][0]/n_iter, times[j][1]/n_iter, times[j][2]/n_iter, (double) n_iter,0.);
            fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f;%09.2f\n", ("-- "+P3DNames[j]).c_str(), 100.*times[j][0]/timeRef, 100*times[j][0]/timeRef, times[j][0], 0., times[j][0]/n_iter, times[j][1]/n_iter, times[j][2]/n_iter, (double) n_iter,0.);
        }
        fclose(file);
    }   
#else
    p3dfft::timers.print(MPI_COMM_WORLD);
#endif
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
    p3dfft::cleanup();

    MPI_Finalize();
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
                const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem,2);
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
                const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem,1);
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                    
                    if (std::fabs(A[id+i0]) > 1.25e-4) {
                        printf("FLU(%d %d %d) %lg \n", i0 + gstart[ax0], i1 + gstart[ax1], i2 + gstart[ax2], A[id+i0]);
                    }
                }
            }
        }
    }
}

void print_res(p3dfft::complex_double *A,int *sdims,int *gstart)
{
  int x,y,z;
  p3dfft::complex_double *p;
  int imo[3],i,j;
  p = A;

  for(z=0;z < sdims[2];z++)
    for(y=0;y < sdims[1];y++)
      for(x=0;x < sdims[0];x++) {
        if(std::fabs(p->real())+std::fabs(p->imag()) > 1.25e-4) 
            printf("P3D(%d %d %d) %lg +i1* %lg\n",x+gstart[0],y+gstart[1],z+gstart[2],p->real(),p->imag());
        p++;
      }
}