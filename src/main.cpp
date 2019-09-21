/**
 * @file main.cpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include <cmath>
#include <iostream>

#include "expint.hpp"

#include "Solver.hpp"
#include "SwitchTopo.hpp"
#include "validation_3d.hpp"

#include "mpi.h"

//default values
const static int                 d_nsolve    = 1;
const static int                 d_nsample   = 1;
const static int                 d_startSize = 16;
const static int                 d_nprocs[3] = {1, 1, 1};
const static double              d_L[3]      = {1., 1., 1.};
const static FLUPS::GreenType    d_kernel    = FLUPS::CHAT_2;
const static FLUPS::BoundaryType d_bcdef     = FLUPS::UNB;

static void print_help(){
    printf("This is FLUPS validation code: \n");
    printf(" --help  , -h :                 print this message\n");
    printf(" --nprocs, -np Ni Nj Nk :       Ni,Nj,Nk is the number of MPI processes in each direction\n");
    printf(" --resolution, -res R :         R is the number of cells per unit length \n");
    printf(" --nresolution, -nres Nr :      Nr is the number of higher resolutions that will be tested, with a resolution (R * 2^[0:Nr-1])\n");
    printf(" --nsolve, -ns Ns :             Ns is the number of times each validation case will be run (for statistics on the profiler) \n");
    printf(" --length, -L Lx Ly Lz :        Lx,Ly,Lz is the dimension of the physical domain \n");
    printf(" --kernel, -k [0-3]:            the Green kernel 0=CHAT2, 1=HEJ2, 2=HEJ4, 3=HEJ6 \n");
    printf(" --boundary-conditions, -bc     \n ");
    printf("     Bxl Bxr Byl Byr Bzl Bzr : the boundary conditions in x/y/z on each side l/r. 0=EVEN, 1=ODD, 3=PERiodic, 4=UNBounded \n");
    printf(" --predefined-test, -pt :       runs a predefined validation test with several combination of UNB BCs and all the Green Kernels (excludes -L, -k and -bc) \n ");
}

int static parse_args(int argc, char *argv[], int nprocs[3], double L[3], FLUPS::BoundaryType bcdef[3][2], int *predef, FLUPS::GreenType *kernel, int *nsample, int **size, int *nsolve){
    BEGIN_FUNC;

    int startSize = d_startSize;

    //assigning default values
    for (int i = 0; i < 3; i++) {
        nprocs[i]   = d_nprocs[i];
        L[i]        = d_L[i];
        bcdef[i][0] = d_bcdef;
        bcdef[i][1] = d_bcdef;
    }
    *nsolve  = d_nsolve;
    *nsample = d_nsample;
    *kernel  = d_kernel;  
    *predef  = 0;

    // modifying if necessary
    if(argc < 1 || ( argc==1 && (!argv[0][0] || strcmp(argv[0],"flups_validation"))) ){
        print_help();
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            print_help();
            return 0;
        } else if ((arg == "-np") || (arg == "--nprocs")) {
            for (int j = 0; j<3;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    nprocs[j] = atoi(argv[i+j+1]); 
                    if(nprocs[j]<1){
                        fprintf(stderr, "nprocs must be >0\n");
                        return 1;
                    }
                } else { //Missing argument
                    fprintf(stderr, "missing argument in --nprocs\n");
                    return 1;
                }  
            }
            i+=3;
        } else if ((arg == "-L") || (arg == "-length") ) {
            for (int j = 0; j<3;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    L[j] = atof(argv[i+j+1]); 
                    if(L[j]<=0.0){
                        fprintf(stderr, "L must be >0\n");
                        return 1;
                    }
                } else { //Missing argument
                    fprintf(stderr, "missing argument in -L\n");
                    return 1;
                }  
            }
            i+=3;
        } else if ((arg == "-nres")|| (arg== "--nresolution") ) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                *nsample = atoi(argv[i+1]); 
                if(*nsample<1){
                    fprintf(stderr, "nresolution must be >0\n");
                    return 1;
                }
            } else { //Missing argument
                fprintf(stderr, "missing -nresolution\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-res")|| (arg== "--resolution") ) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                startSize = atoi(argv[i+1]); 
                if(startSize<1){
                    fprintf(stderr, "resolution must be >0\n");
                    return 1;
                }
            } else { //Missing argument
                fprintf(stderr, "missing --resolution\n");
                return 1;
            }  
            i++;
        }  else if ((arg == "-ns")|| (arg== "--nsolve") ) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                *nsolve = atoi(argv[i+1]); 
                if(*nsolve<1){
                    fprintf(stderr, "nsolve must be >0\n");
                    return 1;
                }
            } else { //Missing argument
                fprintf(stderr, "missing --nsolve\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-k") || (arg == "--kernel")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                *kernel = (FLUPS::GreenType) atoi(argv[i+1]); 
            } else { //Missing argument
                fprintf(stderr, "missing --kernel\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-bc")|| (arg== "--boundary-conditions") ) {
            for (int j = 0; j<6;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    bcdef[j/2][j%2] = (FLUPS::BoundaryType) atoi(argv[i+j+1]); 
                } else { //Missing argument
                    fprintf(stderr, "missing argument in --boundary-conditions\n");
                    return 1;
                }  
            }
            i+=6;
        } else if ((arg == "-pt")|| (arg== "--predefined-test") ) {
            *predef = 1;
        }
    }
    
    // finilizing allocations
    *size =(int*) malloc((*nsample) * sizeof(int));
    
    for (int i = 0; i<*nsample ; i++){
        (*size)[i] = startSize * pow(2,i);
    }
    return 0;
}

int main(int argc, char *argv[]) {

    // Parsing arguments
    int nsolve, nsample, predef;
    int nprocs[3];
    double L[3];
    int *size;
    FLUPS::GreenType kernel;
    FLUPS::BoundaryType bcdef[3][2];
    
    int status = parse_args(argc, argv, nprocs, L, bcdef, &predef, &kernel, &nsample, &size, &nsolve);

    if (status) exit(status);
    if (!size){
        exit(0); //we just printed help
    }     

    // Initialize MPI
    int rank;

    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        FLUPS_ERROR("The MPI-provided thread behavior does not match", LOCATION);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Do the validation
    if (predef == 0){
        // Display
        if (rank == 0) {
            FLUPS_INFO("I will run with:");
            FLUPS_INFO("  --nprocs: %d,%d,%d", nprocs[0], nprocs[1], nprocs[2]);
            FLUPS_INFO("  -L: %lf,%lf,%lf", L[0], L[1], L[2]);
            FLUPS_INFO("  -bc: %d,%d ; %d,%d ; %d,%d", bcdef[0][0], bcdef[0][1], bcdef[1][0], bcdef[1][1], bcdef[2][0], bcdef[2][1]);
            FLUPS_INFO("  --kernel: %d", kernel);
            FLUPS_INFO("  --nsolve: %d", nsolve);
            for (int i = 0; i < nsample; i++) {
                FLUPS_INFO("   -> sample %d: %d", i + 1, size[i]);
            }
        }

        for (int is = 0; is < nsample; is++) {
            struct DomainDescr valCase;

            for (int ip = 0; ip < 3; ip++) {
                valCase.L[ip]       = L[ip];
                valCase.nproc[ip]   = nprocs[ip];
                valCase.nglob[ip]   = size[is] * L[ip];
                valCase.mybc[ip][0] = bcdef[ip][0];
                valCase.mybc[ip][1] = bcdef[ip][1];
            }
            validation_3d(valCase, FLUPS::SRHS, kernel, nsolve);
        }

    } else {
        // Loop over the resolution for convergence study
        for (int is = 0; is < nsample; is++) {
            // Display
            if (rank == 0) {
                printf("I will run predefined test 1 with:");
                printf("  --nprocs: %d,%d,%d", nprocs[0], nprocs[1], nprocs[2]);
                printf("  --nsolve: %d", nsolve);
                for (int i = 0; i < nsample; i++) {
                    printf("   -> sample %d: %d", i + 1, size[i]);
                }
            }

            struct DomainDescr valCase;

            for (int ip = 0; ip < 3; ip++) {
                valCase.nproc[ip]   = nprocs[ip];
                valCase.nglob[ip]   = size[is] * valCase.L[ip];
            }

            if (rank == 0) printf("\n==============================  FULLY UNBOUNDED TEST ============================================\n");
            valCase.mybc[0][0] = FLUPS::UNB;
            valCase.mybc[0][1] = FLUPS::UNB;
            valCase.mybc[1][0] = FLUPS::UNB;
            valCase.mybc[1][1] = FLUPS::UNB;
            valCase.mybc[2][0] = FLUPS::UNB;
            valCase.mybc[2][1] = FLUPS::UNB;

            if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
            validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_2, nsolve);
            if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
            validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_4, nsolve);
            if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
            validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_6, nsolve);
            if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
            validation_3d(valCase, FLUPS::SRHS, FLUPS::CHAT_2, nsolve);

            if (rank == 0) printf("\n==============================  MIX UNBOUNDED TEST ============================================\n");
            for (int id = 0; id < 3; id++) {
                for (int lr = 0; lr < 2; lr++) {
                    valCase.mybc[0][0] = FLUPS::UNB;
                    valCase.mybc[0][1] = FLUPS::UNB;
                    valCase.mybc[1][0] = FLUPS::UNB;
                    valCase.mybc[1][1] = FLUPS::UNB;
                    valCase.mybc[2][0] = FLUPS::UNB;
                    valCase.mybc[2][1] = FLUPS::UNB;
                    // set to even
                    valCase.mybc[id][lr] = FLUPS::EVEN;
                    if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_2);
                    if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_4);
                    if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_6);
                    if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::CHAT_2);
                    // set to odd
                    valCase.mybc[id][lr] = FLUPS::ODD;
                    if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_2);
                    if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_4);
                    if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_6);
                    if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
                    validation_3d(valCase, FLUPS::SRHS, FLUPS::CHAT_2);
                }
            }
        }
    }
    free(size);

    MPI_Finalize();
}
