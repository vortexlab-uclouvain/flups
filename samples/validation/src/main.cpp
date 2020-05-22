/**
 * @file main.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include <cmath>
#include <iostream>

#include "mpi.h"
#include "flups.h"
#include <iostream>
#include <cstring>
#include "validation_3d.hpp"

using namespace std;

//default values
const static int                 d_nsolve    = 1;
const static int                 d_nsample   = 1;
const static int                 d_startSize = 64;
const static int                 d_nprocs[3] = {1, 1, 1};
const static double              d_L[3]      = {1., 1., 1.};
const static FLUPS_GreenType     d_kernel    = CHAT_2;
const static FLUPS_BoundaryType  d_bcdef     = ODD;
const static int                 d_lda       = 1;

static void print_help(){
    printf("This is FLUPS validation code: \n");
    printf(" --help  , -h :                 print this message\n");
    printf(" --nprocs, -np Ni Nj Nk :       Ni,Nj,Nk is the number of MPI processes in each direction\n");
    printf(" --resolution, -res Rx Ry Rz :  Rx,Ry,Rz is the total number of cells in each direction \n");
    printf(" --nresolution, -nres Nr :      Nr is the number of higher resolutions that will be tested, with a resolution (R * 2^[0:Nr-1])\n");
    printf(" --nsolve, -ns Ns :             Ns is the number of times each validation case will be run (for statistics on the profiler) \n");
    printf(" --length, -L Lx Ly Lz :        Lx,Ly,Lz is the dimension of the physical domain \n");
    printf(" --kernel, -k {0-4}:            the Green kernel 0=CHAT2, 1=LGF2, 2=HEJ2, 3=HEJ4, 4=HEJ6 \n");
    printf(" --lda, -l {1,3}:               leading dimension of array, number of components (1=scalar, 3=vector)\n");
    printf(" --boundary-conditions, -bc     \n ");
    printf("     Bxl Bxr Byl Byr Bzl Bzr : the boundary conditions in x/y/z on each side l/r. 0=EVEN, 1=ODD, 3=PERiodic, 4=UNBounded \n");
    printf(" --boundary-conditionv, -bcv     \n ");
    printf("     3 x (Bxl Bxr Byl Byr Bzl Bzr) : the boundary conditions in x/y/z on each side l/r, 3 times for each component. 0=EVEN, 1=ODD, 3=PERiodic, 4=UNBounded \n");
}

int static parse_args(int argc, char *argv[], int nprocs[3], double L[3], FLUPS_BoundaryType bcdef[3][2],  FLUPS_BoundaryType bcdefv[3][2][3], FLUPS_GreenType *kernel, int *lda, int *nsample, int **size, int *nsolve){

    int startSize[3] = {d_startSize,d_startSize,d_startSize};

    //assigning default values
    for (int i = 0; i < 3; i++) {
        nprocs[i]   = d_nprocs[i];
        L[i]        = d_L[i];
        bcdef[i][0] = d_bcdef;
        bcdef[i][1] = d_bcdef;

        bcdefv[i][0][0] = NONE;
        bcdefv[i][0][1] = NONE;
        bcdefv[i][0][2] = NONE;
        bcdefv[i][1][0] = NONE;
        bcdefv[i][1][1] = NONE;
        bcdefv[i][1][2] = NONE;
    }
    *nsolve  = d_nsolve;
    *nsample = d_nsample;
    *kernel  = d_kernel;  
    *lda     = d_lda;

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
                    if(L[j]<0.0){
                        fprintf(stderr, "L must be >=0\n");
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
            for (int j = 0; j<3;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    startSize[j] = atoi(argv[i+j+1]); 
                    if(startSize[j]<=0.0){
                        fprintf(stderr, "res must be >0\n");
                        return 1;
                    }
                } else { //Missing argument
                    fprintf(stderr, "missing argument in -res\n");
                    return 1;
                }  
            }
            i+=3;
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
                *kernel = (FLUPS_GreenType) atoi(argv[i+1]); 
            } else { //Missing argument
                fprintf(stderr, "missing --kernel\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-l") || (arg == "--lda")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                *lda = atoi(argv[i+1]); 
            } else { //Missing argument
                fprintf(stderr, "missing --lda\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-bc")|| (arg== "--boundary-conditions") ) {
            for (int j = 0; j<6;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    bcdef[j/2][j%2] = (FLUPS_BoundaryType) atoi(argv[i+j+1]); 
                } else { //Missing argument
                    fprintf(stderr, "missing argument in --boundary-conditions\n");
                    return 1;
                }  
            }
            i+=6;
        }
         else if ((arg == "-bcv")|| (arg== "--boundary-conditionsv") ) {
            for (int j = 0; j<18;j++){
                if (i + j + 1 < argc) { // Make sure we aren't at the end of argv!
                    bcdefv[(j/2)%3][j%2][j/6] = (FLUPS_BoundaryType) atoi(argv[i+j+1]); 
                } else { //Missing argument
                    fprintf(stderr, "missing argument in --boundary-conditions\n");
                    return 1;
                }  
            }
            i+=18;
        } 
    }
    
    // finilizing allocations
    *size =(int*) malloc((*nsample) * 3 * sizeof(int));
    
    for (int i = 0; i<*nsample*3 ; i+=3){
        (*size)[i]   = startSize[0] * pow(2,i/3);
        (*size)[i+1] = startSize[1] * pow(2,i/3);
        (*size)[i+2] = startSize[2] * pow(2,i/3);
    }
    return 0;
}

int main(int argc, char *argv[]) {

    // Parsing arguments
    int nsolve, nsample, lda;
    int nprocs[3];
    double L[3];
    int *size = NULL;
    FLUPS_GreenType kernel;
    FLUPS_BoundaryType bcdef[3][2];
    FLUPS_BoundaryType bcdefv[3][2][3];
    
    int status = parse_args(argc, argv, nprocs, L, bcdef, bcdefv, &kernel, &lda, &nsample, &size, &nsolve );

    if (status) exit(status);
    if (size==NULL){
        exit(0); //we just printed help
    }     

    // Initialize MPI
    int rank;

    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        printf("The MPI-provided thread behavior does not match\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Do the validation
    struct DomainDescr valCase;
    valCase.dovectorbc = (bcdefv[0][0][0] != NONE);

    // Display
    if (rank == 0) {
        printf("I will run with:\n");
        printf("  --nprocs: %d,%d,%d\n", nprocs[0], nprocs[1], nprocs[2]);
        printf("  -L: %lf,%lf,%lf\n", L[0], L[1], L[2]);
        if (!valCase.dovectorbc) {
            printf("  -bc: %d,%d ; %d,%d ; %d,%d\n", bcdef[0][0], bcdef[0][1], bcdef[1][0], bcdef[1][1], bcdef[2][0], bcdef[2][1]);
        } else {
            printf("  -bcv: [%d,%d ; %d,%d ; %d,%d],[%d,%d ; %d,%d ; %d,%d],[%d,%d ; %d,%d ; %d,%d]\n", bcdefv[0][0][0], bcdefv[0][1][0], bcdefv[1][0][0], bcdefv[1][1][0], bcdefv[2][0][0], bcdefv[2][1][0],
                   bcdefv[0][0][1], bcdefv[0][1][1], bcdefv[1][0][1], bcdefv[1][1][1], bcdefv[2][0][1], bcdefv[2][1][1],
                   bcdefv[0][0][2], bcdefv[0][1][2], bcdefv[1][0][2], bcdefv[1][1][2], bcdefv[2][0][2], bcdefv[2][1][2]);
        }
        printf("  --kernel: %d\n", kernel);
        printf("  --lda: %d\n", lda);
        printf("  --nsolve: %d\n", nsolve);
        for (int i = 0; i < nsample; i++) {
            printf("   -> sample %d: %d %d %d\n", i + 1, size[i * 3], size[i * 3 + 1], size[i * 3 + 2]);
        }
    }

    for (int is = 0; is < nsample; is++) {
        for (int ip = 0; ip < 3; ip++) {
            valCase.L[ip]       = L[ip];
            valCase.nproc[ip]   = nprocs[ip];
            valCase.nglob[ip]   = size[is * 3 + ip];
            valCase.mybc[ip][0] = bcdef[ip][0];
            valCase.mybc[ip][1] = bcdef[ip][1];

            valCase.mybcv[ip][0][0] = bcdefv[ip][0][0];
            valCase.mybcv[ip][0][1] = bcdefv[ip][0][1];
            valCase.mybcv[ip][0][2] = bcdefv[ip][0][2];
            valCase.mybcv[ip][1][0] = bcdefv[ip][1][0];
            valCase.mybcv[ip][1][1] = bcdefv[ip][1][1];
            valCase.mybcv[ip][1][2] = bcdefv[ip][1][2];
        }
        // let's gooo
        validation_3d(valCase, kernel, lda, nsolve);
    }

    free(size);

    MPI_Finalize();
}
