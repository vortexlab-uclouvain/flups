/**
 * @file main_vtube.cpp
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
#include "vtube.hpp"

using namespace std;

//default values
const static int                 d_nsolve    = 1;
const static int                 d_nsample   = 1;
const static int                 d_startSize = 16;
const static int                 d_nprocs[3] = {1, 1, 1};
const static double              d_L[3]      = {1., 1., 1.};
const static FLUPS_GreenType     d_kernel    = CHAT_2;
const static FLUPS_BoundaryType  d_bcdef     = UNB;

static void print_help(){
    printf("This is FLUPS validation code: \n");
    printf(" --help  , -h :                 print this message\n");
    printf(" --nprocs, -np Ni Nj Nk :       Ni,Nj,Nk is the number of MPI processes in each direction\n");
    printf(" --resolution, -res Rx Ry Rz :  Rx,Ry,Rz is the total number of cells in each direction \n");
    printf(" --nresolution, -nres Nr :      Nr is the number of higher resolutions that will be tested, with a resolution (R * 2^[0:Nr-1])\n");
    printf(" --nsolve, -ns Ns :             Ns is the number of times each validation case will be run (for statistics on the profiler) \n");
    printf(" --type, -t t :                 t is the type of the solver: tube=0, ring=1 \n");
    printf(" --length, -L Lx Ly Lz :        Lx,Ly,Lz is the dimension of the physical domain \n");
    printf(" --kernel, -k {0-4}:            the Green kernel 0=CHAT2, 1=LGF2, 2=HEJ2, 3=HEJ4, 4=HEJ6 \n");
    printf(" --lda, -l {1,3}:               leading dimension of array, number of components (1=scalar, 3=vector)\n");
    printf(" --dir, -d {1,2,3}:             direction of the vortex tubes\n");
    printf(" --order, -o {1,2}:             derivation order asked for the velocity field\n");
    printf(" --sym_x, -smx {-1,0,1}:        say if another vortex tube is added in the X direction: -1 = odd symmetry, 0 = none, 1 = even symmetry\n");
    printf(" --sym_y, -smy {-1,0,1}:        say if another vortex tube is added in the Y direction: -1 = odd symmetry, 0 = none, 1 = even symmetry\n");
    
}

int static parse_args(int argc, char *argv[], int nprocs[3], double L[3],int sym[2], FLUPS_GreenType *kernel, int *nsample, int **size, int *nsolve, int* type, int * vdir, FLUPS_DiffType* order){

    int startSize[3] = {d_startSize,d_startSize,d_startSize};

    //assigning default values
    for (int i = 0; i < 3; i++) {
        nprocs[i]   = d_nprocs[i];
        L[i]        = d_L[i];
    }
    *nsolve  = d_nsolve;
    *nsample = d_nsample;
    *kernel  = d_kernel;
    sym[0]   = 0; // no other tube by def
    sym[1]   = 0; // no other tube by def
    *type = 0; // ring by def
    *order = SPE;
    *vdir = 2;

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
            i += 3;
        } else if ((arg == "-nres") || (arg == "--nresolution")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                *nsample = atoi(argv[i + 1]);
                if (*nsample < 1) {
                    fprintf(stderr, "nresolution must be >0\n");
                    return 1;
                }
            } else {  //Missing argument
                fprintf(stderr, "missing -nresolution\n");
                return 1;
            }
            i++;
        } else if ((arg == "-t") || (arg == "--type")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                *type = atoi(argv[i + 1]);
            } else {  //Missing argument
                fprintf(stderr, "missing -nresolution\n");
                return 1;
            }
            i++;
        } else if ((arg == "-res") || (arg == "--resolution")) {
            for (int j = 0; j < 3; j++) {
                if (i + j + 1 < argc) {  // Make sure we aren't at the end of argv!
                    startSize[j] = atoi(argv[i + j + 1]);
                    if (startSize[j] <= 0.0) {
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
        } else if ((arg == "-d") || (arg == "--dir")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                *vdir = atoi(argv[i+1]); 
            } else { //Missing argument
                fprintf(stderr, "missing --dir\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-o") || (arg == "--order")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                *order = (FLUPS_DiffType) atoi(argv[i+1]); 
            } else { //Missing argument
                fprintf(stderr, "missing --order\n");
                return 1;
            }  
            i++;
        } else if ((arg == "-smx") || (arg == "--sym_x")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                sym[0] = atoi(argv[i + 1]);
            } else {  //Missing argument
                fprintf(stderr, "missing --sym_x\n");
                return 1;
            }
            i++;
        } else if ((arg == "-smy") || (arg == "--sym_y")) {
            if (i + 1 < argc) {  // Make sure we aren't at the end of argv!
                sym[1] = atoi(argv[i + 1]);
            } else {  //Missing argument
                fprintf(stderr, "missing --sym_y\n");
                return 1;
            }
            i++;
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
    int nsolve, nsample;
    int nprocs[3];
    double L[3];
    int *size = NULL;
    FLUPS_GreenType kernel;
    FLUPS_BoundaryType bcdef[3][2][3];
    int sym[2];
    int type;
    FLUPS_DiffType order;
    int vdir;

    int status = parse_args(argc, argv, nprocs, L, sym, &kernel, &nsample, &size, &nsolve, &type, &vdir, &order);

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

    if (provided < requested) {
        printf("The MPI-provided thread behavior does not match\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Display
    if (rank == 0) {
        printf("I will run with:\n");
        printf("  --nprocs: %d,%d,%d\n", nprocs[0], nprocs[1], nprocs[2]);
        printf("  -L: %lf,%lf,%lf\n", L[0], L[1], L[2]);
        printf("  --kernel: %d\n", kernel);
        printf("  --nsolve: %d\n", nsolve);
        printf("  --vdir: %d\n", vdir);
        printf("  --order: %d\n", order);
        for (int i = 0; i < nsample; i++) {
            printf("   -> sample %d: %d %d %d\n", i + 1, size[i * 3], size[i * 3 + 1], size[i * 3 + 2]);
        }
    }

    // Do the validation
    struct DomainDescr valCase;
    for (int is = 0; is < nsample; is++) {
        for (int ip = 0; ip < 3; ip++) {
            valCase.L[ip]     = L[ip];
            valCase.nproc[ip] = nprocs[ip];
            valCase.nglob[ip] = size[is * 3 + ip];
            if (type == 0) {// this is the vortex tube
                // the BC are imposed in the Z direction
                int dir2 = vdir;
                int dir0 = (vdir+1)%3;
                int dir1 = (vdir+2)%3;

                if (ip == dir2) {
                    valCase.mybcv[ip][0][dir0] = ODD;   // w_x on the left
                    valCase.mybcv[ip][0][dir1] = ODD;   // w_y on the left
                    valCase.mybcv[ip][0][dir2] = EVEN;  // w_z on the left
                    valCase.mybcv[ip][1][dir0] = ODD;   // w_x on the right
                    valCase.mybcv[ip][1][dir1] = ODD;   // w_y on the right
                    valCase.mybcv[ip][1][dir2] = EVEN;  // w_z on the right
                } else if (ip == dir1) {
                    // the Y bc depends on the symmetry
                    if (sym[1] == 1) {
                        valCase.ycntr           = 0.25;
                        valCase.ysign           = 1.0;
                        valCase.mybcv[ip][0][dir0] = EVEN;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = ODD;   // w_y on the left
                        valCase.mybcv[ip][0][dir2] = EVEN;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;   // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;   // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;   // w_z on the right
                    } else if (sym[1] == 0) {
                        valCase.ycntr           = 0.5;
                        valCase.ysign           = 0.0;
                        valCase.mybcv[ip][0][dir0] = UNB;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = UNB;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = UNB;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;  // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;  // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;  // w_z on the right
                    } else if (sym[1] == -1) {
                        valCase.ycntr           = 0.25;
                        valCase.ysign           = -1.0;
                        valCase.mybcv[ip][0][dir0] = ODD;   // w_x on the left
                        valCase.mybcv[ip][0][dir1] = EVEN;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = ODD;   // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;   // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;   // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;   // w_z on the right
                    }
                } else {
                    if (sym[0] == 1) {
                        valCase.xcntr           = 0.25;
                        valCase.xsign           = 1.0;
                        valCase.mybcv[ip][0][dir0] = ODD;   // w_x on the left
                        valCase.mybcv[ip][0][dir1] = EVEN;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = EVEN;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;   // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;   // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;   // w_z on the right
                    } else if (sym[0] == 0) {
                        valCase.xcntr           = 0.5;
                        valCase.xsign           = 0.0;
                        valCase.mybcv[ip][0][dir0] = UNB;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = UNB;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = UNB;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;  // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;  // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;  // w_z on the right
                    } else if (sym[0] == -1) {
                        valCase.xcntr           = 0.25;
                        valCase.xsign           = -1.0;
                        valCase.mybcv[ip][0][dir0] = EVEN;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = ODD;   // w_y on the left
                        valCase.mybcv[ip][0][dir2] = ODD;   // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;   // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;   // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;   // w_z on the right
                    }
                }
            } else if (type == 1) {  // this is the vortex ring
                                     // the BC are imposed in the Z direction
                valCase.xcntr = 0.5;
                valCase.ycntr = 0.5;
                valCase.zcntr = 0.5;
                valCase.xsign = 1.0;
                valCase.ysign = 1.0;
                valCase.zsign = 1.0;

                valCase.mybcv[ip][0][0] = UNB;  // w_x on the left
                valCase.mybcv[ip][0][1] = UNB;  // w_y on the left
                valCase.mybcv[ip][0][2] = UNB;  // w_z on the left
                valCase.mybcv[ip][1][0] = UNB;  // w_x on the right
                valCase.mybcv[ip][1][1] = UNB;  // w_y on the right
                valCase.mybcv[ip][1][2] = UNB;  // w_z on the right
            }
        }
        // let's gooo
        vtube(valCase, kernel, nsolve, type, order, vdir);
    }

    free(size);

    MPI_Finalize();
}
