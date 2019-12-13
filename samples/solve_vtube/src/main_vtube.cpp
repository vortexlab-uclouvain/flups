/**
 * @file main_vtube.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
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
    printf(" --type, -t t :                 t is the type of the solver: STD=0, RHS=1 \n");
    printf(" --length, -L Lx Ly Lz :        Lx,Ly,Lz is the dimension of the physical domain \n");
    printf(" --kernel, -k {0-4}:            the Green kernel 0=CHAT2, 1=LGF2, 2=HEJ2, 3=HEJ4, 4=HEJ6 \n");
    printf(" --lda, -l {1,3}:               leading dimension of array, number of components (1=scalar, 3=vector)\n");
    printf(" --sym_x, -smx {-1,0,1}:        say if another vortex tube is added in the X direction: -1 = odd symmetry, 0 = none, 1 = even symmetry\n");
    printf(" --sym_y, -smy {-1,0,1}:        say if another vortex tube is added in the Y direction: -1 = odd symmetry, 0 = none, 1 = even symmetry\n");
    
}

int static parse_args(int argc, char *argv[], int nprocs[3], double L[3],int sym[2], FLUPS_GreenType *kernel, int *nsample, int **size, int *nsolve){

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

    int status = parse_args(argc, argv, nprocs, L, sym, &kernel, &nsample, &size, &nsolve);

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
            // the BC are imposed in the Z direction
            if (ip == 2) {
                valCase.mybcv[ip][0][0] = ODD;   // w_x on the left
                valCase.mybcv[ip][0][1] = ODD;   // w_y on the left
                valCase.mybcv[ip][0][2] = EVEN;  // w_z on the left
                valCase.mybcv[ip][1][0] = ODD;   // w_x on the right
                valCase.mybcv[ip][1][1] = ODD;   // w_y on the right
                valCase.mybcv[ip][1][2] = EVEN;  // w_z on the right
            } else if (ip == 1) {
                // the Y bc depends on the symmetry
                if (sym[ip] == 1) {
                    valCase.ycntr           = 0.25;
                    valCase.ysign           = 1.0;
                    valCase.mybcv[ip][0][0] = EVEN;  // w_x on the left
                    valCase.mybcv[ip][0][1] = ODD;   // w_y on the left
                    valCase.mybcv[ip][0][2] = EVEN;  // w_z on the left
                    valCase.mybcv[ip][1][0] = UNB;   // w_x on the right
                    valCase.mybcv[ip][1][1] = UNB;   // w_y on the right
                    valCase.mybcv[ip][1][2] = UNB;   // w_z on the right
                } else if (sym[ip] == 0) {
                    valCase.ycntr           = 0.5;
                    valCase.ysign           = 0.0;
                    valCase.mybcv[ip][0][0] = UNB;  // w_x on the left
                    valCase.mybcv[ip][0][1] = UNB;  // w_y on the left
                    valCase.mybcv[ip][0][2] = UNB;  // w_z on the left
                    valCase.mybcv[ip][1][0] = UNB;  // w_x on the right
                    valCase.mybcv[ip][1][1] = UNB;  // w_y on the right
                    valCase.mybcv[ip][1][2] = UNB;  // w_z on the right
                } else if (sym[ip] == -1) {
                    valCase.ycntr           = 0.25;
                    valCase.ysign           = -1.0;
                    valCase.mybcv[ip][0][0] = ODD;   // w_x on the left
                    valCase.mybcv[ip][0][1] = EVEN;  // w_y on the left
                    valCase.mybcv[ip][0][2] = ODD;   // w_z on the left
                    valCase.mybcv[ip][1][0] = UNB;   // w_x on the right
                    valCase.mybcv[ip][1][1] = UNB;   // w_y on the right
                    valCase.mybcv[ip][1][2] = UNB;   // w_z on the right
                }
            } else {
                if (sym[ip] == 1) {
                    valCase.xcntr           = 0.25;
                    valCase.xsign           = 1.0;
                    valCase.mybcv[ip][0][0] = ODD;   // w_x on the left
                    valCase.mybcv[ip][0][1] = EVEN;  // w_y on the left
                    valCase.mybcv[ip][0][2] = EVEN;  // w_z on the left
                    valCase.mybcv[ip][1][0] = UNB;   // w_x on the right
                    valCase.mybcv[ip][1][1] = UNB;   // w_y on the right
                    valCase.mybcv[ip][1][2] = UNB;   // w_z on the right
                } else if (sym[ip] == 0) {
                    valCase.xcntr           = 0.5;
                    valCase.xsign           = 0.0;
                    valCase.mybcv[ip][0][0] = UNB;  // w_x on the left
                    valCase.mybcv[ip][0][1] = UNB;  // w_y on the left
                    valCase.mybcv[ip][0][2] = UNB;  // w_z on the left
                    valCase.mybcv[ip][1][0] = UNB;  // w_x on the right
                    valCase.mybcv[ip][1][1] = UNB;  // w_y on the right
                    valCase.mybcv[ip][1][2] = UNB;  // w_z on the right
                } else if (sym[ip] == -1) {
                    valCase.xcntr           = 0.25;
                    valCase.xsign           = -1.0;
                    valCase.mybcv[ip][0][0] = EVEN;  // w_x on the left
                    valCase.mybcv[ip][0][1] = ODD;   // w_y on the left
                    valCase.mybcv[ip][0][2] = ODD;   // w_z on the left
                    valCase.mybcv[ip][1][0] = UNB;   // w_x on the right
                    valCase.mybcv[ip][1][1] = UNB;   // w_y on the right
                    valCase.mybcv[ip][1][2] = UNB;   // w_z on the right
                }
            }
        }
        // let's gooo
        vtube(valCase, kernel, nsolve);
    }

    free(size);

    MPI_Finalize();
}
