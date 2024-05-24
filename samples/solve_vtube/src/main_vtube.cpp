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
#include <iostream>
#include <cstring>

#include "h3lpr/profiler.hpp"
#include "h3lpr/parser.hpp"
#include "flups.h"


#include "vtube.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    //--------------------------------------------------------------------------
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

    //--------------------------------------------------------------------------
    // Parsing
    //--------------------------------------------------------------------------
    H3LPR::Parser parser(argc, (const char **)argv);

    // Retreive the int values
    bool arg_compact   = parser.GetFlag("--compact", "uses the compact gaussian instead of the infinite one");
    auto arg_center    = parser.GetValue<int>("--center", "Indicate the location of the data: 0=Node-centred, 1=Cell-centred", 1);
    auto arg_nsolve    = parser.GetValue<int>("--nsolve", "the number of solves to perform", 1);
    auto arg_kernel    = parser.GetValue<int>("--kernel", "the Green kernel 0=CHAT2, 1=LGF2, 2=HEJ2, 3=HEJ4, 4=HEJ6", 0);
    auto arg_sol_t     = parser.GetValue<int>("--case", "the testcase considered: tube=0, ring=1, compact blob=3, periodic=4", 0);
    auto arg_nsample   = parser.GetValue<int>("--nres", "Nr is the number of higher resolutions that will be tested, with a resolution (R * 2^[0:Nr-1])", 1);
    auto arg_direction = parser.GetValue<int>("--dir", "direction of the vortex tubes (1=x,2=y,3=z)", 2);
    auto arg_order     = parser.GetValue<int>("--order", "derivation order asked for the velocity field (1=spectral, 2=FD2, 4=FD4, 6=FD6) ", 1);
    auto arg_sym_x     = parser.GetValue<int>("--sym_x", "say if another vortex tube is added in the X direction: -1 = odd symmetry, 0 = none, 1 = even symmetry ", 0);
    auto arg_sym_y     = parser.GetValue<int>("--sym_y", "ay if another vortex tube is added in the Y direction: -1 = odd symmetry, 0 = none, 1 = even symmetry", 0);
    auto arg_rc        = parser.GetValue<double>("--rc", "the cutoff distance for compact fields", 0.2);
    auto arg_rad       = parser.GetValue<double>("--rad", "radius of the ring case", 0.25);
    auto arg_sigma     = parser.GetValue<double>("--sigma", "the sigma of the exponential", 0.05);

    // Retreive the vector value
    auto arg_nprocs = parser.GetValues<int, 3>("--np", "the number of processes in each direction", {1, 1, 1});
    auto arg_nres   = parser.GetValues<int, 3>("--res", "the number of unknowns each direction", {16, 16, 16});
    auto arg_L      = parser.GetValues<double, 3>("--dom", "the size of the domain each direction", {1., 1., 1.});
    parser.Finalize();

    int *size = (int *)malloc((arg_nsample)*3 * sizeof(int));
    if (size == NULL) {
        exit(0);  // we just printed help
    }

    for (int i = 0; i < arg_nsample * 3; i += 3){
        size[i]   = arg_nres[0] * pow(2,i/3);
        size[i+1] = arg_nres[1] * pow(2,i/3);
        size[i+2] = arg_nres[2] * pow(2,i/3);
    }

    //..........................................................................
    // Display
    if (rank == 0) {
        printf("I will run with:\n");
        printf("  --nprocs: %d,%d,%d\n", arg_nprocs[0], arg_nprocs[1], arg_nprocs[2]);
        printf("  -L: %lf,%lf,%lf\n", arg_L[0], arg_L[1], arg_L[2]);
        printf("  --kernel: %d\n", arg_kernel);
        printf("  --center: %d\n", arg_center);
        printf("  --nsolve: %d\n", arg_nsolve);
        printf("  --case: %d\n", arg_sol_t);
        printf("  --dir: %d\n", arg_direction);
        printf("  --order: %d\n", arg_order);
        printf("  --sym_x: %d\n", arg_sym_x);
        printf("  --sym_y: %d\n", arg_sym_y);
        printf("  --rc: %e\n", arg_rc);
        printf("  --rad: %e\n", arg_rad);
        printf("  --sigma: %e\n", arg_sigma);
        printf("  --compact? %d\n", arg_compact);
        for (int i = 0; i < arg_nsample; i++) {
            printf("   -> sample %d: %d %d %d\n", i + 1, size[i * 3], size[i * 3 + 1], size[i * 3 + 2]);
        }
    }

    //--------------------------------------------------------------------------
    // Do the validation
    //--------------------------------------------------------------------------
    struct DomainDescr valCase;
    for (int is = 0; is < arg_nsample; is++) {
        valCase.rc      = arg_rc;
        valCase.sigma   = arg_sigma;
        valCase.rad     = arg_rad;
        valCase.compact = arg_compact;
        for (int ip = 0; ip < 3; ip++) {
            valCase.L[ip]      = arg_L[ip];
            valCase.nproc[ip]  = arg_nprocs[ip];
            valCase.nglob[ip]  = arg_nres[is * 3 + ip];
            valCase.center[ip] = (FLUPS_CenterType)arg_center;

            if (arg_sol_t == 0) {  // this is the vortex tube
                // the BC are imposed in the Z direction
                int dir2 = arg_direction;
                int dir0 = (arg_direction + 1) % 3;
                int dir1 = (arg_direction + 2) % 3;

                if (ip == dir2) {
                    valCase.mybcv[ip][0][dir0] = ODD;   // w_x on the left
                    valCase.mybcv[ip][0][dir1] = ODD;   // w_y on the left
                    valCase.mybcv[ip][0][dir2] = EVEN;  // w_z on the left
                    valCase.mybcv[ip][1][dir0] = ODD;   // w_x on the right
                    valCase.mybcv[ip][1][dir1] = ODD;   // w_y on the right
                    valCase.mybcv[ip][1][dir2] = EVEN;  // w_z on the right
                } else if (ip == dir1) {
                    // the Y bc depends on the symmetry
                    if (arg_sym_y == 1) {
                        valCase.ycntr           = 0.25;
                        valCase.ysign           = 1.0;
                        valCase.mybcv[ip][0][dir0] = EVEN;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = ODD;   // w_y on the left
                        valCase.mybcv[ip][0][dir2] = EVEN;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;   // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;   // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;   // w_z on the right
                    } else if (arg_sym_y == 0) {
                        valCase.ycntr           = 0.5;
                        valCase.ysign           = 0.0;
                        valCase.mybcv[ip][0][dir0] = UNB;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = UNB;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = UNB;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;  // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;  // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;  // w_z on the right
                    } else if (arg_sym_y == -1) {
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
                    if (arg_sym_x == 1) {
                        valCase.xcntr           = 0.25;
                        valCase.xsign           = 1.0;
                        valCase.mybcv[ip][0][dir0] = ODD;   // w_x on the left
                        valCase.mybcv[ip][0][dir1] = EVEN;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = EVEN;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;   // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;   // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;   // w_z on the right
                    } else if (arg_sym_x == 0) {
                        valCase.xcntr           = 0.5;
                        valCase.xsign           = 0.0;
                        valCase.mybcv[ip][0][dir0] = UNB;  // w_x on the left
                        valCase.mybcv[ip][0][dir1] = UNB;  // w_y on the left
                        valCase.mybcv[ip][0][dir2] = UNB;  // w_z on the left
                        valCase.mybcv[ip][1][dir0] = UNB;  // w_x on the right
                        valCase.mybcv[ip][1][dir1] = UNB;  // w_y on the right
                        valCase.mybcv[ip][1][dir2] = UNB;  // w_z on the right
                    } else if (arg_sym_x == -1) {
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
            } else if (arg_sol_t == 1) {  // we have tubes in all 3 directions
                                          // the BC are imposed in the Z direction
                valCase.xcntr = 0.5;
                valCase.ycntr = 0.5;
                valCase.zcntr = 0.5;
                valCase.xsign = 1.0;
                valCase.ysign = 1.0;
                valCase.zsign = 1.0;

                // here x means the normal distance to the face
                valCase.mybcv[ip][0][(ip + 0) % 3] = EVEN;  // w_x on the left
                valCase.mybcv[ip][0][(ip + 1) % 3] = ODD;   // w_y on the left
                valCase.mybcv[ip][0][(ip + 2) % 3] = ODD;   // w_z on the left
                valCase.mybcv[ip][1][(ip + 0) % 3] = EVEN;  // w_x on the right
                valCase.mybcv[ip][1][(ip + 1) % 3] = ODD;   // w_y on the right
                valCase.mybcv[ip][1][(ip + 2) % 3] = ODD;   // w_z on the right

            } else if (arg_sol_t == 2) {  // gaussian blob
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
            
            } else if (arg_sol_t == 3) {  // this is the vortex ring
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
            } else if (arg_sol_t == 4) {        // this is the fully periodic case
                valCase.mybcv[ip][0][0] = PER;  
                valCase.mybcv[ip][0][1] = PER;  
                valCase.mybcv[ip][0][2] = PER;  
                valCase.mybcv[ip][1][0] = PER;  
                valCase.mybcv[ip][1][1] = PER;  
                valCase.mybcv[ip][1][2] = PER;  
            } else if (arg_sol_t == 5) {        // this is the even-even
                valCase.mybcv[ip][0][0] = PER; 
                valCase.mybcv[ip][0][1] = EVEN; 
                valCase.mybcv[ip][0][2] = PER; 
                valCase.mybcv[ip][1][0] = PER; 
                valCase.mybcv[ip][1][1] = EVEN; 
                valCase.mybcv[ip][1][2] = PER; 
            } else if (arg_sol_t == 6) {        // this is the odd-even
                // bc for direction ip
                valCase.mybcv[ip][0][0] = ODD;  // w_x on the left
                valCase.mybcv[ip][0][1] = ODD;  // w_y on the left
                valCase.mybcv[ip][0][2] = ODD;  // w_z on the left
                valCase.mybcv[ip][1][0] = EVEN;  // w_x on the right
                valCase.mybcv[ip][1][1] = EVEN;  // w_y on the right
                valCase.mybcv[ip][1][2] = EVEN;  // w_z on the right
            } 
        }
        // let's gooo
        vtube(valCase, (FLUPS_GreenType) arg_kernel, arg_nsolve, arg_sol_t, (FLUPS_DiffType) arg_order, arg_direction);
    }

    free(size);

    MPI_Finalize();
}
