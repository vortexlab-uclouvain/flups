/**
 * @file main.cpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
 */

#include <cmath>
#include <iostream>

#include "mpi.h"
#include "h3lpr/profiler.hpp"
#include "h3lpr/macros.hpp"
#include "h3lpr/parser.hpp"
#include "h3lpr/ptr.hpp"
#include "flups.h"
#include <iostream>
#include <cstring>
#include "validation_3d.hpp"

using namespace std;

#define validation_calloc(size)                                                    \
    ({                                                                             \
        H3LPR::m_ptr<H3LPR::H3LPR_ALLOC_POSIX, void *, FLUPS_ALIGNMENT> ptr(size); \
                                                                                   \
        ptr();                                                                     \
    })

int main(int argc, char *argv[]) {
    //--------------------------------------------------------------------------
    // Initialize MPI
    int rank;
    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    // int requested = MPI_THREAD_FUNNELED;
    int requested = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        printf("The MPI-provided thread behavior does not match\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //--------------------------------------------------------------------------
    // Parsing 
    //--------------------------------------------------------------------------
    H3LPR::Parser parser(argc, (const char**) argv);

    // Retreive the int values 
    auto arg_center  = parser.GetValue<int>("--center", "Indicate the location of the data: 0=Node-centred, 1=Cell-centred", 1);
    auto arg_nsolve  = parser.GetValue<int>("--nsolve", "the number of solves to perform", 1);
    auto arg_kernel  = parser.GetValue<int>("--kernel", "the Green kernel 0=CHAT2, 1=LGF2, 2=HEJ2, 3=HEJ4, 4=HEJ6, 5=HEJ8, 6=HEJ10, 7=HEJ0", 0);
    auto arg_lda     = parser.GetValue<int>("--lda", "the leading dimension fo array, number of component  (1=scalar, 3=vector)", 1);
    auto arg_nsample = parser.GetValue<int>("--nres", "Nr is the number of higher resolutions that will be tested, with a resolution (R * 2^[0:Nr-1])", 1);
    auto arg_outputdir = parser.GetValue<std::string>("--outdir", "the output directory for the error","./");

    // Retreive the vector value
    auto arg_nprocs = parser.GetValues<int, 3>("--np", "the number of processes in each direction", {1, 1, 1});
    auto arg_nres   = parser.GetValues<int, 3>("--res", "the number of unknowns each direction", {16, 16, 16});
    auto arg_L = parser.GetValues<double, 3>("--dom", "the size of the domain each direction", {1., 1., 1.});
    auto arg_bc = parser.GetValues<int, 6>("--bc",
                                           "Bxl Bxr Byl Byr Bzl Bzr : the boundary conditions in x/y/z on each side l/r. 0=EVEN, 1=ODD, 3=PERiodic, 4=UNBounded", {4, 4, 4, 4, 4, 4});
    auto arg_bcv = parser.GetValues<int, 18>("--bcv",
                                             "3x(Bxl Bxr Byl Byr Bzl Bzr) : the boundary conditions in x/y/z on each side l/r, 3 times for each component. 0=EVEN, 1=ODD, 3=PERiodic, 4=UNBounded", {9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9});
    
    
    parser.Finalize();

    int *size = (int *)validation_calloc((arg_nsample)*3 * sizeof(int));
    if (size == NULL) {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);
    }

    for (int i = 0; i < arg_nsample * 3; i += 3){
        size[i]   = arg_nres[0] * pow(2,i/3);
        size[i+1] = arg_nres[1] * pow(2,i/3);
        size[i+2] = arg_nres[2] * pow(2,i/3);
    }

    //--------------------------------------------------------------------------
    // Do the validation
    //--------------------------------------------------------------------------
    struct DomainDescr valCase;
    valCase.dovectorbc = (arg_bcv[0] != 9); //(NONE = 9 ) -> If there are no vector Boundary condition provided, then we deal with a scalar test case
    
    //..........................................................................
    // Display
    flups_info(argc, argv);
    if (rank == 0) {
        printf("I will run with:\n");
        printf("  --nprocs: %d,%d,%d\n", arg_nprocs[0], arg_nprocs[1], arg_nprocs[2]);
        printf("  -L: %lf,%lf,%lf\n", arg_L[0], arg_L[1], arg_L[2]);
        if (!valCase.dovectorbc) {
            printf("  -bc: %d,%d ; %d,%d ; %d,%d\n", arg_bc[0], arg_bc[1], arg_bc[2], arg_bc[3], arg_bc[4], arg_bc[5]);
        } else {
            printf("  -bcv: [%d,%d ; %d,%d ; %d,%d],[%d,%d ; %d,%d ; %d,%d],[%d,%d ; %d,%d ; %d,%d]\n", arg_bcv[0], arg_bcv[1], arg_bcv[2], arg_bcv[3], arg_bcv[4], arg_bcv[5],
                   arg_bcv[6], arg_bcv[7], arg_bcv[8], arg_bcv[9], arg_bcv[10], arg_bcv[11],
                   arg_bcv[12], arg_bcv[13], arg_bcv[14], arg_bcv[15], arg_bcv[16], arg_bcv[17]);
        }
        printf("  --kernel: %d\n", arg_kernel);
        printf("  --lda: %d\n", arg_lda);
        printf("  --center: %d\n", arg_center);
        printf("  --nsolve: %d\n", arg_nsolve);
        for (int i = 0; i < arg_nsample; i++) {
            printf("   -> sample %d: %d %d %d\n", i + 1, size[i * 3], size[i * 3 + 1], size[i * 3 + 2]);
        }
        printf(" --outdir:%s\n",arg_outputdir.c_str());
    }

    //..........................................................................
    for (int is = 0; is < arg_nsample; is++) {
        for (int ip = 0; ip < 3; ip++) {
            valCase.L[ip]           = arg_L[ip];
            valCase.nproc[ip]       = arg_nprocs[ip];
            valCase.nglob[ip]       = size[is * 3 + ip];
            valCase.center_type[ip] = (FLUPS_CenterType)   arg_center;
            valCase.mybc[ip][0]     = (FLUPS_BoundaryType) arg_bc[2*ip];
            valCase.mybc[ip][1]     = (FLUPS_BoundaryType) arg_bc[2*ip + 1];
            
            valCase.mybcv[ip][0][0] = (FLUPS_BoundaryType) arg_bcv[0 + 2*ip];
            valCase.mybcv[ip][1][0] = (FLUPS_BoundaryType) arg_bcv[0 + 2*ip + 1];
            
            valCase.mybcv[ip][0][1] = (FLUPS_BoundaryType) arg_bcv[6 + 2*ip];
            valCase.mybcv[ip][1][1] = (FLUPS_BoundaryType) arg_bcv[6 + 2*ip + 1];
            
            valCase.mybcv[ip][0][2] = (FLUPS_BoundaryType) arg_bcv[12+ 2*ip];
            valCase.mybcv[ip][1][2] = (FLUPS_BoundaryType) arg_bcv[12+ 2*ip + 1];
        }

        // let's gooo
        validation_3d(valCase, (FLUPS_GreenType) arg_kernel, arg_lda, arg_nsolve, arg_outputdir);
    }

    

    free(size);

    MPI_Finalize();
    //--------------------------------------------------------------------------
}
