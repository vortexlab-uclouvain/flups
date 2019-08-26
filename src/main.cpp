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
#include "tools.hpp"

#include "FFTW_Solver.hpp"
#include "SwitchTopo.hpp"
#include "validation_2d.hpp"
#include "validation_3d.hpp"

#include "mpi.h"


int main(int argc, char *argv[]) {
#if (DIM == 2)
    // int nsample = 2; int size[2] = {64,128};
    int nsample = 3;
    int size[3] = {64, 128, 256};
    // int nsample = 4; int size[4] = {64,128,256,512};

    validation_2d_UU_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_UU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UU_UE(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_UE(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_UU_EU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_EU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UE_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UE_UU(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_EU_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_EU_UU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UU_UO(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_UO(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_UU_OU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_OU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UO_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UO_UU(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_OU_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_OU_UU(nsample, size, UP_SRHS, HEJ_4);
#elif (DIM == 3)

    int rank, is;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // int nsample = 3; int size[3] = {16, 32, 64};
    // int nsample = 4; int size[4] = {16, 32, 64,128};
    // int nsample = 3; int size[3] = {64,128,256};
    int nsample = 1; int size[1] = {16};

    // loop over the resolution for convergence study
    for (int is = 0; is < nsample; is++) {

        struct DomainDescr valCase ;
        valCase.nglob[0] = size[is] * valCase.L[0];
        valCase.nglob[1] = size[is] * valCase.L[1];
        valCase.nglob[2] = size[is] * valCase.L[2];

        valCase.center[0] = .2; valCase.center[1] = .8;

        if(rank==0) printf("\n==============================  UNB UNB - UNB UNB - UNB UNB ============================================\n");
        
        if(rank==0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        validation_3d(valCase,UP_SRHS,HEJ_2);
        // printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_4);
        // printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_6);
        // printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,CHAT_2);

        if(rank==0) printf("\n==============================  EVEN UNB - UNB UNB - UNB UNB ============================================\n");
        valCase.mybc[0][0] = EVEN;

        if(rank==0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        validation_3d(valCase,UP_SRHS,HEJ_2);
        // printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_4);
        // printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_6);
        // printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,CHAT_2);

        if(rank==0) printf("\n==============================  EVEN UNB - UNB ODD - UNB UNB ============================================\n");
        valCase.mybc[1][1] = ODD;

        if(rank==0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        validation_3d(valCase,UP_SRHS,HEJ_2);


    }
    MPI_Finalize();
#endif
}