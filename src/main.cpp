/**
 * @file main.cpp
 * @author Thomas Gillis
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
#include "validation_2d.hpp"
#include "validation_3d.hpp"
#include "SwitchTopo.hpp"

#include "mpi.h"

typedef double pos[4];

int main(int argc, char *argv[])
{

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

    MPI_Init(&argc,&argv);
    // int nsample = 3; int size[3] = {16, 32, 64};
    int nsample = 4; int size[4] = {16, 32, 64,128};
    // int nsample = 3; int size[3] = {64,128,256};
    // int nsample = 1; int size[1] = {8};
    validation_3d_UU_UU(nsample,size,UP_SRHS,HEJ_2);
    validation_3d_UU_UU(nsample,size,UP_SRHS,HEJ_4);

    MPI_Finalize();
#endif
}