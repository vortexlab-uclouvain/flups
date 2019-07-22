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
#include "validation.hpp"


int main(int argc, char* argv[])
{

    int nsample = 4;
    int size[4] = {64,128,256,512};

    validation_2d_UU_UU(nsample,size,UP_SRHS,HEJ_2);
    validation_2d_UU_UU(nsample,size,UP_SRHS,HEJ_4);

    validation_2d_UU_UE(nsample,size,UP_SRHS,HEJ_2);
    validation_2d_UU_UE(nsample,size,UP_SRHS,HEJ_4);

    validation_2d_UU_UO(nsample,size,UP_SRHS,HEJ_2);
    validation_2d_UU_UO(nsample,size,UP_SRHS,HEJ_4);
    validation_2d_UU_OU(nsample,size,UP_SRHS,HEJ_4);
}