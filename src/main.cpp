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
    validation_2d_UU_UU(UP_SRHS,HEJ_2);
    validation_2d_UU_UU(UP_SRHS,HEJ_4);

    validation_2d_UU_UE(UP_SRHS,HEJ_2);
    validation_2d_UU_UE(UP_SRHS,HEJ_4);

    validation_2d_UU_UO(UP_SRHS,HEJ_2);
    validation_2d_UU_UO(UP_SRHS,HEJ_4);
    validation_2d_UU_OU(UP_SRHS,HEJ_4);
}