/**
 * @file validation.hpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-19
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#ifndef VALIDATION_3D_HPP
#define VALIDATION_3D_HPP

#include "FFTW_Solver.hpp"
#include "defines.hpp"
#include "expint.hpp"

#include "hdf5_io.hpp"
#include "mpi.h"

/**
 * @brief everything required to characterize the computational domain and initial condition 
 * 
 */
struct DomainDescr {
    int          nglob[3]   = {64, 64, 64};
    int          nproc[3]   = {2, 2, 1};
    double       L[3]       = {1.0, 1.0, 1.0};
    double       sigma      = 0.1;              // smoothing length scale of the blob
    double       center[3]  = {0.5, 0.5, 0.5};  //center of the blob (in % of L)
    BoundaryType mybc[3][2] = {{UNB, UNB}, {UNB, UNB}, {UNB, UNB}};
};

/**
 * @name validation of the solver using a gaussian blob
 * 
 */
/**@{ */
void validation_3d(const DomainDescr myCase, const SolverType type, const GreenType typeGreen);
/**@} */

#endif