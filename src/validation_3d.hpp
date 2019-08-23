/**
 * @file validation.hpp
 * @author Thomas Gillis
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
 * @name Fully unbounded
 * 
 */
/**@{ */
void validation_3d_UU_UU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
/**@} */

/**
 * @name Unbounded Even
 * 
 */
/**@{ */
void validation_3d_UU_UE_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
void validation_3d_UU_EU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
void validation_3d_EU_UU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
void validation_3d_UE_UU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
/**@} */

/**
 * @name Unbounded Even
 * 
 */
/**@{ */
void validation_3d_UU_UO_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
void validation_3d_UU_OU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
void validation_3d_OU_UU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
void validation_3d_UO_UU_UU(const int nsample, const int* size, const SolverType type, const GreenType typeGreen);
/**@} */

#endif