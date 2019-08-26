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
#ifndef VALIDATION_2D_HPP
#define VALIDATION_2D_HPP

#include "defines.hpp"
#include "FFTW_Solver.hpp"
#include "expint.hpp"

/**
 * @name Fully unbounded
 * 
 */
/**@{ */
void validation_2d_UU_UU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
/**@} */

/**
 * @name Unbounded Even
 * 
 */
/**@{ */
void validation_2d_UU_UE(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
void validation_2d_UU_EU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
void validation_2d_EU_UU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
void validation_2d_UE_UU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
/**@} */

/**
 * @name Unbounded Even
 * 
 */
/**@{ */
void validation_2d_UU_UO(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
void validation_2d_UU_OU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
void validation_2d_OU_UU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
void validation_2d_UO_UU(const int nsample, const int* size, const SolverType type, const GreenType GreenType);
/**@} */

#endif