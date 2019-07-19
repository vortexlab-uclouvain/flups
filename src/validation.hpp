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
#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include "defines.hpp"
#include "FFTW_Solver.hpp"
#include "expint.hpp"

/**
 * @name Fully unbounded
 * 
 */
/**@{ */
void validation_2d_UU_UU(const SolverType type, const OrderDiff orderdiff);
/**@} */

/**
 * @name Unbounded Even
 * 
 */
/**@{ */
void validation_2d_UU_UE(const SolverType type, const OrderDiff orderdiff);;
void validation_2d_UU_EU(const SolverType type, const OrderDiff orderdiff);;
void validation_2d_EU_UU(const SolverType type, const OrderDiff orderdiff);;
void validation_2d_UE_UU(const SolverType type, const OrderDiff orderdiff);;
/**@} */

/**
 * @name Unbounded Even
 * 
 */
/**@{ */
void validation_2d_UU_UO(const SolverType type, const OrderDiff orderdiff);;
void validation_2d_UU_OU(const SolverType type, const OrderDiff orderdiff);;
void validation_2d_OU_UU(const SolverType type, const OrderDiff orderdiff);;
void validation_2d_UO_UU(const SolverType type, const OrderDiff orderdiff);;
/**@} */

#endif