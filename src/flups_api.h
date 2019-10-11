/**
 * @file flups_api.h
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-10-09
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief the boundary condition can be EVEN, ODD, PERiodic or UNBounded
 * 
 */
enum FLUPS_BoundaryType {
    EVEN = 0, /**< EVEN boundary condition = zero flux  */
    ODD  = 1, /**< ODD boundary condition = zero value */
    PER  = 3, /**< PERiodic boundary conditions */
    UNB  = 4  /**< UNBounded boundary condition */
};

/**
 * @brief The type of Green's function used for the Poisson solver
 * 
 */
enum FLUPS_GreenType {
    CHAT_2 = 0, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2  = 1, /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2  = 2, /**< @brief regularized in zero, order 2, Hejlesen et al. (2015)*/
    HEJ_4  = 3, /**< @brief regularized in zero, order 4, Hejlesen et al. (2015)*/
    HEJ_6  = 4, /**< @brief regularized in zero, order 6, Hejlesen et al. (2015)*/
};

/**
 * @brief Type of Poisson equation solved
 * 
 */
enum FLUPS_SolverType {
    SRHS, /**<@brief scalar \f$ \nabla^2 f = rhs \f$ */
    VRHS, /**<@brief vectorial \f$ \nabla^2 f = rhs \f$ */
    ROT,  /**<@brief vectorial \f$ \nabla^2 f = \nabla \times rhs \f$ */
    DIV   /**<@brief scalar \f$ \nabla^2 f = \nabla \cdot rhs \f$ */
};

#define FLUPS_FORWARD -1  // = FFTW_FORWARD
#define FLUPS_BACKWARD 1  // = FFTW_BACKWARD

typedef struct Solver   FLUPS_Solver;
typedef struct Topology FLUPS_Topology;
typedef struct Profiler FLUPS_Profiler;

/** *********************************************************************
 * @name TOPOLOGIES
 * @{
 ********************************************************************* */

FLUPS_Topology* flups_new_topo(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment);
void            flups_free_topo(FLUPS_Topology* t);

bool flups_topo_get_isComplex(FLUPS_Topology* t);
int  flups_topo_get_axis(FLUPS_Topology* t);
int  flups_topo_get_nglob(FLUPS_Topology* t, const int dim);
int  flups_topo_get_nloc(FLUPS_Topology* t, const int dim);
int  flups_topo_get_nmem(FLUPS_Topology* t, const int dim);
int  flups_topo_get_nproc(FLUPS_Topology* t, const int dim);
void flups_topo_get_istartGlob(FLUPS_Topology* t, int istart[3]);

/**
 * @brief returns the local size of on this proc
 * 
 * @return long 
 */

unsigned long long flups_topo_get_locsize(FLUPS_Topology* t);
/**
 * @brief returns the memory size of on this proc
 * 
 * @return long 
 */
unsigned long long flups_topo_get_memsize(FLUPS_Topology* t);



/**@} *****************************************************************/


/** *********************************************************************
 * @name SOLVER
 * @{
 ********************************************************************* */

// get a new solver
#ifndef PROF
FLUPS_Solver* flups_new_solver(FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]);
#else
FLUPS_Solver* flups_new_solver_timed(FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3],Profiler* prof);
#endif

// destroy the solver
void flups_free_solver(FLUPS_Solver* s);

// setup the solver
void flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type);
void flups_setup(FLUPS_Solver* s);

// solve
//topo may be different from the one used to create the solver, but must be compatible !
void flups_solve(FLUPS_Solver* s, const FLUPS_Topology* t, double* field, double* rhs, const FLUPS_SolverType type);

/**@} *****************************************************************/

/** *********************************************************************
 * @name SOLVER (Advanced)
 * @{
 ********************************************************************* */

void flups_set_alpha(FLUPS_Solver* s, const double alpha);   //must be done before setup
void flups_set_OrderDiff(FLUPS_Solver* s, const int order);  //must be done before setup

FLUPS_Topology* flups_get_topo_physical(FLUPS_Solver* s);
FLUPS_Topology* flups_get_topo_spectral(FLUPS_Solver* s);

void flups_do_copy(FLUPS_Solver* s, const FLUPS_Topology* topo, double* data, const int sign);
void flups_do_FFT(FLUPS_Solver* s, double* data, const int sign);
void flups_do_mult(FLUPS_Solver* s, double* data, const FLUPS_SolverType type);

/**@} *****************************************************************/


//**********************************************************************
//  PROFILER - TIMERS
//**********************************************************************

FLUPS_Profiler* profiler_new();
FLUPS_Profiler* profiler_new_n(char name[]);
void            profiler_free(FLUPS_Profiler* p);
void            profiler_disp();

#ifdef __cplusplus
}
#endif