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
    LGF_2  = 1,  /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2  = 2,  /**< @brief regularized in zero, order 2, Hejlesen et al. (2015)*/
    HEJ_4  = 3,  /**< @brief regularized in zero, order 4, Hejlesen et al. (2015)*/
    HEJ_6  = 4,  /**< @brief regularized in zero, order 6, Hejlesen et al. (2015)*/
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

typedef struct Solver FLUPS_solver;
typedef struct Topology FLUPS_topo;


/** *********************************************************************
 * @name TOPOLOGIES
 * @{
 ********************************************************************* */

#define flups_new_topo(_1,_2,_3,_4,_5) flups_new_topo_(__VA_ARGS__) 
FLUPS_topo* flups_new_topo_(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3]);
#define flups_new_topo(_1,_2,_3,_4,_5,_6) flups_new_topo_i(__VA_ARGS__) 
FLUPS_topo* flups_new_topo_a(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment);
void        flups_free_topo(FLUPS_topo* t);

int flups_topo_axis(FLUPS_topo* t);
// int flups_topo_nf() const ;
bool flups_topo_isComplex(FLUPS_topo* t);
// int flups_topo_switch_toComplex(FLUPS_topo* t);
// int flups_topo_switch_toReal(FLUPS_topo* t);
int flups_topo_get_nglob(FLUPS_topo* t, const int dim);
int flups_topo_get_nloc(FLUPS_topo* t, const int dim);
int flups_topo_get_nmem(FLUPS_topo* t, const int dim);
int flups_topo_get_nproc(FLUPS_topo* t, const int dim);
// int flups_topo_get_rankd(FLUPS_topo* t,const int dim) ;
// int flups_topo_get_nbyproc(FLUPS_topo* t,const int dim);
// int flups_topo_get_axproc(FLUPS_topo* t,const int dim) ;
void flups_topo_get_istartGlob(FLUPS_topo* t, int istart[3]);

/**
 * @brief returns the local size of on this proc
 * 
 * @return long 
 */

long flups_topo_locsize(FLUPS_topo* t);
/**
 * @brief returns the memory size of on this proc
 * 
 * @return long 
 */
long flups_topo_memsize(FLUPS_topo* t);



/**@} *****************************************************************/


/** *********************************************************************
 * @name SOLVER
 * @{
 ********************************************************************* */

// get a new solver
FLUPS_solver* flups_new_solver(FLUPS_topo* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]);
#ifdef PROF
FLUPS_solver* flups_new_solver_timed(FLUPS_topo* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]);
#endif

// destroy the solver
void flups_free_solver(FLUPS_solver* s);

// setup the solver
void flups_set_greenType(FLUPS_solver* s, const FLUPS_GreenType type);
void flups_setup(FLUPS_solver* s);

// solve
void flups_solve(FLUPS_solver* s, const FLUPS_topo* tIn, const FLUPS_topo* tOut);


// -- ADVANCED FEATURES --

void flups_set_alpha(const double alpha); //must be done before setup
void flups_set_OrderDiff(const int order); //must be done before setup

// void flups_do_FFTfwd(FLUPS_solver* s);
// void flups_do_FFTbck(FLUPS_solver* s);
// void flups_do_mult(FLUPS_solver* s);

/**@} *****************************************************************/


// #ifdef PROF
// //**********************************************************************
// //  PROFILER - TIMERS
// //**********************************************************************

// typedef struct Profiler Profiler;

// Profiler* profiler_new();
// Profiler* profiler_new(char* name);
// void    profiler_free(Profiler* p);
// void profiler_disp();
// #endif


#ifdef __cplusplus
}
#endif