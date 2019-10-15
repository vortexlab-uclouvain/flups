/**
 * @file flups.h
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief This is the external API of the FLUPS library
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

#ifndef FLUPS_H
#define FLUPS_H


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief type to descive the number of data per proc.
 * 
 * If your application is C++, you can replace this by size_t.
 * 
 */
#define FLUPS_SIZE unsigned long long

//=============================================================================
/**
 * @name Common definitions
 * @{
 */
//=============================================================================


/**
 * @brief List of supported boundary conditions
 * 
 * The boundary condition can be EVEN, ODD, PERiodic or UNBounded.
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
    VRHS, /**<@brief vector \f$ \nabla^2 f = rhs \f$ */
    ROT,  /**<@brief vector \f$ \nabla^2 f = \nabla \times rhs \f$ */
    DIV   /**<@brief scalar \f$ \nabla^2 f = \nabla \cdot rhs \f$ */
};

/**
 * @brief to be used as "sign" for all of the FORARD tranform
 * 
 */
#define FLUPS_FORWARD -1  // = FFTW_FORWARD

/**
 * @brief to be used as "sign" for all of the BACKWARD tranform
 * 
 */
#define FLUPS_BACKWARD 1  // = FFTW_BACKWARD

/**
 * @brief Memory alignment constant in bytes.
 * 
 */
#define FLUPS_ALIGNMENT 16

/**
 * @brief FFTW planner flag
 * 
 */
#define FFTW_FLAG FFTW_PATIENT

typedef struct Solver   FLUPS_Solver;
typedef struct Topology FLUPS_Topology;
typedef struct Profiler FLUPS_Profiler;

typedef enum FLUPS_BoundaryType FLUPS_BoundaryType;
typedef enum FLUPS_GreenType    FLUPS_GreenType;
typedef enum FLUPS_SolverType   FLUPS_SolverType;

/**@} */

//=============================================================================
/**
 * @name MEMORY MANAGEMENT
 * @{
 */
//=============================================================================

/**
 * 
 * @brief Allocate the memory aligned on FLUPS_ALIGNMENT
 * 
 * @param size the data to be allocated
 */
void * flups_malloc(FLUPS_SIZE size);

/**
 * 
 * @brief Free the memory allocated with flups_malloc
 * 
 * @param data the data to be freed
 */
void flups_free(void* data);

/**
 * @brief compute the memory local index for a point (i0,i1,i2) in axsrc-indexing in a memory.
 * The returned value is in the axtrg-indexing
 * 
 * For example if going through a topology following the standard indexing:
 * @code{.cpp}
 *  const int ax0     = flups_topo_get_axis(topo);
    const int nmem[3] = {flups_topo_get_nmem(topo,0),flups_topo_get_nmem(topo,1), flups_topo_get_nmem(topo,2)};
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo,2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo,1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo,0); i0++) {
                const FLUPS_SIZE id = flups_locID(0, i0, i1, i2, ax0, nmem, 1);
 *                  
 *              data[id] = ...;
 *          }
 *      }
 *  }
 * @endcode
 * 
 * @param axsrc the FRI for the point (i0,i1,i2)
 * @param i0 the index in the axsrc direction
 * @param i1 the index in the (axsrc+1)%3 direction
 * @param i2 the index in the (axsrc+2)%3 direction
 * @param axtrg the topology FRI, i.e. the way the memory is aligned in the current topology
 * @param size the size of the memory (012-indexing)
 * @param nf the number of unknows in one element
 * @return FLUPS_SIZE 
 */
static inline FLUPS_SIZE flups_locID(const int axsrc, const int i0, const int i1, const int i2, const int axtrg, const int size[3], const int nf) {
    const int i[3] = {i0, i1, i2};
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;

    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * i[dax2]);
}

/**@} */

//=============================================================================
/**
 * @name TOPOLOGIES
 * @{
 */
//=============================================================================

/**
 * @brief Create and returns a topology.
 * 
 * @param axis The direction which is aligned with the fastest rotating index
 * @param nglob The global number of points in each direction of the domain
 * @param nproc The number of processors per direction.
 * @param isComplex The state of the topo: real (false) or complex (true)
 * @param axproc The correspondance between the physical dimensions and the dimensions in memory, otherwise NULL.
 * @param alignment Memory alignement constant: the memsize are adapted so that . See FLUPS_ALIGNMENT, or by default 
 * @return FLUPS_Topology* pointer to the topology
 */
const FLUPS_Topology* flups_topo_new(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment);

/**
 * @brief Clean and free the topo.
 * 
 * @param t topo to be freed
 */
void flups_topo_free(const FLUPS_Topology* t);

/**
 * @brief Determines if the topo works on real or complex numbers
 * 
 * @param t 
 * @return true if the topo is on complex numbers
 * @return false if the topo is on real numbers
 */
bool flups_topo_get_isComplex(const FLUPS_Topology* t);

/**
 * @brief Determines the physical direction aligned in memory
 * 
 * @param t 
 * @return int 
 */
int  flups_topo_get_axis(const FLUPS_Topology* t);
int  flups_topo_get_nglob(const FLUPS_Topology* t, const int dim);
int  flups_topo_get_nloc(const FLUPS_Topology* t, const int dim);
int  flups_topo_get_nmem(const FLUPS_Topology* t, const int dim);
int  flups_topo_get_nproc(const FLUPS_Topology* t, const int dim);
void flups_topo_get_istartGlob(const FLUPS_Topology* t, int istart[3]);

/**
 * @brief returns the local size of on this proc
 * 
 * @return long 
 */

FLUPS_SIZE flups_topo_get_locsize(const FLUPS_Topology* t);
/**
 * @brief returns the memory size of on this proc
 * 
 * @return long 
 */
FLUPS_SIZE flups_topo_get_memsize(const FLUPS_Topology* t);

/**@} */

//=============================================================================
/**
 * @name SOLVER
 * @{
 */

// get a new solver
FLUPS_Solver* flups_init(const FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]);
FLUPS_Solver* flups_init_timed(const  FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3],FLUPS_Profiler* prof);

// destroy the solver
void flups_cleanup(FLUPS_Solver* s);

// setup the solver
void    flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type);
double* flups_setup(FLUPS_Solver* s);

// solve
void flups_solve(FLUPS_Solver* s, double* field, double* rhs, const FLUPS_SolverType type);

/**@} */

//=============================================================================
/**
 * @name SOLVER (Advanced)
 * @{
 */
FLUPS_SIZE flups_get_allocSize(FLUPS_Solver* s);

void flups_set_alpha(FLUPS_Solver* s, const double alpha);   //must be done before setup
void flups_set_OrderDiff(FLUPS_Solver* s, const int order);  //must be done before setup

const FLUPS_Topology* flups_get_topo_physical(FLUPS_Solver* s);
const FLUPS_Topology* flups_get_topo_spectral(FLUPS_Solver* s);

void flups_do_copy(FLUPS_Solver* s, const FLUPS_Topology* topo, double* data, const int sign);
void flups_do_FFT(FLUPS_Solver* s, double* data, const int sign);
void flups_do_mult(FLUPS_Solver* s, double* data, const FLUPS_SolverType type);

/**@} */


//=============================================================================
/**
 * @name PROFILER - TIMERS
 * @{
 */

FLUPS_Profiler* flups_profiler_new();
FLUPS_Profiler* flups_profiler_new_n(const char name[]);
void            flups_profiler_free(FLUPS_Profiler* p);
void            flups_profiler_disp(FLUPS_Profiler* p);
void            flups_profiler_disp_root(FLUPS_Profiler* p,const char name[]);

/**@} */

//=============================================================================
/**
 * @name HDF5 exports
 * @{
 */

void flups_hdf5_dump(const FLUPS_Topology *topo, const char filename[], const double *data);

/**@} */

#ifdef __cplusplus
}
#endif

#endif