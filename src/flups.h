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

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#define MAX(a,b) std::max(a,b)
#else
#include "stdlib.h"
#define MAX(a,b) a>b?a:b;
#endif

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
void * flups_malloc(size_t size);

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
                const size_t id = flups_locID(0, i0, i1, i2, ax0, nmem, 1);
 *                  
 *              data[id] = ...;
 *          }
 *      }
 *  }
 * @endcode
 * 
 * @param axsrc the FRI, reference axis aligned with index i0
 * @param i0 the index in the axsrc direction
 * @param i1 the index in the (axsrc+1)%3 direction
 * @param i2 the index in the (axsrc+2)%3 direction
 * @param axtrg the topology FRI, i.e. the way the memory is aligned in the current topology
 * @param size the size of the memory (012-indexing)
 * @param nf the number of unknows in one element
 * @return size_t 
 */
static inline size_t flups_locID(const int axsrc, const int i0, const int i1, const int i2, const int axtrg, const int size[3], const int nf) {
    const int i[3] = {i0, i1, i2};
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;

    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * i[dax2]);
}

/**
 * @brief compute the local k-index in spectral coordinates for a point (i0,i1,i2) in axsrc-indexing.
 * The returned value is in the axtrg-indexing.
 * 
 * For example if going through a topology following the standard indexing:
 * @code{.cpp}
 *  const int ax0     = flups_topo_get_axis(topoSpec);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nf      = 2; //topo is complex
    
    int nmemSpec[3];
    for(int i=0;i<3;i++){
        nmemSpec[i] = flups_topo_get_nmem(topoSpec,i);
    }
        
    for (int i2 = 0; i2 < flups_topo_get_nloc(topoSpec,ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topoSpec,ax1); i1++) {
            //local indexes start
            const size_t id = flups_locID(ax0, 0, i1, i2, ax0, nmemSpec,nf);
            for (int i0 = 0; i0 < flups_topo_get_nloc(topoSpec,ax0); i0++) {
                int is[3];
                flups_symID(ax0, i0, i1, i2, istartSpec, symstart, 0, is);

                // the (symmetrized) wave numbers:
                const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];

                data[id + i0 * nf] = ...; //REAL part
                data[id + i0 * nf + 1] = ...; //COMPLEX part
            }
        }
    }
 * @endcode
 * 
 * @param axsrc the FRI, reference axis aligned with index i0
 * @param i0 the index in the axsrc direction
 * @param i1 the index in the (axsrc+1)%3 direction
 * @param i2 the index in the (axsrc+2)%3 direction
 * @param istart start index of the local block (as provided by @flups_get_istartGlob)
 * @param symstart indexes where the symmetry starts (as provided by @flups_get_spectralInfo)
 * @param axtrg the FRI of the target topology, i.e. the way the memory is aligned in the current topology
 * @param is the spectral index
 */
static inline void flups_symID(const int axsrc, const int i0, const int i1, const int i2, const int istart[3], const double symstart[3], const int axtrg, int is[3]) {
    // get the global indexes in the axsrc configuration
    const int ie[3] = {(istart[axsrc] + i0), (istart[(axsrc + 1) % 3] + i1), (istart[(axsrc + 2) % 3] + i2)};
    // cmpt the shift in axis and the axis for the symstart
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;
    const int ax2  = (ax0 + 2) % 3;
    // fill the array in the axtrg configuration
    is[0] = (symstart[ax0] == 0.0 || ie[dax0] <= symstart[ax0]) ? ie[dax0] : MAX((int)fabs(2.0 * symstart[ax0] - ie[dax0]), 1);
    is[1] = (symstart[ax1] == 0.0 || ie[dax1] <= symstart[ax1]) ? ie[dax1] : MAX((int)fabs(2.0 * symstart[ax1] - ie[dax1]), 1);
    is[2] = (symstart[ax2] == 0.0 || ie[dax2] <= symstart[ax2]) ? ie[dax2] : MAX((int)fabs(2.0 * symstart[ax2] - ie[dax2]), 1);
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
 * @param axproc The correspondance between the physical dimensions and the rank decomposition. NULL for the default behavior (0 1 2).
 * @param alignment Memory alignement constant: the memsize are adapted so that . See FLUPS_ALIGNMENT, or by default 
 * @return FLUPS_Topology* pointer to the topology
 */
FLUPS_Topology* flups_topo_new(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm);

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
 * @return int the direction [0-3]
 */
int  flups_topo_get_axis(const FLUPS_Topology* t);
/**
 * @brief Determines the total number of points in the domain in a given direction
 * 
 * @param t 
 * @param dim 
 * @return int the number of elements (real or complex)
 */
int  flups_topo_get_nglob(const FLUPS_Topology* t, const int dim);
/**
 * @brief Determines the local number of points in the domain (on this process) in a given direction
 * 
 * @param t 
 * @param dim 
 * @return int the number of elements (real or complex)
 */
int  flups_topo_get_nloc(const FLUPS_Topology* t, const int dim);
/**
 * @brief Determines the local memory usage per direction
 * 
 * @param t 
 * @param dim 
 * @return int the memory occupation in double/float
 */
int  flups_topo_get_nmem(const FLUPS_Topology* t, const int dim);
/**
 * @brief Determines the number of processes in a given direction
 * 
 * @param t 
 * @param dim 
 * @return int 
 */
int  flups_topo_get_nproc(const FLUPS_Topology* t, const int dim);
/**
 * @brief Determines the start index of this process in all 3 directions
 * 
 * @param t 
 * @param istart 
 */
void flups_topo_get_istartGlob(const FLUPS_Topology* t, int istart[3]);

/**
 * @brief returns the local size of on this proc
 * 
 * @return long 
 */
size_t flups_topo_get_locsize(const FLUPS_Topology* t);

/**
 * @brief returns the memory size of on this proc
 * 
 * @return long 
 */
size_t flups_topo_get_memsize(const FLUPS_Topology* t);

/**
 * @brief returns the communicator of the topology
 * 
 * @param t the Topology of interest
 * @param comm the communicator
 */
MPI_Comm flups_topo_get_comm(FLUPS_Topology* t);

/**@} */

//=============================================================================
/**
 * @name SOLVER
 * @{
 */

// get a new solver
FLUPS_Solver* flups_init(FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]);
FLUPS_Solver* flups_init_timed(FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3],FLUPS_Profiler* prof);

// destroy the solver
void flups_cleanup(FLUPS_Solver* s);

// setup the solver
void    flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type);

/**
 * @brief setup the solver
 * 
 * @warning after this call the solver cannot change anymore!
 * 
 * @warning if changeComm is true, you need to update MPI rank based on the new communicator that is provided by @flups_topo_get_comm 
 * 
 * @param s 
 * @param changeComm indicate if FLUPS is allowed to change the communicator of the Topology used to initialize the solver (only if compiled with RORDER_RANKS)
 * @return double* 
 */
double* flups_setup(FLUPS_Solver* s,const bool changeComm);

// solve
void flups_solve(FLUPS_Solver* s, double* field, double* rhs, const FLUPS_SolverType type);

/**@} */

//=============================================================================
/**
 * @name SOLVER (Advanced)
 * @{
 */
size_t flups_get_allocSize(FLUPS_Solver* s);

void flups_get_spectralInfo(FLUPS_Solver* s, double kfact[3], double koffset[3], double symstart[3]);

void flups_set_alpha(FLUPS_Solver* s, const double alpha);   //must be done before setup
void flups_set_OrderDiff(FLUPS_Solver* s, const int order);  //must be done before setup

const FLUPS_Topology* flups_get_innerTopo_physical(FLUPS_Solver* s);
const FLUPS_Topology* flups_get_innerTopo_spectral(FLUPS_Solver* s);

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
void            flups_profiler_disp_root(FLUPS_Profiler* p);
void            flups_profiler_disp(FLUPS_Profiler* p,const char name[]);

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