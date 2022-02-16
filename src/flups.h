/**
 * @file flups.h
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief This is the external API of the FLUPS library
 * @version
 * 
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef FLUPS_H
#define FLUPS_H

#include "mpi.h"

#ifdef __cplusplus
#include <cmath>
#include <iostream>
extern "C" {
#define MAX(a,b) std::max(a,b)
#else
#include "stdlib.h"
#define MAX(a,b) a>b?a:b
#endif

#ifdef PB_DBG
#define DEBUG_ST 1 
#endif
//=============================================================================
/**
 * @name STRUCTURES AND DEFINITIONS
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
    UNB  = 4,  /**< UNBounded boundary condition */
    NONE = 9 /**< No boundary condition = dimension not used */
};

/**
 * @brief The type of Green's function used for the Poisson solver
 * 
 */
enum FLUPS_GreenType {
    CHAT_2 = 0, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2  = 1, /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2  = 2, /**< @brief regularized, order 2, Hejlesen et al. (2015)*/
    HEJ_4  = 3, /**< @brief regularized, order 4, Hejlesen et al. (2015)*/
    HEJ_6  = 4, /**< @brief regularized, order 6, Hejlesen et al. (2015)*/
    HEJ_8  = 5, /**< @brief regularized, order 8, Hejlesen et al. (2015)*/
    HEJ_10 = 6, /**< @brief regularized, order 10, Hejlesen et al. (2015)*/
    HEJ_0  = 7, /**< @brief Fourier cutoff, spectral-like, Hejlesen et al. (2019)*/
};

/**
 * @brief The type of possible solvers
 * 
 * When solving for Biot-Savart, the Green's kernel \f$ G \f$ in Fourier space is adapted so that the Fourier
 * transform of the solution is obtained as
 *      \f[ \hat{\phi} = \hat{K} \times \hat{f} \f]
 * where \f$ \hat{K} \f$ is the spectral equivalent of the gradient of \f$ G \f$. Indeed,
 * the derivation is performed directly in Fourier space, according to the parameter @ref FLUPS_DiffType.
 */
enum FLUPS_SolverType {
    STD = 0, /**< @brief the standard poisson solver: \f$ \nabla^2(\phi) = (rhs) \f$ */
    ROT = 1 /**< @brief the Bio-Savart poisson solver: \f$ \nabla^2(\phi) = \nabla \times (rhs) \f$ */
};

/**
 * @brief The type of derivative to be used with @ref FLUPS_SolverType ROT.
 *  
 */
enum FLUPS_DiffType {
    NOD = 0, /**< @brief Default parameter to be used with the STD type solve */
    SPE = 1, /**< @brief Spectral derivation, \f$ \hat{K} = i \, k \, \hat{G} \f$ */
    FD2 = 2 /**< @brief Spectral equivalent of 2nd order finite difference, \f$ \hat{K} = i \, \sin(k) \, \hat{G} \f$ */
};


/**
 * @brief List of supported data center
 * 
 * The data can be either Cell-centered or Node-centered
 */
enum FLUPS_CenterType {
    NODE_CENTER = 0, /**< NODE centered data  */
    CELL_CENTER = 1 /**< CELL centered data */
};


/**
 * @brief to be used as "sign" for all of the FORWARD tranform
 * 
 */
#define FLUPS_FORWARD -1  // equivalent to FFTW_FORWARD

/**
 * @brief to be used as "sign" for all of the BACKWARD tranform
 * 
 */
#define FLUPS_BACKWARD 1  // equivalen to FFTW_BACKWARD

/**
 * @brief to be used as "sign" for all of the BACKWARD tranform
 * 
 */
#define FLUPS_BACKWARD_DIFF 2

/**
 * @brief Memory alignment in bytes.
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
typedef enum FLUPS_DiffType     FLUPS_DiffType;
typedef enum FLUPS_CenterType   FLUPS_CenterType;

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
void* flups_malloc(size_t size);

/**
 * 
 * @brief Free the memory allocated with flups_malloc
 * 
 * @warning You must free the memory allocate using flups_malloc using this function.
 * 
 * @param data the data to be freed
 */
void flups_free(void* data);

/**
 * @brief compute the memory local index for a point (i0,i1,i2) in axsrc-indexing in a memory.
 * The returned value is in the axtrg-indexing
 * 
 * For example if going through a complex topology following the standard indexing:
 * @code{.cpp}
    // the topology is complex
    const int nf = 2;
    // get the topology indexing
    const int ax0     = flups_topo_get_axis(topo);
    // the memory size is given in the 012 order
    const int nmem[3] = {flups_topo_get_nmem(topo,0),flups_topo_get_nmem(topo,1), flups_topo_get_nmem(topo,2)};
    for (int lia = 0; lia < flups_topo_get_lda(topo); lia++){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,1); i1++) {
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,0); i0++) {
                    // the i0, i1 and i2 are given in a 0-indexing
                    // the id is aimed for an array in the ax0-indexing
                    const size_t id = flups_locID(0, i0, i1, i2, lia, ax0, nmem, nf);
                        
                    data[id+0] = ...;
                    data[id+1] = ...;
                }
            }
        }
    }
 * @endcode
 * 
 * @param axsrc the FRI, reference axis aligned with index i0
 * @param i0 the index in the axsrc direction
 * @param i1 the index in the (axsrc+1)%3 direction
 * @param i2 the index in the (axsrc+2)%3 direction
 * @param lia the index of the vector component (leading index of array)
 * @param axtrg the topology FRI, i.e. the way the memory is aligned in the current topology
 * @param size the size of the memory (given in the 012-order)
 * @param nf the number of unknows in one element
 * @return size_t 
 */
static inline size_t flups_locID(const int axsrc, const int i0, const int i1, const int i2, const int lia, const int axtrg, const int size[3], const int nf) {
    const int i[3] = {i0, i1, i2};
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;
    const int ax2  = (ax0 + 2) % 3;

    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * (i[dax2] + lia * size[ax2]));
}

/**
 * @brief compute them symmetrized local index for a point (i0,i1,i2) in axsrc-indexing in an extended topology (e.g. spectral topologies).
 * The returned value is in the axtrg-indexing.
 * 
 * For example if going through a complex topology following the standard indexing, one can get the spectral indexing:
 * @code{.cpp}
    const int ax0     = flups_topo_get_axis(topoSpec);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nf      = 2; //topo is complex
    
    // get the memory size of the spectral array
    int nmemSpec[3];
    for(int i=0;i<3;i++){
        nmemSpec[i] = flups_topo_get_nmem(topoSpec,i);
    }
        
    for (int lia = 0; lia < flups_topo_get_lda(topoSpec); lia++){
        for (int i2 = 0; i2 < flups_topo_get_nloc(topoSpec,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topoSpec,ax1); i1++) {
                //local indexes start
                const size_t id = flups_locID(ax0, 0, i1, i2, lia, ax0, nmemSpec, nf);
                for (int i0 = 0; i0 < flups_topo_get_nloc(topoSpec,ax0); i0++) {
                    int is[3];
                    // get the symmetrized ID
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
    }
 * @endcode
 * 
 * @param axsrc the FRI, reference axis aligned with index i0
 * @param i0 the index in the axsrc direction
 * @param i1 the index in the (axsrc+1)%3 direction
 * @param i2 the index in the (axsrc+2)%3 direction
 * @param istart start index of the local block (as provided by @ref flups_get_istartGlob)
 * @param symstart indexes where the symmetry starts, i.e. the first index which is symmetrized (as provided by @ref flups_get_spectralInfo)
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
 * @brief Creates and returns a topology.
 * 
 * @warning Once specified, the fastest rotating inde defines the memory layout.
 * We assume fortran memory layout, i.e. if the FRI is 2, the next dimension is 0 and the last one is 1.
 * This is opposed to the C indexing: when the FRI is 2, the next dimension is 1 and the last one is 0.
 * 
 * @param axis The direction which is aligned with the fastest rotating index
 * @param lda leading dimension of the array, i.e. the number of components for a vector field
 * @param nglob The global number of points in each direction of the domain
 * @param nproc The number of processors per direction.
 * @param isComplex The state of the topo: real (false) or complex (true)
 * @param axproc The correspondance between the physical dimensions and the rank decomposition. NULL for the default behavior (0 1 2).
 * @param alignment Memory alignement constant: the memsize are adapted so that . See FLUPS_ALIGNMENT, or by default 
 * @return FLUPS_Topology* pointer to the topology
 */
FLUPS_Topology* flups_topo_new(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm);

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
 * @brief Determines the leading dimension of arrays, i.e. the number of vector components
 * 
 * @param t 
 * @return int
 */
int flups_topo_get_lda(const FLUPS_Topology* t);

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
 * @brief Determines the local number of points in the domain (on this rank) in a given direction
 * 
 * @warning due to some memory padding to ensure memory alignement for the FFTs, @ref flups_topo_get_nloc may
 * not return the same result as @ref flups_topo_get_nmem
 * 
 * @param t 
 * @param dim 
 * @return int the number of elements (real or complex)
 */
int  flups_topo_get_nloc(const FLUPS_Topology* t, const int dim);
/**
 * @brief Determines the local memory size per direction
 * 
 * @warning due to some memory padding to ensure memory alignement for the FFTs, @ref flups_topo_get_nloc may
 * not return the same result as @ref flups_topo_get_nmem
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
 * @brief returns the scalar local size of on this rank, i.e. the number of unknowns in this rank for one dimension
 * 
 * @return long 
 */
size_t flups_topo_get_locsize(const FLUPS_Topology* t);

/**
 * @brief returns the memory size of on this proc, i.e. the number of bytes in this proc, including padded memory and all the vector components
 * 
 * @return long 
 */
size_t flups_topo_get_memsize(const FLUPS_Topology* t);

/**
 * @brief returns the MPI-communicator of the topology
 * 
 * @param t the Topology of interest
 * @param comm the communicator
 */
MPI_Comm flups_topo_get_comm(FLUPS_Topology* t);



/**
 * @brief split the rank into rank per dimensions
 * 
 * axproc is not used if comm is of type MPI_CART.
 * 
 * @param rank the rank of the proc (from MPI, in the current communicator of the topo)
 * @param nproc the number of procs along each direction
 * @param comm a communicator. If it is of type MPI_CART, we use that information in the returned rankd.
 * @param rankd the rank per dimension in XYZ format
 */
static inline void flups_ranksplit(const int rank, const int axproc[3], const int nproc[3], MPI_Comm comm, int rankd[3]);

/**
 * @brief get the rank from the rank per dimension
 * 
 * @param rankd the rank in XYZ format
 * @param topo the topology
 * @return int 
 */
static inline int flups_rankindex(const int rankd[3], const Topology *topo);

/**@} */

//=============================================================================
/**
 * @name SOLVER
 * @{
 */

/**
 * @brief Creates a solver for the specified domain.
 * 
 * @param t user-determined topology of data in physical space, describing the data that will be provided to the solver
 * @param bc boundary conditions of the domain for the right hand side
 * @param h physical space increment in each direction
 * @param L physical length of the domain in each direction
 * @param orderdiff order of the derivatives for ROT solver (SPE = spectral, FD2 = 2nd order final differences). Can be set to NONE if only STD solve are called.
 * @return FLUPS_Solver* the new solver
 */
FLUPS_Solver* flups_init(FLUPS_Topology* t, FLUPS_BoundaryType* bc[3][2], const double h[3], const double L[3], FLUPS_DiffType orderDiff, const FLUPS_CenterType center_type[3]);
/**
 * @brief Same as @ref flups_init, with a profiler for the timing of the code (if compiled with PROF, if not, it will not use the profiler).
 * 
 * @param prof 
 */
FLUPS_Solver* flups_init_timed(FLUPS_Topology* t, FLUPS_BoundaryType* bc[3][2], const double h[3], const double L[3], const FLUPS_DiffType orderDiff, const FLUPS_CenterType center_type[3], FLUPS_Profiler* prof);

/**
 * @brief must be called before execution terminates as it frees the memory used by the solver
 * 
 * @param s 
 */
void flups_cleanup(FLUPS_Solver* s);

/**
 * @brief sets the type of the Green's function used by the solver
 * 
 * @warning must be done before @ref flups_setup
 * 
 * @param s 
 * @param type the type of the solver (CHAT2 by default)
 */
void    flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type);

/**
 * @brief setup the solver and do the memory allocation
 * 
 * @warning after this call the solver cannot been change anymore!
 * 
 * @warning if changeComm is true, you need to update MPI rank based on the new communicator that is provided by @ref flups_topo_get_comm
 * 
 * @param s 
 * @param changeComm indicate if FLUPS is allowed to change the communicator of the Topology used to inilialize the solver (only valid if compiled with RORDER_RANKS)
 * @return double* 
 */
double* flups_setup(FLUPS_Solver* s,const bool changeComm);

/**
 * @brief solve the Poisson equation on rhs, and returns the solution in field (can be done in-place)
 * 
 * FLUPS uses his own allocated memory to do the operations. 
 * 
 * @param s 
 * @param field 
 * @param rhs 
 */
void flups_solve(FLUPS_Solver* s, double* field, double* rhs, const FLUPS_SolverType type);

/**@} */

//=============================================================================
/**
 * @name SOLVER (Advanced)
 * Optimized use of the library may require the use of advanced features. This is especially true if you
 * want to use the memory allocated by flups, and avoid copies.
 * 
 * The functions @ref flups_do_copy, @ref flups_do_FFT, @ref flups_do_mult give access to atomic operations 
 * that are done in @ref flups_solve.
 * @{
 */

/**
 * @brief get the maximun amount of memory required by FLUPS
 * 
 * @param s 
 * @return size_t 
 */
size_t flups_get_allocSize(FLUPS_Solver* s);

/**
 * @brief get information required to compute the spectral mode associated with each spectral field entry
 * 
 * The spectral mode in direction i is given by (index[i] + koffset[i])*kfact[i]
 * 
 * @param s the FLUPS solver
 * @param kfact returns the multiplication factor to used to get 
 * @param koffset returns the spectral offeset given the type of boundary condition used
 * @param symstart the first point which is symmetrized, to use with @ref flups_symID
 */
void flups_get_spectralInfo(FLUPS_Solver* s, double kfact[3], double koffset[3], double symstart[3]);

/**
 * @brief while using regularized Hejlesen kernels, set the alpha factor, i.e. the number of grid points in the smoothing Gaussian
 * Notice: this parameter only affect kernels: HEJ2,HEJ4,HEJ6,HEJ8,HEJ10
 * 
 * @param s 
 * @param alpha (default value is 2.0)
 */
void flups_set_alpha(FLUPS_Solver* s, const double alpha);   //must be done before setup

// /**
//  * @brief sets the order of derivative while using divergence or rotational formulation
//  * 
//  * @param s 
//  * @param order 
//  */
// void flups_set_OrderDiff(FLUPS_Solver* s, const int order);  //must be done before setup

/**
 * @brief returns the physical topology, i.e. the one used for rhs and solution
 * 
 * @param s 
 * @return const FLUPS_Topology* 
 */
const FLUPS_Topology* flups_get_innerTopo_physical(FLUPS_Solver* s);
/**
 * @brief returns the spectral topology, i.e. the one which is fully spectral
 * 
 * @param s 
 * @return const FLUPS_Topology* 
 */
const FLUPS_Topology* flups_get_innerTopo_spectral(FLUPS_Solver* s);

/**
 * @brief do the copy from the data provided by the user to FLUPS owned data arrays
 * 
 * @param s 
 * @param topo 
 * @param data 
 * @param sign 
 */
void flups_do_copy(FLUPS_Solver* s, const FLUPS_Topology* topo, double* data, const int sign);
/**
 * @brief compute the FFT, go from the physical space to the spectral one
 * 
 * @param s 
 * @param data 
 * @param sign 
 */
void flups_do_FFT(FLUPS_Solver* s, double* data, const int sign);
/**
 * @brief compute the multiplication between the Green's function and the field
 * 
 * @param s 
 * @param data 
 * @param type 
 */
void flups_do_mult(FLUPS_Solver* s, double* data, const FLUPS_SolverType type);

/**@} */


//=============================================================================
/**
 * @name PROFILER - TIMERS
 * @{
 */

/**
 * @brief create a timer using the default name "default".
 * 
 * @return FLUPS_Profiler* 
 */
FLUPS_Profiler* flups_profiler_new();
/**
 * @brief create a timer with a name "name"
 * 
 * @param name 
 * @return FLUPS_Profiler* 
 */
FLUPS_Profiler* flups_profiler_new_n(const char name[]);
/**
 * @brief free the profiler created
 * 
 * @param p 
 */
void            flups_profiler_free(FLUPS_Profiler* p);
/**
 * @brief display the profiler using the "root" as a reference
 * 
 * @param p 
 */
void            flups_profiler_disp_root(FLUPS_Profiler* p);
/**
 * @brief display the profiler using "name" as reference
 * 
 * @param p 
 * @param name 
 */
void            flups_profiler_disp(FLUPS_Profiler* p,const char name[]);

/**@} */

//=============================================================================
/**
 * @name HDF5 exports
 * @{
 */

/**
 * @brief dumps data in a hdf5 format with a xdmf compatible description, in the folder ./data
 * 
 * @param topo 
 * @param filename 
 * @param data 
 */
void flups_hdf5_dump(const FLUPS_Topology *topo, const char filename[], const double *data);

/**@} */

void flups_print_data(const FLUPS_Topology *topo, double* data);


#ifdef __cplusplus
}
#endif

#endif
