/**
 * @file Solver.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
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

#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <cstring>
#include <map>
#include "FFTW_plan_dim.hpp"
#include "defines.hpp"
#include "green_functions.hpp"
#include "hdf5_io.hpp"

#include "SwitchTopo.hpp"
#include "SwitchTopo_a2a.hpp"
#include "SwitchTopo_nb.hpp"

#include "Profiler.hpp"
#include "omp.h"

#ifdef HAVE_METIS
#include "metis.h"
#endif

using namespace std;

/**
 * @brief The Poisson solver
 * 
 * A collection of 3 FFTW_plan_dim for the forward and backward FFT transform of
 * data, plus the transformed required for the Green's function to solve the Poisson equation.
 * The tranforms are done in-place in each direction successively. Between each transform, data
 * are remapped (i.e. transposed) in order to have the memory aligned with the direction of the
 * transform. This is done using SwitchTopo which changes the layout of data between 2 topos.
 * 
 * @warning
 * The memory alignement follows the rules explained on the mainpage.
 * Yet, a transposition of the data is required to perfom the transfroms in the correct order.
 * For an element (ix,iy,iz) its tranposed location is computed as\code{.cpp}
 * const size_t id_transposed = iz*dim_multfact_[2] + iy*dim_multfact_[1] + ix*dim_multfact_[0] + offset_;
 * \endcode
 *  
 * 
 */
class Solver {
   protected:
    int            lda_           = 1;     /**@brief the number of components of the problem, i.e. 2D or 3D */
    int            ndim_          = 3;     /**@brief the dimension of the problem, i.e. 2D or 3D */
    int            fftwalignment_ = 0;     /**< @brief alignement assumed by the FFTW Solver  */
    FLUPS_DiffType odiff_         = NOD;  /**< @brief the order of derivative (spectral = SPE, 2nd order FD = FD2) */
    double         normfact_      = 1.0;   /**< @brief normalization factor so that the forward/backward FFT gives output = input */
    double         volfact_       = 1.0;   /**< @brief volume factor due to the convolution computation */
    double         hgrid_[3]      = {0.0}; /**< @brief grid spacing in the tranposed directions */
    double*        data_          = NULL;  /**< @brief data pointer to the transposed memory */

#if DEBUG_ST
    double*        datad_[4]      = {NULL, NULL, NULL, NULL};  /**< @brief data pointer to the transposed memory */
    void           copydata_(const Topology * const topo, double *datad);
    void           comparedata_(const Topology * const topo, double *datad);
#endif
    

    /**
     * @name Forward and backward 
     * transforms related objects
     */
    /**@{ */
    FFTW_plan_dim* plan_forward_[3];       /**< @brief map containing the plans for the forward fft transforms */
    FFTW_plan_dim* plan_backward_[3];      /**< @brief map containing the plans for the backward fft transforms */
    FFTW_plan_dim* plan_backward_diff_[3]; /**< @brief map containing the plans for the backward fft transforms */

    Topology*      topo_phys_     = NULL;
    Topology*      topo_hat_[3]   = {NULL, NULL, NULL}; /**< @brief map containing the topologies (i.e. data memory layout) corresponding to each transform */
    SwitchTopo*    switchtopo_[3] = {NULL, NULL, NULL}; /**< @brief switcher of topologies for the forward transform (phys->topo[0], topo[0]->topo[1], topo[1]->topo[2]).*/
    opt_double_ptr sendBuf_       = NULL;               /**<@brief The send buffer for switchtopo_ */
    opt_double_ptr recvBuf_       = NULL;               /**<@brief The recv buffer for switchtopo_ */
    /**@} */

    /**
     * @name Green's function (and corresponding forward transform) related vars and objects
     * 
     */
    /**@{ */
    double    alphaGreen_ = 2.0;    /**< @brief regularization parameter for HEJ_* Green's functions */
    double*   green_      = NULL;   /**< @brief data pointer to the transposed memory for Green */
    GreenType typeGreen_  = CHAT_2; /**< @brief the type of Green's function */

    FFTW_plan_dim* plan_green_[3];                            /**< @brief map containing the plan for the Green's function */
    Topology*      topo_green_[3]       = {NULL, NULL, NULL}; /**< @brief list of topos dedicated to Green's function */
    SwitchTopo*    switchtopo_green_[3] = {NULL, NULL, NULL}; /**< @brief switcher of topos for the Green's forward transform*/
    /**@} */

    // time the solve
    Profiler* prof_ = NULL;

   protected:
    /**
     * @name Data management
     * 
     * @{
     */
    void allocate_data_(const Topology *const topo[3], const Topology *topo_phys, double **data);
    void delete_switchtopos_(SwitchTopo* switchtopo[3]);
    void delete_topologies_(Topology* topo[3]);
    /**@}  */

    /**
     * @name Plan management
     * 
     * @{
     */
    void sort_plans_(FFTW_plan_dim* plan[3]);
    void init_plansAndTopos_(const Topology* topo, Topology* topomap[3], SwitchTopo* switchtopo[3], FFTW_plan_dim* planmap[3], bool isGreen);
    void allocate_plans_(const Topology* const topo[3], FFTW_plan_dim* planmap[3], double* data);
    void delete_plans_(FFTW_plan_dim* planmap[3]);
    /**@} */

    /**
     * @name SwitchTopo management
     * 
     * @{
     */
    void allocate_switchTopo_(const int ntopo, SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    void deallocate_switchTopo_(SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    void reorder_metis_(MPI_Comm comm, int *sources, int *sourcesW, int *dests, int *destsW, int *order);
    /**@} */

    /**
     * @name Do the magic
     * 
     * @{
     */
    void dothemagic_std_real(double *data);
    void dothemagic_std_complex(double *data);
    void dothemagic_rot_real_o1(double *data,const double koffset[3],const double kfact[3][3][2], const double symstart[3]);
    void dothemagic_rot_complex_o1(double *data,const double koffset[3],const double kfact[3][3][2], const double symstart[3]);
    void dothemagic_rot_real_o2(double *data,const double koffset[3],const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    void dothemagic_rot_complex_o2(double *data,const double koffset[3],const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void cmptGreenFunction_(Topology* topo[3], double* green, FFTW_plan_dim* planmap[3]);
    void cmptGreenSymmetry_(const Topology* topo, const int sym_idx, double* data, const bool isComplex);
    void scaleGreenFunction_(const Topology* topo, double* data, bool killModeZero);
    void finalizeGreenFunction_(Topology* topo_field, double* green, const Topology* topo, FFTW_plan_dim* planmap[3]);
    /**@} */

   public:
    Solver(Topology* topo, BoundaryType* rhsbc[3][2], const double h[3], const double L[3], const FLUPS_DiffType orderDiff, const FLUPS_CenterType centerType[3], Profiler* prof);
    ~Solver();

    double* setup(const bool changeTopoComm);
    const Topology* get_innerTopo_physical() ;
    const Topology* get_innerTopo_spectral() ;


    /**
     * @brief Get the total allocated size of the pointer data (returned by setup)
     * 
     * @return size_t 
     */
    size_t get_allocSize() {
        size_t size_tot = 1;
        for (int id = 0; id < ndim_; id++) {
            size_tot = std::max(topo_hat_[id]->memsize(), size_tot);
        }
        return size_tot;
    };

    /**
     * @brief Get the spectral information to compute the modes k in full spectral space
     * 
     * @param kfact  multiply the index by this factor to obtain the wave number (1/2/3 corresponds to x/y/z )
     * @param koffset  add this to the index to obtain the wave number (1/2/3 corresponds to x/y/z )
     * @param symstart  returns the first index of the symmetry
     */
    void get_spectralInfo(double kfact[3], double koffset[3], double symstart[3]) {
        for (int ip = 0; ip < 3; ip++) {
            const int dimID = plan_forward_[ip]->dimID();
            kfact[dimID]    = plan_forward_[ip]->kfact();
            symstart[dimID] = plan_forward_[ip]->symstart();
            koffset[dimID]  = plan_forward_[ip]->koffset();
        }
    }

    /**
     * @name Solver use 
     * 
     * @{
     */
    void solve(double *field, double *rhs,const FLUPS_SolverType type);
    /**@} */

    /**
     * @name Solver use (advanced)
     * 
     * @{
     */
    void do_copy(const Topology *topo, double *data, const int sign );
    void do_FFT(double *data, const int sign);
    void do_mult(double *data,const FLUPS_SolverType type);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void set_GreenType(const GreenType type) { typeGreen_ = type; }
    void set_alpha(const double alpha) { alphaGreen_ = alpha; }
    /**@} */
};

// /**
//  * @brief compute the pencil layout given the pencil direction
//  * 
//  * The pencil layout is computed so as to obtain pencils with an aspect
//  * ratio close to 1, i.e. the same number points per proc in the the 2 other directions than id.
//  * 
//  * @param id the pencil direction
//  * @param nproc the number of proc in each direction
//  * @param comm_size the total communicator size
//  * @param nglob the domain size in each direction
//  */
// static inline void pencil_nproc(const int id, int nproc[3], const int comm_size, const int nglob[3]) {
//     int id1 = (id + 1) % 3;
//     int id2 = (id + 2) % 3;

//     nproc[id] = 1;

//     double       n1       = 1;
//     double       n2       = (double) comm_size;
//     //invert indexes so that id1 is the dimension where nglob is the smallest
//     if( nglob[id1] > nglob[id2]){
//         const int tmp = id2;
//         id2 = id1;
//         id1 = tmp;
//     }
//     double       np1      = (double) nglob[id1];
//     double       np2      = (double) nglob[id2]/ comm_size;
//     const double npsquare = sqrt((double)(nglob[id1] * nglob[id2]) / comm_size);  //target number of points per dimension

//     //keep on deviding as long as ncurr/2>nsquare
//     //we want to leave n1=1, and we do not want to reach n2=1
//     while ( (np1 > npsquare) && std::floor(n2*.5) == n2*.5) {
//         n1  *= 2.0;
//         np1 *= 0.5;
//         n2  *= 0.5;
//         np2 *= 2.0;
//     }
//     nproc[id1] = (int)n1;
//     nproc[id2] = (int)n2;

//     FLUPS_INFO("my proc repartition is %d %d %d",nproc[0],nproc[1],nproc[2]);
//     if(nproc[0] * nproc[1] * nproc[2] != comm_size){
//         FLUPS_ERROR("the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size, LOCATION);
//     }
//     if(comm_size>8 && (n1==1||n2==1)){
//         FLUPS_WARNING("A slab decomposition was used instead of a pencil decomposition in direction %d. This may increase communication time.",id, LOCATION);
//         //Loss of performance may originate in slab decompositions, as an actual All2All communication is required, whereas with the pencils,
//         // we manage to do All2All communications in subcoms of size sqrt(comm_size).
//         //We could prevent this to happen by doing something like:
//         // if(n2==1){
//         //     n2*=2;
//         //     n1*=0.5;
//         // }
//     }
// }

/**
 * @brief compute the pencil layout given the pencil direction, compatible with another pencil decoposition given as a hint
 * 
 * @param id the pencil direction
 * @param nproc the number of proc in each direction
 * @param comm_size the total communicator size
 * @param id_hint the axis where we allow the proc decomposition to change
 * @param nproc_hint the number of procs in the other decomposition we want to be compatible with
 * 
 */
static inline void pencil_nproc_hint(const int id, int nproc[3], const int comm_size, const int id_hint, const int nproc_hint[3]) {
    // get the id shared between the hint topo
    int sharedID = 0;
    for (int i = 0; i < 3; i++) {
        if (i != id && i != id_hint) {
            sharedID = i;
            break;
        }
    }
    nproc[id]       = 1;
    nproc[sharedID] = nproc_hint[sharedID];
    nproc[id_hint]  = comm_size / nproc[sharedID];

    FLUPS_INFO("My proc repartition in this topo is %d %d %d",nproc[0],nproc[1],nproc[2]);
    FLUPS_CHECK(nproc[0] * nproc[1] * nproc[2] == comm_size, "the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size, LOCATION);

    if(comm_size>8 && (nproc[sharedID]==1||nproc[id_hint]==1)){
        FLUPS_WARNING("A slab decomposition was used instead of a pencil decomposition in direction %d. This may increase communication time.",id, LOCATION);
    }
}


#endif
