/**
 * @file Solver.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
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

#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <cstring>
#include <map>
#include "FFTW_plan_dim.hpp"
#include "defines.hpp"
#include "green_functions_3d.hpp"
#include "hdf5_io.hpp"

#include "SwitchTopo.hpp"
#include "SwitchTopo_a2a.hpp"
#include "SwitchTopo_nb.hpp"

#include "Profiler.hpp"
#include "omp.h"

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
 * const size_t id_transposed = iz*_dim_multfact[2] + iy*_dim_multfact[1] + ix*_dim_multfact[0] + _offset;
 * \endcode
 *  
 * 
 */
class Solver {
    // the memory allocation is assumed to be data[iz][iy][ix]
    // so the fastest running index is n[0] then n[1] then n[2]
    // even is the dimension is 2, we allocate arrays of dimension 3

   protected:
    int _fftwalignment = 0; /**< @brief alignement assumed by the FFTW Solver  */
    int _orderdiff     = 0; /**< @brief the order of derivative (spectral = 0)  */
    int _nbr_imult     = 0; /**< @brief the number of time we have applied a DST transform */
    int _nbr_spectral  = 0; /** @brief the number of spectral directions involved     */

    double  _normfact = 1.0;   /**< @brief normalization factor so that the forward/backward FFT gives output = input */
    double  _volfact  = 1.0;   /**< @brief volume factor due to the convolution computation */
    double  _hgrid[3] = {0.0}; /**< @brief grid spacing in the tranposed directions */
    double* _data     = NULL;  /**< @brief data pointer to the transposed memory */

    /**
     * @name Forward and backward 
     * transforms related objects
     */
    /**@{ */
    FFTW_plan_dim*  _plan_forward[3];  /**< @brief map containing the plans for the forward fft transforms */
    FFTW_plan_dim*  _plan_backward[3]; /**< @brief map containing the plans for the backward fft transforms */
    const Topology* _topo_phys     = NULL;
    Topology*       _topo_hat[3]   = {NULL, NULL, NULL}; /**< @brief map containing the topologies (i.e. data memory layout) corresponding to each transform */
    SwitchTopo*     _switchtopo[3] = {NULL, NULL, NULL}; /**< @brief switcher of topologies for the forward transform (phys->topo[0], topo[0]->topo[1], topo[1]->topo[2]).*/
    
    opt_double_ptr _sendBuf = NULL; /**<@brief The send buffer for _switchtopo */
    opt_double_ptr _recvBuf = NULL; /**<@brief The recv buffer for _switchtopo */
    /**@} */

    /**
     * @name Green's function (and corresponding forward transform) related vars and objects
     * 
     */
    /**@{ */
    // int       _shiftgreen[3]   = {0, 0, 0}; /**< @brief the shift in the Green's function which chose to take the flip-flop mode or not */
    double    _alphaGreen      = 2.0;       /**< @brief regularization parameter for HEJ_* Green's functions */
    double*   _green           = NULL;      /**< @brief data pointer to the transposed memory for Green */
    GreenType _typeGreen       = CHAT_2;    /**< @brief the type of Green's function */

    FFTW_plan_dim* _plan_green[3]; /**< @brief map containing the plan for the Green's function */
    Topology*   _topo_green[3]       = {NULL, NULL, NULL}; /**< @brief list of topos dedicated to Green's function */
    SwitchTopo* _switchtopo_green[3] = {NULL, NULL, NULL}; /**< @brief switcher of topos for the Green's forward transform*/
    /**@} */

    // time the solve
    Profiler* _prof = NULL;

   protected:
    /**
     * @name Data management
     * 
     * @{
     */
    void _allocate_data(const Topology* const topo[3], double** data);
    void _delete_switchtopos(SwitchTopo* switchtopo[3]);
    void _delete_topologies(Topology* topo[3]);
    /**@}  */

    /**
     * @name Plan management
     * 
     * @{
     */
    void _sort_plans(FFTW_plan_dim* plan[3]);
    void _init_plansAndTopos(const Topology* topo, Topology* topomap[3], SwitchTopo* switchtopo[3], FFTW_plan_dim* planmap[3], bool isGreen);
    void _allocate_plans(const Topology* const topo[3], FFTW_plan_dim* planmap[3], double* data);
    void _delete_plans(FFTW_plan_dim* planmap[3]);
    /**@} */

    /**
     * @name SwitchTopo management
     * 
     * @{
     */
    void _allocate_switchTopo(const int ntopo, SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    void _deallocate_switchTopo(SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    /**@} */

    /**
     * @name Do the magic
     * 
     * @{
     */
    void dothemagic_rhs_real(double *data);
    void dothemagic_rhs_complex_nmult0(double *data);
    void dothemagic_rhs_complex_nmult1(double *data);
    void dothemagic_rhs_complex_nmult2(double *data);
    void dothemagic_rhs_complex_nmult3(double *data);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void _cmptGreenFunction(Topology* topo[3], double* green, FFTW_plan_dim* planmap[3]);
    void _cmptGreenSymmetry(const Topology* topo, const int sym_idx, double* data, const bool isComplex);
    void _scaleGreenFunction(const Topology* topo, double* data, bool killModeZero);
    void _finalizeGreenFunction(Topology* topo_field[3],double* green, Topology* topo[3], FFTW_plan_dim* plans[3]);
    /**@} */

   public:
    Solver(const Topology* topo, const BoundaryType mybc[3][2], const double h[3], const double L[3],Profiler* prof);
    // Solver(const Topology* topo_glob,const BoundaryType mybc[3][2]);
    ~Solver();

    void setup();
    void set_OrderDiff(const int order) { _orderdiff = order; }
    Topology* get_innerTopo_physical() ;
    Topology* get_innerTopo_spectral() ;

    size_t get_allocSize() {
        size_t size_tot = 1;
        for (int id = 0; id < 3; id++) {
            size_tot = std::max(_topo_hat[id]->memsize(), size_tot);
        }
        return size_tot;
    };

    /**
     * @name Solver use 
     * 
     * @{
     */
    void solve(double* field, double* rhs, const SolverType type);
    /**@} */

    /**
     * @name Solver use (advanced)
     * 
     * @{
     */
    void do_copy(const Topology *topo, double *data, const int sign );
    void do_FFT(double *data, const int sign);
    void do_mult(double *data, const SolverType type);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void set_GreenType(const GreenType type) { _typeGreen = type; }
    void set_alpha(const double alpha) { _alphaGreen = alpha; }
    /**@} */
};

/**
 * @brief compute the pencil layout given the pencil direction
 * 
 * The pencil layout is computed so as to obtain pencils with an aspect
 * ratio close to 1, i.e. the same number points per proc in the the 2 other directions than id.
 * 
 * @param id the pencil direction
 * @param nproc the number of proc in each direction
 * @param comm_size the total communicator size
 * @param nglob the domain size in each direction
 */
static inline void pencil_nproc(const int id, int nproc[3], const int comm_size, const int nglob[3]) {
    int id1 = (id + 1) % 3;
    int id2 = (id + 2) % 3;

    nproc[id] = 1;

    double       n1       = 1;
    double       n2       = (double) comm_size;
    //invert indexes so that id1 is the dimension where nglob is the smallest
    if( nglob[id1] > nglob[id2]){
        const int tmp = id2;
        id2 = id1;
        id1 = tmp;
    }
    double       np1      = (double) nglob[id1];
    double       np2      = (double) nglob[id2]/ comm_size;
    const double npsquare = sqrt((double)(nglob[id1] * nglob[id2]) / comm_size);  //target number of points per dimension

    //keep on deviding as long as ncurr/2>nsquare
    //we want to leave n1=1, and we do not want to reach n2=1
    while ( (np1 > npsquare) && std::floor(n2*.5) == n2*.5) {
        n1  *= 2.0;
        np1 *= 0.5;
        n2  *= 0.5;
        np2 *= 2.0;
    }
    nproc[id1] = (int)n1;
    nproc[id2] = (int)n2;

    FLUPS_INFO("my proc repartition is %d %d %d\n",nproc[0],nproc[1],nproc[2]);
    if(nproc[0] * nproc[1] * nproc[2] != comm_size){
        FLUPS_ERROR("the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size, LOCATION);
    }
    if(comm_size>8 && (n1==1||n2==2)){
        FLUPS_WARNING("A slab decomposition was used instead of a pencil decomposition in direction %d. This may increase communication time.",id, LOCATION);
        //Loss of performance may originate in slab decompositions, as an actual All2All communication is required, whereas with the pencils,
        // we manage to do All2All communications in subcoms of size sqrt(comm_size).
        //We could prevent this to happen by doing something like:
        // if(n2==1){
        //     n2*=2;
        //     n1*=0.5;
        // }
    }
}

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

    FLUPS_INFO("my proc repartition is %d %d %d\n",nproc[0],nproc[1],nproc[2]);
    FLUPS_CHECK(nproc[0] * nproc[1] * nproc[2] == comm_size, "the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size, LOCATION);
}

#endif
