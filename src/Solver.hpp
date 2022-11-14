/**
 * @file Solver.hpp
 * @copyright Copyright © UCLouvain 2020
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <cstring>
#include <map>
#include <array>
#include <tuple>

#include "FFTW_plan_dim.hpp"
#include "defines.hpp"
#include "green_functions.hpp"
#include "hdf5_io.hpp"

#if (FLUPS_MPI_AGGRESSIVE)
#include "SwitchTopoX_a2a.hpp"
#include "SwitchTopoX_isr.hpp"
#include "SwitchTopoX_nb.hpp"
#else
#include "SwitchTopo.hpp"
#include "SwitchTopo_a2a.hpp"
#include "SwitchTopo_nb.hpp"
#endif

#include "omp.h"

#ifdef HAVE_METIS
#include "metis.h"
#endif

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
    bool     skip_st0_      = false;  //!< dictates the solver to skip the first Switchtopo and that the "physical info" are given according to topo_hat_[0]
    int      lda_           = 1;      //!< the number of components of the problem, i.e. 2D or 3D */
    int      ndim_          = 3;      //!< the dimension of the problem, i.e. 2D or 3D */
    int      fftwalignment_ = 0;      //!< alignement assumed by the FFTW Solver  */
    DiffType odiff_         = NOD;    //!< the order of derivative (spectral = SPE, 2nd order FD = FD2) */
    double   normfact_      = 1.0;    //!< normalization factor so that the forward/backward FFT gives output = input */
    double   volfact_       = 1.0;    //!< volume factor due to the convolution computation */
    double   hgrid_[3]      = {0.0};  //!< grid spacing in the tranposed directions */
    double*  data_          = NULL;   //!< data pointer to the transposed memory */

    /**
     * @name Forward and backward
     * transforms related objects
     */
    /**@{ */
    FFTW_plan_dim* plan_forward_[3];       /**< @brief map containing the plans for the forward fft transforms */
    FFTW_plan_dim* plan_backward_[3];      /**< @brief map containing the plans for the backward fft transforms */
    FFTW_plan_dim* plan_backward_diff_[3]; /**< @brief map containing the plans for the backward fft transforms */

    Topology* topo_phys_   = NULL;
    Topology* topo_hat_[3] = {NULL, NULL, NULL}; /**< @brief map containing the topologies (i.e. data memory layout) corresponding to each transform */

#if (FLUPS_MPI_AGGRESSIVE)
    SwitchTopoX* switchtopo_[3] = {NULL, NULL, NULL}; /**< @brief switcher of topologies for the forward transform (phys->topo[0], topo[0]->topo[1], topo[1]->topo[2]).*/
#else
    SwitchTopo*    switchtopo_[3]       = {NULL, NULL, NULL}; /**< @brief switcher of topologies for the forward transform (phys->topo[0], topo[0]->topo[1], topo[1]->topo[2]).*/
#endif

#if (FLUPS_MPI_AGGRESSIVE)
    m_ptr_t sendBuf_;
    m_ptr_t recvBuf_;
#else
    opt_double_ptr sendBuf_             = NULL;               /**<@brief The send buffer for switchtopo_ */
    opt_double_ptr recvBuf_             = NULL;               /**<@brief The recv buffer for switchtopo_ */
#endif
    /**@} */

    /**
     * @name Green's function (and corresponding forward transform) related vars and objects
     *
     */
    /**@{ */
    double    alphaGreen_ = 2.0;    /**< @brief regularization parameter for HEJ_* Green's functions */
    double*   green_      = NULL;   /**< @brief data pointer to the transposed memory for Green */
    GreenType typeGreen_  = CHAT_2; /**< @brief the type of Green's function */

    FFTW_plan_dim* plan_green_[3];                      /**< @brief map containing the plan for the Green's function */
    Topology*      topo_green_[3] = {NULL, NULL, NULL}; /**< @brief list of topos dedicated to Green's function */

#if (FLUPS_MPI_AGGRESSIVE)
    SwitchTopoX* switchtopo_green_[3] = {NULL, NULL, NULL}; /**< @brief switcher of topos for the Green's forward transform*/
#else
    SwitchTopo*    switchtopo_green_[3] = {NULL, NULL, NULL}; /**< @brief switcher of topos for the Green's forward transform*/
#endif
    /**@} */

    // time the solve
    H3LPR::Profiler* prof_ = NULL;

   protected:
    /**
     * @name Data management
     *
     * @{
     */
    void allocate_data_(const Topology* const topo[3], const Topology* topo_phys, double** data);
#if (FLUPS_MPI_AGGRESSIVE)
    void delete_switchtopos_(SwitchTopoX* switchtopo[3]);
#else
    void           delete_switchtopos_(SwitchTopo* switchtopo[3]);
#endif
    void delete_topologies_(Topology* topo[3]);
    /**@}  */

    /**
     * @name Plan management
     *
     * @{
     */
#if (FLUPS_MPI_AGGRESSIVE)
    void init_plansAndTopos_(const Topology* topo, Topology* topomap[3], SwitchTopoX* switchtopo[3], FFTW_plan_dim* planmap[3], bool isGreen);
#else
    void           init_plansAndTopos_(const Topology* topo, Topology* topomap[3], SwitchTopo* switchtopo[3], FFTW_plan_dim* planmap[3], bool isGreen);
#endif
    void allocate_plans_(const Topology* const topo[3], FFTW_plan_dim* planmap[3], double* data);
    void delete_plans_(FFTW_plan_dim* planmap[3]);
    /**@} */

    /**
     * @name SwitchTopo management
     *
     * @{
     */
#if (FLUPS_MPI_AGGRESSIVE)
    void allocate_switchTopo_(const int ntopo, SwitchTopoX** switchtopo, m_ptr_t* send_buff, m_ptr_t* recv_buff);
    void deallocate_switchTopo_(SwitchTopoX** switchtopo, m_ptr_t* send_buff, m_ptr_t* recv_buff);
#else
    void           allocate_switchTopo_(const int ntopo, SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    void           deallocate_switchTopo_(SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
#endif
    void reorder_metis_(MPI_Comm comm, int* sources, int* sourcesW, int* dests, int* destsW, int* order);
    /**@} */

    /**
     * @name Do the magic
     *
     * @{
     */
    void dothemagic_std_real(double* data);
    void dothemagic_std_complex(double* data);
    void dothemagic_rot_real_o1(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3]);
    void dothemagic_rot_complex_o1(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3]);
    void dothemagic_rot_real_o2(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    void dothemagic_rot_complex_o2(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    void dothemagic_rot_real_o4(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    void dothemagic_rot_complex_o4(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    void dothemagic_rot_real_o6(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
    void dothemagic_rot_complex_o6(double* data, const double koffset[3], const double kfact[3][3][2], const double symstart[3], const double hgrid[3]);
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
    Solver(Topology* topo, BoundaryType* rhsbc[3][2], const double h[3], const double L[3], const DiffType orderDiff, const CenterType centerType[3], H3LPR::Profiler* prof);
    ~Solver();

    void      setup(const bool changeTopoComm);
    Topology* get_innerTopo_physical();
    Topology* get_innerTopo_spectral();

    double* get_innerBuffer() { return data_; };

    void skip_firstSwitchtopo() { skip_st0_ = true; };

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
    void solve(double* field, double* rhs, const SolverType type);
    /**@} */

    /**
     * @name Solver use (advanced)
     *
     * @{
     */
    void do_copy(const Topology* topo, double* data, const int sign);
    void do_FFT(double* data, const int sign);
    void do_mult(double* data, const SolverType type);
    /**@} */

    /**
     * @name Green's function
     *
     * @{
     */
    void set_GreenType(const GreenType type) { typeGreen_ = type; }
    void set_alpha(const double alpha) { alphaGreen_ = alpha; }
    /**@} */

    /**
     * @name Print MPI info of the Switchtopos
     *
     */
    void get_switchtopo_info() {
        for (int i = 0; i < lda_; i++) {
            switchtopo_[i]->print_info();
        }
    };
};

/**
 * @brief compute the pencil layout given the pencil direction, compatible with another pencil decoposition given as a hint
 *
 * @param id the desired pencil direction
 * @param nproc the targeted number of proc in each direction (output)
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

    FLUPS_INFO("My proc repartition in this topo is %d %d %d", nproc[0], nproc[1], nproc[2]);
    FLUPS_CHECK(nproc[0] * nproc[1] * nproc[2] == comm_size, "the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size);

    if (comm_size > 8 && (nproc[sharedID] == 1 || nproc[id_hint] == 1)) {
        FLUPS_WARNING("A slab decomposition was used instead of a pencil decomposition in direction %d. This may increase communication time.", id);
    }
}

#endif
