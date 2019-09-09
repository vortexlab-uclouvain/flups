/**
 * @file Solver.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <cstring>
#include <map>
#include "FFTW_plan_dim.hpp"
#include "defines.hpp"
#include "green_functions_2d.hpp"
#include "green_functions_3d.hpp"
#include "hdf5_io.hpp"

#include "SwitchTopo.hpp"
#include "tools.hpp"

#include "Profiler.hpp"

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
class FLUPS::Solver {
    // the memory allocation is assumed to be data[iz][iy][ix]
    // so the fastest running index is n[0] then n[1] then n[2]
    // even is the dimension is 2, we allocate arrays of dimension 3

   protected:
    int _orderdiff    = 0; /**< @brief the order of derivative (spectral = 0)  */
    int _nbr_imult    = 0; /**< @brief the number of time we have applied a DST transform */
    int _nbr_spectral = 0; /** @brief the number of spectral directions involved     */

    double  _normfact = 1.0;   /**< @brief normalization factor so that the forward/backward FFT gives output = input */
    double  _volfact  = 1.0;   /**< @brief volume factor due to the convolution computation */
    double  _hgrid[3] = {0.0}; /**< @brief grid spacing in the tranposed directions */
    double* _data     = NULL;  /**< @brief data pointer to the transposed memory */

    /**
     * @name Forward and backward 
     * transforms related objects
     */
    /**@{ */
    FFTW_plan_dim* _plan_forward[3];  /**< @brief map containing the plans for the forward fft transforms */
    FFTW_plan_dim* _plan_backward[3]; /**< @brief map containing the plans for the backward fft transforms */
    Topology*      _topo_hat[3]   = {NULL, NULL, NULL}; /**< @brief map containing the topologies (i.e. data memory layout) corresponding to each transform */
    SwitchTopo*    _switchtopo[3] = {NULL, NULL, NULL}; /**< @brief switcher of topologies for the forward transform (phys->topo[0], topo[0]->topo[1], topo[1]->topo[2]).*/
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
     * @name Do the magic
     * 
     * @{
     */
    void dothemagic_rhs_real();
    void dothemagic_rhs_complex_nmult0();
    void dothemagic_rhs_complex_nmult1();
    void dothemagic_rhs_complex_nmult2();
    void dothemagic_rhs_complex_nmult3();
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void _cmptGreenFunction(Topology* topo[3], double* green, FFTW_plan_dim* planmap[3]);
    void _cmptGreenSymmetry(const Topology* topo, const int sym_idx, double* data, const bool isComplex);
    void _scaleGreenFunction(const Topology* topo, double* data, bool killModeZero);
    /**@} */

   public:
    Solver(const Topology* topo, const BoundaryType mybc[3][2], const double h[3], const double L[3]);
    // Solver(const Topology* topo_glob,const BoundaryType mybc[3][2]);
    ~Solver();

    void setup();
    void set_OrderDiff(const int order) { _orderdiff = order; }

    /**
     * @name Solver use
     * 
     * @{
     */
    void solve(const Topology* topo, double* field, double* rhs, const SolverType type);
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
 * @param id the pencil direction
 * @param nproc the number of proc in each direction
 * @param comm_size the total communicator size
 */
static inline void _pencil_nproc(const int id, int nproc[3], const int comm_size) {
    const int id1 = (id + 1) % 3;
    const int id2 = (id + 2) % 3;

    nproc[id] = 1;

    double n1 = 1;
    double n2 = (double)comm_size;
    while (n1 < n2 && std::floor(n2) == n2) {
        n1 *= 2.0;
        n2 /= 2.0;
    }
    nproc[id1] = (int)n1;
    nproc[id2] = (int)n2;

    FLUPS_CHECK(nproc[0] * nproc[1] * nproc[2] == comm_size, "the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size);
}

#endif
