/**
 * @file FFTW_Solver.hpp
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

#include <map>
#include "defines.hpp"
#include "FFTW_plan_dim.hpp"
#include "green_functions.hpp"
#include "fftw3.h"

#include "tools.hpp"


using namespace std;


/**
 * @brief Type of Poisson equation solved
 * 
 */
enum SolverType{
    UP_SRHS, /**<@brief scalar \f$ \nabla^2 f = rhs \f$ */
    UP_VRHS, /**<@brief vectorial \f$ \nabla^2 f = rhs \f$ */
    UP_ROT, /**<@brief vectorial \f$ \nabla^2 f = \nabla \times rhs \f$ */
    UP_DIV /**<@brief scalar \f$ \nabla^2 f = \nabla \cdot rhs \f$ */
};

/**
 * @brief The Poisson solver
 * 
 * 
 * A collection of 2 or 3 FFTW_plan_dim and a Green's function to solve the Poisson equation.
 * The tranformation are done in-place in each direction successively
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
class FFTW_Solver{
    // the memory allocation is assumed to be data[iz][iy][ix]
    // so the fastest running index is n[0] then n[1] then n[2]
    // even is the dimension is 2, we allocate arrays of dimension 3

protected:
    int    _type; /**< @brief the type of the solver, see #SolverType */
    int    _nbr_imult       = 0;        /**< @brief the number of time we have applied a DST transform */
    int    _dimorder    [3] = {0,1,2};  /**< @brief the transposed order of the dimension as used throughout the transforms*/
    int _size_field  [3] = {1,1,1};  /**< @brief the size of the field that comes in in double indexing unit  */
    int _fieldstart  [3] = {0,0,0};  /**< @brief the place in memory (in double indexing unit) where we start copying the rhs and the source terms in the extended domain */
    int _size_hat    [3] = {1,1,1};  /**< @brief the size of the transform in fftw_complex indexing unit if the #_isComplex is true, in double indexing unit otherwize */
    int _dim_multfact[3] = {1,1,1};  /**< @brief the multiplication factors used to transpose the data. */
    int _shiftgreen[3] = {0,0,0}; 
    size_t _offset          = 0;        /**< @brief the offset in memory in double indexing unit due to #_fieldstart (only used for a R2R transfrom!)*/
    double _hgrid       [3] = {0.0}; /**< @brief grid spacing in the tranposed directions */

    bool _isComplex = false;    /**< @brief boolean to indicate if the transfrom data is complex (true) or not */
    double* _data   = NULL;     /**< @brief data pointer to the transposed memory */

    double _normfact = 1.0; /**< @brief normalization factor so that the forward/backward FFT gives output = input */
    double _volfact  = 1.0; /**< @brief volume factor due to the convolution computation */

    multimap<int, FFTW_plan_dim*> _plan_forward;    /**< @brief map containing the plan forward  */
    multimap<int, FFTW_plan_dim*> _plan_backward;   /**< @brief map containing the plan backward */
    
    /**
     * @name Green's function 
     * 
     */
    /**@{ */
    int _size_hat_green [3] = {1,1,1};  /**< @brief the size of the Green's transformed in fftw_complex indexing unit if the #_isComplex is true, in double indexing unit otherwize */
    double* _green           = NULL; /**< @brief data pointer to the transposed memory for Green */
    OrderDiff _greenorder = CHAT_2; /**< @brief order and type of the Green function, see #OrderGreen */
    OrderDiff _greendiff = DIF_2; /**< @brief order of the spectral differentiation, see #OrderGreen */
    double _greenalpha = 2.0; /**< @brief regularization parameter for HEJ_* Green's functions */
    multimap<int, FFTW_plan_dim*> _plan_green;      /**< @brief map containing the plan for the Green's function */
    /**@} */

    
protected:    
    
    /**
     * @name Data management
     * 
     * @{
     */
    void _allocate_data(const int size[3],double** data);
    void _deallocate_data(double* data);
    /**@}  */

    /**
     * @name Plan management
     * 
     * @{
     */
    void _init_plan_map(int sizeorder[3], int fieldstart[3], int dimorder[3], bool* isComplex, multimap<int,FFTW_plan_dim* > *planmap);
    void _allocate_plan(const int size[3],const size_t offset, const bool isComplex,double* data, multimap<int,FFTW_plan_dim* > *planmap);
    void _delete_plan(multimap<int,FFTW_plan_dim* > *planmap);
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
    void _compute_Green(const int size_green[3],double* green, multimap<int,FFTW_plan_dim* >* planmap);
    /**@} */

public:
    FFTW_Solver(const int field_size[DIM],const double h[DIM],const double L[DIM],const BoundaryType mybc[DIM][2]);
    ~FFTW_Solver();

    void setup(const SolverType mytype);

    /**
     * @name Solver use
     * 
     * @{
     */
    void solve(double* field, double* rhs);
    // void solve_rotrhs(double* field, double* rhs);
    // void solve_div_rhs(double* field, double* rhs);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void set_GreenType(const OrderDiff order);
    void set_GreenDiff(const OrderDiff order);
    void set_alpha(const double alpha);
    /**@} */
};

#endif