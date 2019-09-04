/**
 * @file FFTW_plan_dim.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#ifndef FFTW_PLAN_DIM_HPP
#define FFTW_PLAN_DIM_HPP

#include "defines.hpp"
#include "fftw3.h"

/**
 * @brief the boundary condition can be EVEN, ODD, PERiodic or UNBounded
 * 
 */
enum BoundaryType {
    EVEN = 0, /**< EVEN boundary condition = zero flux  */
    ODD  = 1, /**< ODD boundary condition = zero value */
    PER  = 3, /**< PERiodic boundary conditions */
    UNB  = 4  /**< UNBounded boundary condition */
};

/**
 * @brief PlanType is the type of plan considered and is computed as the sum of both BoundaryType variables
 * 
 * The integer value associated gives is the priority of processing.
 * We first have to do the real to real transforms, then the padded real to real (mix direction = unbounded + boundary condition),
 * then the periodic (DFT) directions and finally the padded periodic boundary condition.
 * This order is chosen in order to reduce the computational cost.
 */
enum PlanType {
    SYMSYM = 2, /**< type real 2 real (DCT / DST) : EE (0) , EO/OE (1) , OO (2) */
    MIXUNB = 5, /**< type unbounded and a symetry condition: UE/EU (4) , UO/OU (5) */
    PERPER = 6, /**< type periodic - periodic: PERPER (6) */
    UNBUNB = 8  /**< type fully unbounded UU (8) */
};

/**
 * @brief A FFTW plan in one dimension
 * 
 */
class FFTW_plan_dim {
   protected:
    const bool _isGreen; /**< @brief boolean is true if this plan is for a Green's function */
    const int  _dimID;   /**< @brief the dimension of the plan in the field reference */
    const int  _sign;    /**< @brief FFT_FORWARD (-1) or FFT_BACKWARD(+1) */

    bool   _isInputReal= true;  /**< @brief is the input of this plan real?*/
    bool   _isr2c      = false; /**< @brief is this plan the one that changes to complex?*/
    bool   _imult      = false; /**< @brief boolean to determine if we have to multiply by (i=sqrt(-1)) or not*/
    bool   _isSpectral = false; /**< @brief indicate if the Green's function has to be done spectrally (leading to a helmolz problem) */
    int    _howmany    = -1;    /**< @brief number of transfroms to perfom */
    int    _fieldstart = 0;     /**< @brief the starting index for the field copy in the direction of the plan*/
    int    _n_in       = 0;     /**< @brief the number of element in the transform*/
    int    _n_out      = 0;     /**< @brief the number of element coming out of the transform*/
    int    _symstart   = 0;     /**< @brief the starting index for the symmetry of the Green's function, set to 0 if no symmetry is needed*/
    int    _shiftgreen = 0;     /**< @brief the shift to set in the Green's function when doing the convolution*/
    double _normfact   = 0.0;   /**< @brief factor you need to multiply to get the transform on the right scaling*/
    double _volfact    = 0.0;   /**< @brief volume factor*/
    double _kfact      = 0.0;   /**< @brief multiplication factor to have the correct k numbers*/

    PlanType     _type;  /**< @brief type of this plan, see #PlanType*/
    BoundaryType _bc[2]; /**< @brief boundary condition [0]=LEFT/MIN - [1]=RIGHT/MAX*/

    fftw_r2r_kind _kind;        /**< @brief kind of transfrom to perform (used by r2r and mix plan only)*/
    fftw_plan     _plan = NULL; /**< @brief the actual FFTW plan*/

   public:
    FFTW_plan_dim(const int dimID, const double h[DIM], const double L[DIM], const BoundaryType mybc[2], const int sign, const bool isGreen);
    ~FFTW_plan_dim();

    void init(const int size[DIM], const bool isComplex, const bool isRefR2C);

    void allocate_plan(const int size_plan[DIM], double* data);
    void execute_plan();

    /**
     * @name Getters - return the value
     * 
     */
    /**@{ */
    inline bool   isSpectral() const { return _isSpectral; }
    inline bool   isInputReal() const { return _isInputReal; }
    inline bool   isr2c() const { return _isr2c; }
    inline int    dimID() const { return _dimID; }
    inline int    imult() const { return _imult; }
    inline int    shiftgreen() const { return _shiftgreen; }
    inline int    type() const { return _type; }
    inline int    symstart() const { return _symstart; }
    inline double normfact() const { return _normfact; }
    inline double volfact() const { return _volfact; }
    inline double kfact() const { return _kfact; }

    inline void get_outsize(int* size) const { size[_dimID] = _n_out; };
    inline void get_fieldstart(int* start) const { start[_dimID] = _fieldstart; };
    inline void get_isNowComplex(bool* isComplex) const { (*isComplex) = (*isComplex) || _isr2c; };
    /**@} */

    void disp();

   protected:
    /**
     * @name Initialization
     */
    /**@{ */
    void _init_real2real(const int size[DIM], bool isComplex);
    void _init_mixunbounded(const int size[DIM], bool isComplex);
    void _init_periodic(const int size[DIM], bool isComplex, bool isRefR2C);
    void _init_unbounded(const int size[DIM], bool isComplex);
    /**@} */

    /**
     * @name Plan allocation
     */
    /**@{ */
    void _allocate_plan_real(const int size_ordered[DIM], double* data);
    void _allocate_plan_complex(const int size_ordered[DIM], double* data);
    /**@} */
};

#endif