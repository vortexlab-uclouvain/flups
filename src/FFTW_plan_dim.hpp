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
enum BoundaryType
{
    EVEN    = 0, /**< EVEN boundary condition = zero flux  */
    ODD     = 1, /**< ODD boundary condition = zero value */
    PER     = 3, /**< PERiodic boundary conditions */
    UNB     = 4  /**< UNBounded boundary condition */
};

/**
 * @brief SolverType is the type of solver considered and is computed as the sum of both BoundaryType variables
 * 
 * The integer value associated gives is the priority of processing.
 * We first have to do the real to real transforms, then the padded real to real (mix direction = unbounded + boundary condition),
 * then the periodic (DFT) directions and finally the padded periodic boundary condition.
 * This order is chosen in order to reduce the computational cost.
 * 
 */
enum SolverType
{
    R2R     = 2, /**< type real 2 real (DCT / DST) : EE (0) , EO/OE (1) , OO (2) */ 
    MIX     = 5, /**< type unbounded and a symetry condition: UE/EU (4) , UO/OU (5) */
    PERPER  = 6, /**< type periodic - periodic: PERPER (6) */
    UNBUNB  = 8  /**< type fully unbounded UU (8) */
};


/**
 * @brief A FFTW plan in one dimension
 * 
 */
class FFTW_plan_dim{

protected:
    const int _dimID; /**< the dimension of the plan in the field reference */
    const int _sign; /**< FFT_FORWARD (-1) or FFT_BACKWARD(+1) */
    const bool _isGreen ;


    bool   _imult; // boolean to determine if we have to multiply by (i=sqrt(-1)) or not
    double _normfact;  // factor you need to multiply to get the transform on the right scaling
    double _volfact; 
    double _k_fact;
    bool _isDataComplex; // is the data allocated complex? (yes if any plan is complex)
    bool _switch2Complex; // is this plan the one that changes to complex?
    int _orderID; // my ID in the !!ordered!! dimension

    SolverType _type;
    BoundaryType _bc[2]; // boundary condition [0]=LEFT/MIN - [1]=RIGHT/MAX
    size_t _fieldstart;
    
    size_t _n_in; // the number of element in the transform
    size_t _n_out; // the number of element coming out of the transform
    
    fftw_r2r_kind _kind;
    fftw_plan _plan = NULL;
    int _howmany;


    // // data's for Green to be destroyed one performed
    // size_t _n_green;
    // fftw_plan* _plan_green;
    // fftw_r2r_kind _kind_green;
    // int _howmany;
    // int _istride;
    // int _ostride;
    // int _idist;
    // int _odist;
    
public:

    FFTW_plan_dim(const int dimID,const double h[DIM],const double L[DIM],const BoundaryType mybc[2],const int sign, const bool isGreen);
    ~FFTW_plan_dim();

    void init(const size_t size[DIM],const bool isComplex);
    void allocate_plan(const size_t size_plan[DIM],const size_t offset,const bool isComplex, void* data);

    void execute_plan();

    int  get_type()      const;
    bool get_isComplex() const;
    double get_normfact() const;
    double get_volfact() const;
    void get_outsize  (size_t* size ) const;
    void get_fieldstart (size_t* start) const;
    void get_isComplex(bool* isComplex) const;
    void get_dimID   (const int id, int    dimID[DIM]) const;
    void get_outsize (const int id, size_t size [DIM]) const;
    void get_fieldstart(const int id, size_t start[DIM]) const;
    
    void set_order(const int id);

    void disp();

protected:
    void _init_real2real (const size_t size[DIM],bool isComplex);
    void _init_mixpoisson(const size_t size[DIM],bool isComplex);
    void _init_periodic  (const size_t size[DIM],bool isComplex);
    void _init_unbounded (const size_t size[DIM],bool isComplex);

    void _allocate_plan_real   (const size_t size_ordered[DIM],const size_t offset, void* data);
    void _allocate_plan_complex(const size_t size_ordered[DIM],void* data);

};

#endif