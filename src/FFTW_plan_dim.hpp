#ifndef FFTW_PLAN_DIM_HPP
#define FFTW_PLAN_DIM_HPP

#include "defines.hpp"
#include "fftw3.h"

enum BoundaryType
{
    //The sum of the BoundaryType gives the priority. The lowest sum is to perfom first
    // DCT/DST < mix < unbounded < periodic
    EVEN    = 0, // even boundary condition = zero flux
    ODD     = 1, // odd boundary condition = zero value
    PER     = 3, // periodic boundary conditions
    UNB     = 4  // unbounded boundary condition
};

enum SolverType
{
    R2R     = 2, // real 2 real (DCT / DST) : EE (0) , EO/OE (1) , OO (2)
    MIX     = 5, // unbounded and a symetry condition: UE/EU (4) , UO/OU (5)
    PERPER  = 6, // periodic - periodic: PERPER (6)
    UNBUNB  = 8  // fully unbounded UU (8)
};


class FFTW_plan_dim{

protected:
    const int _dimID; // my ID in the dimension - used to fill the dimension arrays
    const int _sign; // FFT_FORWARD (-1) or FFT_BACKWARD(+1)
    const bool _isGreen ;

    bool   _imult; // boolean to determine if we have to multiply by (i=sqrt(-1)) or not
    double _fact;  // factor you need to multiply to get the transform on the right scaling
    double _k_fact;
    bool _isComplex; // is this plan complex? 
    bool _isDataComplex; // is the data allocated complex? (yes if any plan is complex)
    int _orderID; // my ID in the !!ordered!! dimension

    SolverType _type;
    BoundaryType _bc[2]; // boundary condition [0]=LEFT/MIN - [1]=RIGHT/MAX
    size_t _padstart;
    
    size_t _n_in; // the number of element in the transform
    size_t _n_out; // the number of element coming out of the transform
    
    fftw_r2r_kind _kind;
    fftw_plan* _plan = NULL;
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
    void allocate_plan(const size_t size_plan[DIM],const bool isComplex, void* data);

    int  get_type()      const;
    bool get_isComplex() const;
    void get_outsize  (size_t* size ) const;
    void get_padstart (size_t* start) const;
    void get_isComplex(bool* isComplex) const;
    void get_dimID   (const int id, int    dimID[DIM]) const;
    void get_outsize (const int id, size_t size [DIM]) const;
    void get_padstart(const int id, size_t start[DIM]) const;

    void set_order(const int id);
    
    void disp();

protected:
    void _init_real2real (const size_t size[DIM],bool isComplex);
    void _init_mixpoisson(const size_t size[DIM],bool isComplex);
    void _init_periodic  (const size_t size[DIM],bool isComplex);
    void _init_unbounded (const size_t size[DIM],bool isComplex);

    void _allocate_plan_real2real (const size_t size_plan[DIM],void* data);
    // void _allocate_plan_mixpoisson(const size_t size_plan[DIM],double* data);
    // void _allocate_plan_periodic  (const size_t size_plan[DIM],double* data);
    // void _allocate_plan_unbounded (const size_t size_plan[DIM],double* data);

};

#endif