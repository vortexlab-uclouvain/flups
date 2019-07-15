#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <map>
#include "defines.hpp"
#include "FFTW_plan_dim.hpp"
#include "fftw3.h"

#include "tools.hpp"


using namespace std;

class FFTW_Solver{
    // the memory allocation is assumed to be data[iz][iy][ix]
    // so the fastest running index is n[0] then n[1] then n[2]
    // even is the dimension is 2, we allocate arrays of dimension 3
    int    _dimorder   [3] = {0,1,2};
    size_t _size_field [3] = {1,1,1};
    size_t _fieldstart [3] = {0,0,0};
    size_t _size_hat   [3] = {1,1,1};

    size_t _dim_multfact[3] = {1,1,1};
    
    size_t _offset = 0;

    // we do not know apriori if the data will be real or complex
    bool _isComplex;
    void* _data=NULL;

    

protected:
    multimap<int, FFTW_plan_dim*> _plan_forward;
    multimap<int, FFTW_plan_dim*> _plan_backward;

    void _init_plan_map(size_t sizeorder[3], size_t fieldstart[3], int dimorder[3], bool* isComplex, multimap<int,FFTW_plan_dim* > *planmap);
    
    void _allocate_data(const size_t size[DIM],void** data);
    void _deallocate_data(void* data);

    void _allocate_plan(const size_t size[DIM],const bool isComplex,void* data, multimap<int,FFTW_plan_dim* > *planmap) const;
    void _delete_plan(multimap<int,FFTW_plan_dim* > *planmap);

public:
    FFTW_Solver(const size_t field_size[DIM],const double h[DIM],const double L[DIM],const BoundaryType mybc[DIM][2]);
    ~FFTW_Solver();

    void solve(double* field, double* rhs);

};

#endif