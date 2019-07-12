#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <map>
#include "defines.hpp"
#include "FFTW_plan_dim.hpp"
#include "fftw3.h"


using namespace std;

class FFTW_Solver{
    // the memory allocation is assumed to be data[iz][iy][ix]
    // so the fastest running index is n[0] then n[1] then n[2]
    int _dimorder[DIM];
    size_t _size_field[DIM];
    size_t _size_hat  [DIM];
    size_t _padstart  [DIM];

    // we do not know apriori if the data will be real or complex
    bool _isComplex;
    void* _data=NULL;

    

protected:
    multimap<int, FFTW_plan_dim*> _plan_forward;
    multimap<int, FFTW_plan_dim*> _plan_backward;

    void _init_plan_map(size_t sizeorder[DIM],int dimorder[DIM], bool* isComplex, multimap<int,FFTW_plan_dim* > *planmap) const;
    void _allocate_data(const size_t size[DIM],void** data) const;
    void _allocate_plan(const size_t size[DIM],const bool isComplex,void* data, multimap<int,FFTW_plan_dim* > *planmap) const;

public:
    FFTW_Solver(const size_t field_size[DIM],const double h[DIM],const double L[DIM],const BoundaryType mybc[DIM][2]);

};

#endif