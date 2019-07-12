#include "FFTW_Solver.hpp"


FFTW_Solver::FFTW_Solver(const size_t size_field[DIM],const double h[DIM],const double L[DIM],const BoundaryType mybc[DIM][2])
{

    //-------------------------------------------------------------------------
    // store the field size
    //-------------------------------------------------------------------------
    for(int id=0; id<DIM; id++) _size_field[id] = size_field[id];

    //-------------------------------------------------------------------------
    // for each dim, compute the plan and its type
    //-------------------------------------------------------------------------
    for(int id=0; id<DIM; id++)
    {
        // forward
        FFTW_plan_dim* myplan_forward = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_FORWARD,false);
        _plan_forward. insert(pair<int,FFTW_plan_dim*>(myplan_forward->get_type(),myplan_forward));
        // backward
        FFTW_plan_dim* myplan_backward = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_BACKWARD,false);
        _plan_backward.insert(pair<int,FFTW_plan_dim*>(myplan_backward->get_type(),myplan_backward));
    }
    
    //-------------------------------------------------------------------------
    // initialise the plans and get the sizes + dimorder
    //-------------------------------------------------------------------------
    // forward, store the size and dim order for the object
    _init_plan_map(_size_hat,_dimorder,&_isComplex,&_plan_forward);

    // backward use temporary sizes and check the output
    size_t size_tmp[DIM];
    int dimorder_tmp[DIM];
    bool isComplex_tmp;
    _init_plan_map(size_tmp,dimorder_tmp,&isComplex_tmp,&_plan_backward);
    for(int id=0; id<DIM; id++){
        assert(size_tmp[id]      == _size_hat[id]);
        assert(dimorder_tmp[id]  == _dimorder[id]);
        assert(isComplex_tmp == _isComplex);
    }

    _allocate_data(_size_hat,&_data);
    _allocate_plan(_size_hat,_isComplex,_data,&_plan_forward);
}

void FFTW_Solver::_init_plan_map(size_t sizeorder[DIM],int dimorder[DIM], bool* isComplex, multimap<int,FFTW_plan_dim* > *planmap) const
{
    //-------------------------------------------------------------------------
    // get the plan info
    //-------------------------------------------------------------------------
    (*isComplex) = false;
    // get the sizes
    size_t size_tmp[DIM];
    for(int id=0; id<DIM; id++) size_tmp[id] = _size_field[id];
    // init the plans
    int count = 0;
    for(multimap<int,FFTW_plan_dim* >::iterator it = (*planmap).begin(); it != (*planmap).end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        // initialize the plan - read only
        myplan->init(size_tmp,*isComplex);
        // update the size to the new one starting from the slowest index
        myplan->get_outsize(DIM-1-count,sizeorder);
        myplan->get_dimID  (count      ,dimorder);
        myplan->get_isComplex(isComplex);
        myplan->set_order(count);
        // display it
        myplan->disp();
        // update the counter
        count++;
    }
    
    printf("solver final order %d %d\n",dimorder[0],dimorder[1]);
    printf("solver final size %ld %ld\n",sizeorder[0],sizeorder[1]);
}

void FFTW_Solver::_allocate_data(const size_t size[DIM],void** data) const
{
    //-------------------------------------------------------------------------
    // Sanity checks
    //-------------------------------------------------------------------------
    assert((*data) == NULL);

    //-------------------------------------------------------------------------
    // Do the memory allocation
    //-------------------------------------------------------------------------
    size_t size_tot = 1;
    for(int id=0; id<DIM; id++) size_tot *= size[id];

    if(_isComplex){
        (*data) =(void*) fftw_malloc(size_tot*sizeof(fftw_complex));
    }
    else{
        (*data) =(void*) fftw_malloc(size_tot*sizeof(double));
    }

    INFOLOG("memory allocation finished\n");
}


void  FFTW_Solver::_allocate_plan(const size_t size[DIM],const bool isComplex,void* data, multimap<int,FFTW_plan_dim* > *planmap) const
{
    INFOLOG("start plan allocation\n");
    for(multimap<int,FFTW_plan_dim* >::iterator it = (*planmap).begin(); it != (*planmap).end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        // initialize the plan - read only
        myplan->allocate_plan(size,isComplex,data);
    }
}
