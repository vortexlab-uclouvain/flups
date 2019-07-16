/**
 * @file FFTW_Solver.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "FFTW_Solver.hpp"

/*!
 * @brief Construct a new fftw solver
 * 
 * @param[in] size_field the size of the field (real) that has to be solved
 * @param[in] h the grid spacing uniform in each direction
 * @param[in] L the length of the computational domain
 * @param[in] mybc the boundary conditions asked
 */
FFTW_Solver::FFTW_Solver(const size_t size_field[DIM],const double h[DIM],const double L[DIM],const BoundaryType mybc[DIM][2])
{
    BEGIN_FUNC
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
    _init_plan_map(_size_hat,_fieldstart,_dimorder,&_isComplex,&_plan_forward);

    // backward uses temporary sizes and check the output
    size_t size_tmp      [3] = {1,1,1};
    size_t fieldstart_tmp[3] = {0,0,0};
    int    dimorder_tmp  [3] = {0,1,2};
    bool   isComplex_tmp     = false;
    _init_plan_map(size_tmp,fieldstart_tmp,dimorder_tmp,&isComplex_tmp,&_plan_backward);
    // sanity checks
    for(int id=0; id<3; id++){
        assert(size_tmp[id]       == _size_hat[id]);
        assert(fieldstart_tmp[id] == _fieldstart[id]);
        assert(dimorder_tmp[id]   == _dimorder[id]);
        assert(isComplex_tmp      == _isComplex);
    }

    //-------------------------------------------------------------------------
    // Store some usefull factors in 'double' index calculus
    //-------------------------------------------------------------------------
    // _dim_multfact is used to compute the tranposed index location = 
    // _dim_multfact[0] * ix +  _dim_multfact[1] * iy + _dim_multfact[2] * iz
    size_t acc = 1;
    for(int id=0; id<3; id++){
        _dim_multfact[_dimorder[id]] = acc;
        acc *= _size_hat[id];
    }
    // if we are complex, we have to double the indexes to use them with doubles, except the first ones
    if(_isComplex) for(int id=1; id<3; id++) _dim_multfact[_dimorder[id]] *=2;

    if(_isComplex){
        _offset = _fieldstart[2]*_size_hat[1]*(_size_hat[0]*2)  + _fieldstart[1] * (_size_hat[0]*2) + _fieldstart[0];
    }
    else{
        _offset = _fieldstart[2]*_size_hat[1]*_size_hat[0]      + _fieldstart[1] * _size_hat[0]     + _fieldstart[0];
    }

    //-------------------------------------------------------------------------
    // Get the normalization factor
    //-------------------------------------------------------------------------
    _normfact = 1.0;
    _volfact  = 1.0;
    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_forward.begin(); it != _plan_forward.end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        _normfact *= myplan->get_normfact();
        _volfact  *= myplan->get_volfact ();
    }

    //-------------------------------------------------------------------------
    // allocate the data
    //-------------------------------------------------------------------------
    _allocate_data(_size_hat,&_data);

    //-------------------------------------------------------------------------
    // allocate the plans
    //-------------------------------------------------------------------------
    _allocate_plan(_size_hat,_offset,_isComplex,_data,&_plan_forward );
    _allocate_plan(_size_hat,_offset,_isComplex,_data,&_plan_backward);
}

/**
 * @brief Destroy the fftw solver
 * 
 */
FFTW_Solver::~FFTW_Solver(){
    BEGIN_FUNC
    _delete_plan(&_plan_forward);
    _delete_plan(&_plan_backward);

    _deallocate_data(_data);
}
/**
 * @brief delete the FFTW_plan_dim stored in planmap
 * 
 * @param planmap 
 */
void FFTW_Solver::_delete_plan(multimap<int,FFTW_plan_dim* > *planmap){
    BEGIN_FUNC
    // deallocate the plans
    for(multimap<int,FFTW_plan_dim* >::iterator it = (*planmap).begin(); it != (*planmap).end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        delete myplan;
    }
}

/**
 * @brief initialize a map containing plan
 * 
 * @param[out] sizeorder the size of the tranformed field in the correct order
 * @param[out] dimorder the order of 
 * @param[out] isComplex indicate if the data needed by the plan is complex or not
 * @param[in/out] planmap the 2 (or 3) plans for each direction to initialize
 */
void FFTW_Solver::_init_plan_map(size_t sizeorder[3], size_t fieldstart[3], int dimorder[3], bool* isComplex, multimap<int,FFTW_plan_dim* > *planmap)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // get the plan info
    //-------------------------------------------------------------------------
    (*isComplex) = false;
    // get the sizes
    size_t size_tmp[DIM];
    for(int id=0; id<DIM; id++) size_tmp[id] = _size_field[id];
    // init the plans
    int count = 0; // index from 0 to DIM-1 following the priority given by the type
    int orderID = 0;  // the order ID to give to the plan in the memory
    int max_count = DIM-1; // max bound to compute the order in memory

    // by default, we put the first plans on the slowest rotating index
    // if we have a R2C, it HAS to be on the faster rotating index
    for(multimap<int,FFTW_plan_dim* >::iterator it = (*planmap).begin(); it != (*planmap).end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        // initialize the plan - read only
        myplan->init(size_tmp,*isComplex);
        myplan->get_outsize(size_tmp); // update the size for the next plans, keep order unchanged
        
        // get the orderID of the plan
        if(!(*isComplex)){ // if we are not complex yet
            myplan->get_isComplex(isComplex);
            // if we just changed to complex = R2C, we have to put the plan on the fastest index = 0
            // we also increment the max_count to prevent any other plan reaching this 0 index
            if((*isComplex)){
                max_count += 1;
                orderID    = 0;
            }
            else{
                // if we didn't became complex, put the latest plan in the fastest indexing direction
                orderID = max_count-count;
            }
        }
        else{ // already complex, put the latest plan in the fastest rotating index
            orderID = max_count-count;
        }
        myplan->get_outsize    (orderID,sizeorder); // store the size in the correct order
        myplan->get_dimID      (orderID,dimorder); // store the correspondance of the transposition
        myplan->get_fieldstart (orderID,fieldstart);
        // set the orderID
        myplan->set_order(orderID);
        // display it
        myplan->disp();
        // update the counter
        count++;
    }
    
    printf("solver final order %d %d\n",dimorder[0],dimorder[1]);
    printf("solver final size %ld %ld\n",sizeorder[0],sizeorder[1]);
}

/**
 * @brief allocate the data associated to the solver
 * 
 * @param size  the size to allocate
 * @param data  the adress of the data to allocate
 */
void FFTW_Solver::_allocate_data(const size_t size[DIM],void** data)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // Sanity checks
    //-------------------------------------------------------------------------
    assert((*data) == NULL);

    //-------------------------------------------------------------------------
    // Do the memory allocation
    //-------------------------------------------------------------------------
    printf("are we complex? %d\n",_isComplex);
    size_t size_tot = 1;
    for(int id=0; id<DIM; id++) size_tot *= size[id];

    if(_isComplex){
        INFOLOG2("Complex memory allocation, size = %ld\n",size_tot);
        (*data) =(void*) fftw_malloc(size_tot*sizeof(fftw_complex));
    }
    else{
        INFOLOG2("Real memory allocation, size = %ld\n",size_tot);
        (*data) =(void*) fftw_malloc(size_tot*sizeof(double));
    }   
}

/**
 * @brief deallocate the data associated with the solver
 * 
 * @param data the data to deallocate
 */
void FFTW_Solver::_deallocate_data(void* data)
{
    BEGIN_FUNC
    if(data != NULL){
        fftw_free(data);
    }
}


void  FFTW_Solver::_allocate_plan(const size_t size[DIM],const size_t offset, const bool isComplex,void* data, multimap<int,FFTW_plan_dim* > *planmap)
{
    BEGIN_FUNC
    for(multimap<int,FFTW_plan_dim* >::iterator it = (*planmap).begin(); it != (*planmap).end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        // initialize the plan - read only
        myplan->allocate_plan(size,offset,isComplex,data);
    }
}


 /**
 * @brief Solve the Poisson equation nabla^2 field = rhs
 * 
 * @param field 
 * @param rhs 
 * 
 * -----------------------------------------------
 * We perform the following operations:
 */
void FFTW_Solver::solve(double* field, double* rhs)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    assert(field != NULL);
    assert(rhs   != NULL);
    
    {
        int mysize4[2] = {_size_field[0],_size_field[1]};
        write_array(mysize4,rhs,"rhs_in");
    }

    //-------------------------------------------------------------------------
    /** - clean the data memory */
    //-------------------------------------------------------------------------
    // get the data pointer to double style    
    double* in = (double*) _data;
    
    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){
                const int id = iz*_size_hat[1]*_size_hat[0] + iy * _size_hat[0] + ix;
                in[id] = 0.0 ;
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - copy the rhs in the correct order */
    //-------------------------------------------------------------------------
    INFOLOG ("------------------------------------------\n");
    INFOLOG ("## memory information\n")
    INFOLOG4("- size field   = %ld %ld %ld\n",_size_field[0],_size_field[1],_size_field[2]);
    INFOLOG4("- size hat     = %ld %ld %ld\n",_size_hat[0],_size_hat[1],_size_hat[2]);
    INFOLOG4("- dim order    = %d %d %d\n",_dimorder[0],_dimorder[1],_dimorder[2]);
    INFOLOG4("- field start  = %ld %ld %ld\n",_fieldstart[0],_fieldstart[1],_fieldstart[2]);
    INFOLOG4("- dim multfact = %ld %ld %ld\n",_dim_multfact[0],_dim_multfact[1],_dim_multfact[2]);
    INFOLOG2("- offset       = %ld\n",_offset);
    INFOLOG ("------------------------------------------\n");

    for(int iz=0; iz<_size_field[2]; iz++){
        for(int iy=0; iy<_size_field[1]; iy++){
            for(int ix=0; ix<_size_field[0]; ix++){
                // comnpute the index permutation
                const int id_field   = iz*_size_field[1]*_size_field[0] + iy * _size_field[0] + ix;
                const int id_fourier = iz*_dim_multfact[2] + iy*_dim_multfact[1] + ix*_dim_multfact[0] + _offset;
                // put the data
                in[id_fourier] = rhs[id_field] ;
            }
        }
    }

    
    if(_isComplex){
        int mysize[2] = {_size_hat[0]*2,_size_hat[1]};
        write_array(mysize,in,"rhs_pad");
    }
    else{
        int mysize[2] = {_size_hat[0],_size_hat[1]};
        write_array(mysize,in,"rhs_pad");
    }

    //-------------------------------------------------------------------------
    /** - go to Fourier */
    //-------------------------------------------------------------------------
    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_forward.begin(); it != _plan_forward.end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        // initialize the plan - read only
        myplan->execute_plan();
    }    

    //-------------------------------------------------------------------------
    /** - Perform the magic */
    //-------------------------------------------------------------------------
    printf("doing the magic\n");
    if(!_isComplex){
        double* mydata = (double*) _data;
        for(int iz=0; iz<_size_hat[2]; iz++){
            for(int iy=0; iy<_size_hat[1]; iy++){
                for(int ix=0; ix<_size_hat[0]; ix++){
                    const int id = iz*_size_hat[1]*_size_hat[0] + iy * _size_hat[0] + ix;
                    mydata[id] *= _normfact;
                }
            }
        }
    }
    else{
        fftw_complex* mydata = (fftw_complex*) _data;
        for(int iz=0; iz<_size_hat[2]; iz++){
            for(int iy=0; iy<_size_hat[1]; iy++){
                for(int ix=0; ix<_size_hat[0]; ix++){
                    const int id = iz*_size_hat[1]*_size_hat[0] + iy * _size_hat[0] + ix;
                    mydata[id][0] *= _normfact;
                    mydata[id][1] *= _normfact;
                }
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - go back to reals */
    //-------------------------------------------------------------------------
    for(multimap<int,FFTW_plan_dim* >::reverse_iterator it = _plan_backward.rbegin(); it != _plan_backward.rend(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        myplan->execute_plan();
    }

    if(_isComplex){
        int mysize2[2] = {_size_hat[0]*2,_size_hat[1]};
        write_array(mysize2,in,"sol_pad");
    }
    else{
        int mysize2[2] = {_size_hat[0],_size_hat[1]};
        write_array(mysize2,in,"sol_pad");
    }

    //-------------------------------------------------------------------------
    /** - copy the solution to the field */
    //-------------------------------------------------------------------------
    for(int iz=0; iz<_size_field[2]; iz++){
        for(int iy=0; iy<_size_field[1]; iy++){
            for(int ix=0; ix<_size_field[0]; ix++){
                // comnpute the index permutation
                const int id_field   = iz*_size_field[1]*_size_field[0] + iy * _size_field[0] + ix;
                const int id_fourier = iz*_dim_multfact[2] + iy*_dim_multfact[1] + ix*_dim_multfact[0] + _offset;
                // put the data
                field[id_field] = in[id_fourier] ;
            }
        }
    }

    {
        int mysize4[2] = {_size_field[0],_size_field[1]};
        write_array(mysize4,field,"sol_final");
    }
    
}