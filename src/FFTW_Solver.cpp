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
 * 
 * --------------------------------------
 * We do the following operations:
 */
FFTW_Solver::FFTW_Solver(const int size_field[DIM],const double h[DIM],const double L[DIM],const BoundaryType mybc[DIM][2])
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Initialize MPI if needed */
    //-------------------------------------------------------------------------
    #if (DIM == 3)
        fftw_mpi_init();
    #endif

    //-------------------------------------------------------------------------
    /** - Store the field size */
    //-------------------------------------------------------------------------
    for(int id=0; id<DIM; id++) _size_field[id] = size_field[id];

    //-------------------------------------------------------------------------
    /** - For each dim, compute the plan and its type */
    //-------------------------------------------------------------------------
    for(int id=0; id<DIM; id++)
    {
        // forward
        FFTW_plan_dim* myplan_forward = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_FORWARD,false);
        _plan_forward. insert(pair<int,FFTW_plan_dim*>(myplan_forward->get_type(),myplan_forward));
        // backward
        FFTW_plan_dim* myplan_backward = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_BACKWARD,false);
        _plan_backward.insert(pair<int,FFTW_plan_dim*>(myplan_backward->get_type(),myplan_backward));

        // Green's function
        FFTW_plan_dim* myplan_green = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_FORWARD,true);
        _plan_green.insert(pair<int,FFTW_plan_dim*>(myplan_green->get_type(),myplan_green));
    }

    //-------------------------------------------------------------------------
    /** - Initialise the plans and get the sizes + dimorder */
    //-------------------------------------------------------------------------
    // forward, store the size and dim order for the object
    _init_plan_map(_size_hat,_fieldstart,_dimorder,&_isComplex,&_plan_forward);

    // backward uses temporary sizes and check the output
    int size_tmp      [3] = {1,1,1};
    int fieldstart_tmp[3] = {0,0,0};
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

    // backward uses temporary sizes and check the output
    _init_plan_map(_size_hat_green,fieldstart_tmp,dimorder_tmp,&isComplex_tmp,&_plan_green);
    // sanity checks
    for(int id=0; id<3; id++){
        assert(fieldstart_tmp[id] == 0); // no field start for Green
        assert(dimorder_tmp[id]   == _dimorder[id]);
        assert(isComplex_tmp      == _isComplex);
    }

    //-------------------------------------------------------------------------
    /** - Get the normalization factors and the grid spacing and the #_shiftgreen factor */
    //-------------------------------------------------------------------------
    _normfact  = 1.0;
    _volfact   = 1.0;
    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_forward.begin(); it != _plan_forward.end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        _normfact *= myplan->get_normfact();
        _volfact  *= myplan->get_volfact ();
        // get the order ID
        const int orderID    = myplan->get_order();
        _hgrid[orderID]      = h[_dimorder[orderID]];
        _shiftgreen[orderID] = myplan->get_shiftgreen();
    }

    //-------------------------------------------------------------------------
    /** - Get the imult factor */
    //-------------------------------------------------------------------------
    _nbr_imult = 0;
    bool imult_forward [DIM] = {false};
    bool imult_backward[DIM] = {false};
    bool imult_green   [DIM] = {false};
    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_forward.begin(); it != _plan_forward.end(); ++it)
    {
        FFTW_plan_dim* myplan   = it->second;
        const int orderID       = myplan->get_order();
        imult_forward[orderID]  = myplan->get_imult();
    }
    // store the number of imult
    for(int id=0; id<DIM; id++) if(imult_forward[id]) _nbr_imult++;

    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_backward.begin(); it != _plan_backward.end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        const int orderID       = myplan->get_order();
        imult_backward[orderID] = myplan->get_imult();
    }
    // store the number of imult
    for(int id=0; id<DIM; id++) if(imult_backward[id]) _nbr_imult--;

    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_green.begin(); it != _plan_green.end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        const int orderID       = myplan->get_order();
        imult_green[orderID] = myplan->get_imult();
    }
    // store the number of imult
    for(int id=0; id<DIM; id++) if(imult_green[id]) _nbr_imult++;
}

/**
 * @brief Sets up the Solver
 * 
 * After this function the parameter of the solver (size etc) cannot be changed anymore
 * 
 * -------------------------------------------
 * We do the following operations
 */
void FFTW_Solver::setup(const SolverType mytype)
{
    //-------------------------------------------------------------------------
    /** - Store the solver type */
    //-------------------------------------------------------------------------
    _type = mytype;

    UP_CHECK1(mytype == UP_SRHS,"unsupported solver type %d\n",mytype);

    //-------------------------------------------------------------------------
    /** - Store some usefull factors in 'double' index calculus */
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
    /** - allocate the data */
    //-------------------------------------------------------------------------
    _allocate_data(_size_hat,&_data);
    _allocate_data(_size_hat_green,&_green);

    //-------------------------------------------------------------------------
    /** - allocate the plans */
    //-------------------------------------------------------------------------
    _allocate_plan(_size_hat,_offset,_isComplex,_data,&_plan_forward );
    _allocate_plan(_size_hat,_offset,_isComplex,_data,&_plan_backward);

    //-------------------------------------------------------------------------
    /** - allocate the plan and comnpute the Green's function */
    //-------------------------------------------------------------------------
    _allocate_plan(_size_hat_green,0,_isComplex,_green,&_plan_green);
    _compute_Green(_size_hat_green,_green,&_plan_green);

    //-------------------------------------------------------------------------
    /** - do the forward transfrom*/
    //-------------------------------------------------------------------------
    for(multimap<int,FFTW_plan_dim* >::iterator it = _plan_green.begin(); it != _plan_green.end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        myplan->execute_plan();
    }

    // {
    //     int mysize4[2] = {_size_hat_green[0],_size_hat_green[1]};
    //     write_array(mysize4,(fftw_complex*) _green,"new_green_fourier");
    // }
    
    //-------------------------------------------------------------------------
    /** - delete the plan */
    //-------------------------------------------------------------------------
    _delete_plan(&_plan_green);

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
    _deallocate_data(_green);

    //cleanup
    fftw_cleanup();
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
 * @brief initialize a multimap containing plan
 * 
 * @param[out] sizeorder the size of the tranformed field in the correct order
 * @param[out] dimorder the order of 
 * @param[out] isComplex indicate if the data needed by the plan is complex or not
 * @param[in/out] planmap the 2 (or 3) plans for each direction to initialize
 */
void FFTW_Solver::_init_plan_map(int sizeorder[3], int fieldstart[3], int dimorder[3], bool* isComplex, multimap<int,FFTW_plan_dim* > *planmap)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // get the plan info
    //-------------------------------------------------------------------------
    (*isComplex) = false;
    // get the sizes
    int size_tmp[DIM];
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
}

/**
 * @brief allocate the data associated to the solver
 * 
 * @param size  the size to allocate
 * @param data  the adress of the data to allocate
 */
void FFTW_Solver::_allocate_data(const int size[3],double** data)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0((*data) == NULL,"Pointer has to be NULL for allocation");

    //-------------------------------------------------------------------------
    /** - Do the memory allocation */
    //-------------------------------------------------------------------------
    size_t size_tot = 1;
    for(int id=0; id<DIM; id++) size_tot *= size[id];

    if(_isComplex){
        INFOLOG2("Complex memory allocation, size = %ld\n",size_tot);
        (*data) =(double*) fftw_malloc(size_tot*sizeof(fftw_complex));
    }
    else{
        INFOLOG2("Real memory allocation, size = %ld\n",size_tot);
        (*data) =(double*) fftw_malloc(size_tot*sizeof(double));
    }

    //-------------------------------------------------------------------------
    /** - Check memory alignement */
    //-------------------------------------------------------------------------
    UP_CHECK1(UP_ISALIGNED(*data),"FFTW alignement not compatible with UP_ALIGNMENT (=%d)",UP_ALIGNMENT);
}

/**
 * @brief deallocate the data associated with the solver
 * 
 * @param data the data to deallocate
 */
void FFTW_Solver::_deallocate_data(double* data)
{
    BEGIN_FUNC
    if(data != NULL){
        fftw_free(data);
    }
}


void  FFTW_Solver::_allocate_plan(const int size[3],const size_t offset, const bool isComplex,double* data, multimap<int,FFTW_plan_dim* > *planmap)
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
 * @brief compute the green function and fill the green data
 * 
 * @param size_green 
 * @param green 
 * @param planmap 
 * 
 * @warning
 * We to fill the green function using double unit indexing, so we first require the size in double
 * 
 * --------------------------------------
 * We do the following operations:
 */
void FFTW_Solver::_compute_Green(const int size_green[3], double* green, multimap<int,FFTW_plan_dim* >* planmap){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    assert(green != NULL);

    //-------------------------------------------------------------------------
    /** - get the direction where we need to do spectral diff and count them */
    //-------------------------------------------------------------------------
    int  order[3];
    bool dospectral[3] = {false};

    double hfact   [3] = {_hgrid[0],_hgrid[1],_hgrid[2]};
    double kfact   [3] = {0.0};
    int symstart[3] = {0};

    for(multimap<int,FFTW_plan_dim* >::iterator it = (*planmap).begin(); it != (*planmap).end(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;

        const int id   = myplan->get_order();
        dospectral[id] = myplan->get_dospectral();
        symstart[id]   = myplan->get_symstart();

        if(dospectral[id]){
            hfact[id] = 0.0;
            kfact[id] = myplan->get_kfact();
        }
    }

    // count the number of spectral dimensions
    int nbr_spectral = 0;
    for(int id=0; id<3; id++) if(dospectral[id]) nbr_spectral++;

    //-------------------------------------------------------------------------
    /** - get the green size in double unit indexing to fill it, we double only the fastest rotating index*/
    //-------------------------------------------------------------------------
    int dsize_green[3] = {size_green[0],size_green[1],size_green[2]};
    if(_isComplex) dsize_green[0] *= 2;

    //-------------------------------------------------------------------------
    /** - get the expression of Green in the full domain*/
    //-------------------------------------------------------------------------
    
    if     (nbr_spectral == 0){
        INFOLOG(">> using Green function 3 dir unbounded\n");
        if(DIM==2){
            if     (_greenorder == CHAT_2) Green_2D_3dirunbounded_0dirspectral_chat2(dsize_green,hfact,green);
            else if(_greenorder == HEJ_2 ) Green_2D_3dirunbounded_0dirspectral_hej2 (dsize_green,hfact,green,_greenalpha);
            else if(_greenorder == HEJ_4 ) Green_2D_3dirunbounded_0dirspectral_hej4 (dsize_green,hfact,green,_greenalpha);
            else UP_CHECK2(false,"Green Function = %d  unknow for nbr_spectral = %d",_greenorder,nbr_spectral);
        }
        else if(DIM==3){
            // if     (_greenorder == CHAT_2) Green_3D_3dirunbounded_0dirspectral_chat2(dsize_green,hfact,green);
            // else if(_greenorder == HEJ_2 ) Green_3D_3dirunbounded_0dirspectral_hej2 (dsize_green,hfact,green,_greenalpha);
            // else if(_greenorder == HEJ_4 ) Green_3D_3dirunbounded_0dirspectral_hej4 (dsize_green,hfact,green,_greenalpha);
            // else UP_CHECK2(false,"Green Function = %d  unknow for nbr_spectral = %d",_greenorder,nbr_spectral);
        }
    }
    else if(nbr_spectral == 1){
        INFOLOG(">> using Green function 2 dir unbounded - 1 dir spectral\n");
        // _compute_Green_2dirunbounded_1dirspectral
    } 
    else if(nbr_spectral == 2){
        INFOLOG(">> using Green function 1 dir unbounded - 2 dir spectral\n");
        // _compute_Green_1dirunbounded_2dirspectral
    }
    else if(nbr_spectral == 3){
        INFOLOG(">> using Green function 3 dir spectral\n");
        // _compute_Green_0dirunbounded_3dirspectral
    }

    //-------------------------------------------------------------------------
    /** - do the symmetry */
    //-------------------------------------------------------------------------

    // {int mysize[2] = {dsize_green[0],dsize_green[1]};
    // write_array(mysize,green,"green_ext");}

    // check if we have to symmetrize a direction

    for(int i2=0; i2<dsize_green[2]; i2++){
        for(int i1=0; i1<dsize_green[1]; i1++){
            for(int i0=0; i0<dsize_green[0]; i0++){

                const size_t id = i0 + dsize_green[0]*(i1 + dsize_green[1]*i2);
                // we have to take the symmetry around symstart: symstart - (iy - symstart) = 2 symstart - iy
                // to use the abs we have to go back to integers
                const int is0 = (symstart[0]==0 || i0 <= symstart[0]) ? i0 : abs(2*(int)symstart[0]-i0);
                const int is1 = (symstart[1]==0 || i1 <= symstart[1]) ? i1 : abs(2*(int)symstart[1]-i1);
                const int is2 = (symstart[2]==0 || i2 <= symstart[2]) ? i2 : abs(2*(int)symstart[2]-i2);

                const size_t id_sym = is0 + dsize_green[0]*(is1 + dsize_green[1]*is2);

                UP_CHECK3(id >= 0,"ID is not positive => id = %d, %d, %d",i0,i1,i2);
                UP_CHECK3(id < dsize_green[0]*dsize_green[1]*dsize_green[2],"ID is greater than the max size => id = %d, %d, %d",i0,i1,i2);

                green[id] = green[id_sym];
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - scale the Green data using #_volfact */
    //-------------------------------------------------------------------------
    for(int i2=0; i2<dsize_green[2]; i2++){
        for(int i1=0; i1<dsize_green[1]; i1++){
            for(int i0=0; i0<dsize_green[0]; i0++){
                const size_t id = i0 + dsize_green[0]*(i1 + dsize_green[1]*i2);
                green[id] = green[id] * _volfact;
            }
        }
    }
}

/**
 * @brief sets the Green type #_greenorder, see #OrderDiff
 * 
 * @param order 
 */
void FFTW_Solver::set_GreenType(const OrderDiff order){
    BEGIN_FUNC
    _greenorder = order;
}
/**
 * @brief set the Green spectral order of accuracy #_greendiff, see #OrderDiff
 * 
 * @param diff 
 */
void FFTW_Solver::set_GreenDiff(const OrderDiff diff){
    BEGIN_FUNC
    _greendiff = diff;
}

/**
 * @brief set the alpha parameter for regularization, #_greenalpha
 * 
 * @param alpha 
 */
void FFTW_Solver::set_alpha(const double alpha){
    BEGIN_FUNC
    _greenalpha = alpha;
}


/**
 * @brief Solve the Poisson equation
 * 
 * @param field 
 * @param rhs 
 * 
 * -----------------------------------------------
 * We perform the following operations:
 */
void FFTW_Solver::solve(double * field, double * rhs)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0(field != NULL,"field is NULL");
    UP_CHECK0(rhs   != NULL,"rhs is NULL");
    UP_CHECK1(UP_ISALIGNED(field),"pointer no aligned to UP_ALIGNMENT (=%d)",UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(rhs  ),"pointer no aligned to UP_ALIGNMENT (=%d)",UP_ALIGNMENT);
    

    opt_double_ptr myfield     = field;
    opt_double_ptr mydata      = (double*) _data;
    const opt_double_ptr myrhs = rhs;

    //-------------------------------------------------------------------------
    /** - clean the data memory */
    //-------------------------------------------------------------------------
    // get the data pointer to double style    
   

    // __assume_aligned(mydata, UP_ALIGNMENT);

    // #pragma vector aligned
    // for(int iz=0; iz<_size_hat[2]; iz++){
    //     for(int iy=0; iy<_size_hat[1]; iy++){
    //         for(int ix=0; ix<_size_hat[0]; ix++){
    //             const int id = ix + _size_hat[0] *( iy + _size_hat[1]* iz);
    //             mydata[id] = 0;
    //         }
    //     }
    // }
    if(_isComplex) std::memset(mydata,0,sizeof(fftw_complex)*_size_hat[0]*_size_hat[1]*_size_hat[2]);
    else           std::memset(mydata,0,sizeof(double      )*_size_hat[0]*_size_hat[1]*_size_hat[2]);

    //-------------------------------------------------------------------------
    /** - copy the rhs in the correct order */
    //-------------------------------------------------------------------------
    INFOLOG ("------------------------------------------\n");
    INFOLOG ("## memory information\n")
    INFOLOG4("- size field   = %d %d %d\n",_size_field[0],_size_field[1],_size_field[2]);
    INFOLOG4("- size hat     = %d %d %d\n",_size_hat[0],_size_hat[1],_size_hat[2]);
    INFOLOG4("- dim order    = %d %d %d\n",_dimorder[0],_dimorder[1],_dimorder[2]);
    INFOLOG4("- field start  = %d %d %d\n",_fieldstart[0],_fieldstart[1],_fieldstart[2]);
    INFOLOG4("- dim multfact = %d %d %d\n",_dim_multfact[0],_dim_multfact[1],_dim_multfact[2]);
    INFOLOG2("- offset       = %ld\n",_offset);
    INFOLOG ("------------------------------------------\n");

    __assume_aligned(mydata, UP_ALIGNMENT);
    __assume_aligned(myrhs , UP_ALIGNMENT);

    for(int iz=0; iz<_size_field[2]; iz++){
        for(int iy=0; iy<_size_field[1]; iy++){
            for(int ix=0; ix<_size_field[0]; ix++){
                // comnpute the index permutation
                const int id_field   = ix + _size_field[0] * (iy + _size_field[1] * iz);
                const int id_fourier = iz*_dim_multfact[2] + iy*_dim_multfact[1] + ix*_dim_multfact[0] + _offset;
                // put the data
                mydata[id_fourier] = myrhs[id_field] ;
            }
        }
    }

    
    // if(_isComplex){
    //     int mysize[2] = {_size_hat[0]*2,_size_hat[1]};
    //     write_array(mysize,in,"rhs_pad");
    // }
    // else{
    //     int mysize[2] = {_size_hat[0],_size_hat[1]};
    //     write_array(mysize,in,"rhs_pad");
    // }

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
    if(_type == UP_SRHS){
        if(!_isComplex){
            if     (_nbr_imult == 0) dothemagic_rhs_real();
            else UP_CHECK1(false,"the number of imult = %d is not supported\n",_nbr_imult);
        }
        else{
            if     (_nbr_imult == 0) dothemagic_rhs_complex_nmult0();
            else if(_nbr_imult == 1) dothemagic_rhs_complex_nmult1();
            else if(_nbr_imult == 2) dothemagic_rhs_complex_nmult2();
            else if(_nbr_imult == 3) dothemagic_rhs_complex_nmult3();
            else UP_CHECK1(false,"the number of imult = %d is not supported\n",_nbr_imult);
        }
    }

    // {
    //     int mysize4[2] = {_size_hat[0],_size_hat[1]};
    //     write_array(mysize4,(fftw_complex*) _data,"new_field_fourier");
    // }

    //-------------------------------------------------------------------------
    /** - go back to reals */
    //-------------------------------------------------------------------------
    for(multimap<int,FFTW_plan_dim* >::reverse_iterator it = _plan_backward.rbegin(); it != _plan_backward.rend(); ++it)
    {
        FFTW_plan_dim* myplan = it->second;
        myplan->execute_plan();
    }

    // if(_isComplex){
    //     int mysize2[2] = {_size_hat[0]*2,_size_hat[1]};
    //     write_array(mysize2,in,"sol_pad");
    // }
    // else{
    //     int mysize2[2] = {_size_hat[0],_size_hat[1]};
    //     write_array(mysize2,in,"sol_pad");
    // }

    //-------------------------------------------------------------------------
    /** - copy the solution in the field */
    //-------------------------------------------------------------------------
    __assume_aligned(mydata , UP_ALIGNMENT);
    __assume_aligned(myfield, UP_ALIGNMENT);
    for(int iz=0; iz<_size_field[2]; iz++){
        for(int iy=0; iy<_size_field[1]; iy++){
            for(int ix=0; ix<_size_field[0]; ix++){
                // comnpute the index permutation
                const int id_field   = ix + _size_field[0] * (iy + _size_field[1] * iz);
                const int id_fourier = iz*_dim_multfact[2] + iy*_dim_multfact[1] + ix*_dim_multfact[0] + _offset;
                // put the data
                myfield[id_field] = mydata[id_fourier] ;
            }
        }
    }
}


/**
 * @brief perform the convolution for real to real cases
 * 
 */
void FFTW_Solver::dothemagic_rhs_real()
{
    opt_double_ptr mydata  = _data;
    const opt_double_ptr mygreen = _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                __assume_aligned(mydata,  UP_ALIGNMENT);
                __assume_aligned(mygreen, UP_ALIGNMENT);

                const size_t id       = ix + _size_hat[0] * (iy + _size_hat[1]* iz);
                const size_t id_green = ix + _shiftgreen[0] + _size_hat_green[0] * ( iy+_shiftgreen[1] + _size_hat_green[1] * (iz+_shiftgreen[2]));

                mydata[id] = mydata[id] * _normfact * mygreen[id_green];
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult0()
{
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                __assume_aligned(mydata,  UP_ALIGNMENT);
                __assume_aligned(mygreen, UP_ALIGNMENT);

                const size_t id       = ix + _size_hat[0] * (iy + _size_hat[1] * iz);
                const size_t id_green = (ix+_shiftgreen[0]) + _size_hat_green[0] * ((iy+_shiftgreen[1]) + _size_hat_green[1] * (iz+_shiftgreen[2]));

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                // update the values
                mydata[id][0] = _normfact * temp_real;
                mydata[id][1] = _normfact * temp_imag;
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult1()
{
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                __assume_aligned(mydata,  UP_ALIGNMENT);
                __assume_aligned(mygreen, UP_ALIGNMENT);

                const size_t id       = ix + _size_hat[0]       * (iy + _size_hat[1]       * iz);
                const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                mydata[id][0] = (-1.0) * _normfact * ( temp_imag );
                mydata[id][1] =          _normfact * ( temp_real );
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-1)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult2()
{
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                __assume_aligned(mydata,  UP_ALIGNMENT);
                __assume_aligned(mygreen, UP_ALIGNMENT);

                const size_t id       = ix + _size_hat[0]       * (iy + _size_hat[1]       * iz);
                const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                mydata[id][0] = (-1.0) * _normfact * ( temp_real );
                mydata[id][1] = (-1.0) * _normfact * ( temp_imag );
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult3()
{
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                __assume_aligned(mydata,  UP_ALIGNMENT);
                __assume_aligned(mygreen, UP_ALIGNMENT);    

                const size_t id       = ix + _size_hat[0]       * (iy + _size_hat[1]       * iz);
                const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                mydata[id][0] =          _normfact * ( temp_real );
                mydata[id][1] = (-1.0) * _normfact * ( temp_imag );
            }
        }
    }
}