/**
 * @file FFTW_plan_dim.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#include "FFTW_plan_dim.hpp"


/**
 * @brief Construct a new FFTW_plan_dim object
 * 
 * @param dimID the dimension id in the non-transpose reference = the field reference
 * @param h the grid spacing
 * @param L the lenght of the computational domain
 * @param mybc the boundary condition to use for this plan
 * @param sign the sign of the plan (FFTW_FORWARD or FFTW_BACKWARD)
 * @param isGreen boolean to indicate if the plan is intended for Green's function
 */
FFTW_plan_dim::FFTW_plan_dim(const int dimID,const double h[DIM],const double L[DIM],const BoundaryType mybc[2],const int sign, const bool isGreen):
_dimID(dimID),
_sign(sign),
_isGreen(isGreen)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(dimID <  DIM);
    assert(dimID >= 0);
    
    //-------------------------------------------------------------------------
    // Initialisation of the sizes and types
    //-------------------------------------------------------------------------
    _bc[0] = mybc[0];
    _bc[1] = mybc[1];
    //determine the type of the solver
    int mytype = _bc[0]+_bc[1];
    //-------------------------------------------------------------------------
    // Get type and mult factors
    //-------------------------------------------------------------------------
    if( mytype <= R2R ){
        _type   = R2R;
        _normfact   = 1.0;
        _volfact  = 1.0; // no convolution so no multiplication by h
        _kfact = c_2pi/(2.0*L[_dimID]);
        if(_isGreen) _dospectral = true;
    }
    else if( mytype <= MIX ){
        _type   = MIX;
        _normfact   = 1.0;
        _volfact  = h[_dimID];
        _kfact = c_2pi/(4.0*L[_dimID]);
    }
    else if( mytype == PERPER ){
        _type   = PERPER;
        _normfact   = 1.0; 
        _volfact  = 1.0; // no convolution so no multiplication by h
        _kfact = c_2pi/(L[_dimID]);
        if(_isGreen) _dospectral = true;
    }
    else if( mytype == UNBUNB ){
        _type   = UNBUNB;
        _normfact = 1.0;
        _volfact  = h[_dimID];
        _kfact   = c_2pi/(2.0*L[_dimID]);
    }
}
FFTW_plan_dim::~FFTW_plan_dim(){
    BEGIN_FUNC
    if(_plan != NULL) fftw_destroy_plan(_plan);
}

/**
 * @brief Initialize the FFTW_plan_dim by performing a 'dry run'
 * 
 * The function redirects to one of the init functions depending on the type:
 * - _init_real2real()
 * - _init_mixpoisson()
 * - _init_periodic()
 * - _init_unbounded()
 * 
 * Each of the sub-function initializes the following variables
 * - #_n_in
 * - #_n_out
 * - #_fieldstart
 * - #_switch2Complex
 * - #_howmany (for non-Green functions only)
 * - #_imult
 * - #_kind (for R2R and MIX plans only)
 * - #_symstart (for Green's function only)
 * 
 * @param size the current size that will come in (hence already partially transformed)
 * @param isComplex the current complex state of the data
 */
void FFTW_plan_dim::init(const int size[DIM],const bool isComplex){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(size[_dimID] >= 0);

    //-------------------------------------------------------------------------
    // redirect to the corresponding subfunction
    //-------------------------------------------------------------------------
    if( _type == R2R ){
        _init_real2real(size,isComplex);
    }
    else if( _type == MIX ){
        // _n = 2*size[dimID]; // we have to double the size
        _init_mixpoisson(size,isComplex);
    }
    else if( _type == PERPER ){
        _init_periodic(size,isComplex);
    }
    else if( _type == UNBUNB ){
        _init_unbounded(size,isComplex);
    }
}

/**
 * @brief Initialize for a real to real plan
 * 
 * @param size 
 * @param isComplex 
 * 
 * -------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim::_init_real2real(const int size[DIM],const bool isComplex)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    assert(isComplex == false);
    
    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__switch2Complex)  */
    //-------------------------------------------------------------------------
    _n_in       = size[_dimID];
    _n_out      = size[_dimID];
    _fieldstart = 0;

    // no switch to complex
    _switch2Complex = false;
    _shiftgreen     = 0;

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0; // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - get #_howmany if we are not a Green function */
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _howmany = 1;
        for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
        for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];
    }

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 /(2.0*size[_dimID]);

    //-------------------------------------------------------------------------
    /** - Get the #_kind of Fourier transforms and #_imult */
    //-------------------------------------------------------------------------
    if(_isGreen){
        return;
    }
    else if(_bc[0] == EVEN){ // We have a DCT

        _imult = false; // we do NOT have to multiply by i=sqrt(-1)

        if(_bc[1] == EVEN){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT10; // DCT type II
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT01;
        }
        else if(_bc[1] == ODD){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT01; // DCT type III
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT10;
        }   
    }
    else if(_bc[0] == ODD){ // We have a DST

        _imult = true; // we DO have to multiply by -i=-sqrt(-1)

        if(_bc[1] == ODD){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_RODFT10; // DST type II
            if(_sign == FFTW_BACKWARD) _kind = FFTW_RODFT01;
        }
        else if(_bc[1] == EVEN){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT11; // DST type IV
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT11;
        }
    }
    else{
        UP_ERROR("unable to init the solver required\n")
    }
}

/**
 * @brief Initialize for a mix unbounded-symmetry plan
 * 
 * @param size 
 * @param isComplex 
 * 
 * ----------------------------------------
 * We do the following operations
 */
void FFTW_plan_dim::_init_mixpoisson(const int size[DIM],const bool isComplex){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    assert(isComplex == false);

    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart and #__switch2Complex)  */
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _n_in  =  2*size[_dimID];
        _n_out =  2*size[_dimID];
    }
    else if( _isGreen){
        _n_in  =  2*size[_dimID]+1;
        _n_out =  2*size[_dimID]+1;
    }

    _switch2Complex = false;

    if (_isGreen) _fieldstart = 0;
    else if(_bc[0] == UNB) _fieldstart = size[_dimID];   // padding to the left
    else if(_bc[1] == UNB) _fieldstart = 0;              // padding to the right

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0; // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - get #_howmany if we are not a Green function */
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _howmany = 1;
        for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
        for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];
    }

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0/(4.0*size[_dimID]);
    
    //-------------------------------------------------------------------------
    /** - Get the #_kind of Fourier transforms, #_imult and #_shiftgreen */
    //-------------------------------------------------------------------------
    if(_isGreen){
        _imult      =  false;
        _shiftgreen = 0;
        // The Green function is ALWAYS EVEN - EVEN
        if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT00; // DCT type I
        if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT00;
    }
    else{
        if((_bc[0] == EVEN && _bc[1] == UNB)||(_bc[0] == UNB && _bc[1] == EVEN)){ // We have a DCT - we are EVEN - EVEN
            _imult      = false;
            _shiftgreen = 0;
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT10; // DCT type II
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT01;
        }
        else if((_bc[0] == UNB && _bc[1] == ODD)||(_bc[0] == ODD && _bc[1] == UNB)){ // We have a DCT - we are EVEN - ODD
            _imult      = true;
            _shiftgreen = 1;
            if(_sign == FFTW_FORWARD ) _kind = FFTW_RODFT10; // DST type II
            if(_sign == FFTW_BACKWARD) _kind = FFTW_RODFT01;
        }
        else{
            UP_ERROR("unable to init the solver required\n")     
        }
    }
}

/**
 * @brief Initialize for a periodic plan
 * 
 * @param size 
 * @param isComplex
 * 
 * ----------------------------------------
 * We do the following operations
 */
void FFTW_plan_dim::_init_periodic(const int size[DIM],const bool isComplex){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__switch2Complex)  */
    //-------------------------------------------------------------------------
    if(isComplex){
        _n_in  = 2*size[_dimID]; // takes n complex, return n complex
        _n_out = 2*size[_dimID];

        _switch2Complex = false;
    }
    else{
        _n_in  = size[_dimID]; // takes n real
        _n_out = 2*(_n_in/2+1); // return n_in/2 + 1 complex

        _switch2Complex = true;
    }

    _fieldstart = 0;
    _shiftgreen = 0;

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0; // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - get #_howmany if we are not a Green function */
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _howmany = 1;
        for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
        for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];
    }

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0/(size[_dimID]);

    //-------------------------------------------------------------------------
    /** - Get the #_imult factor */
    //-------------------------------------------------------------------------
    _imult = false;
}

/**
 * @brief 
 * 
 * @param size 
 * @param isComplex 
 * 
 * --------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim::_init_unbounded(const int size[DIM],const bool isComplex){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__switch2Complex)  */
    //-------------------------------------------------------------------------
    if(isComplex){
        _n_in  = 4*size[_dimID]; // takes 2n complex, return 2n complex
        _n_out = 4*size[_dimID];

        _switch2Complex = false;
    }
    else{
        _n_in  = 2*size[_dimID]; // takes 2n real
        _n_out = 2*(_n_in/2+1); // return n_in/2 + 1 complex
        
        _switch2Complex = true;
    }

    _fieldstart = 0;
    _shiftgreen = 0;

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = size[_dimID];
    
    //-------------------------------------------------------------------------
    /** - get howmany if we are not a Green function */
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _howmany = 1;
        for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
        for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];
    }

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0/(2*size[_dimID]);

    //-------------------------------------------------------------------------
    /** - Get the #_imult */
    //-------------------------------------------------------------------------
    _imult = false;
}

/**
 * @brief allocate the plan based on the information computed by _init()
 * 
 * The function redirects to one of the init functions depending on the type:
 * - _allocate_plan_real()
 * - _allocate_plan_complex()
 * 
 * @param size_plan the size of the transposed data
 * @param isComplex if the transpoed data is complex or real
 * @param data the pointer to the transposed data (has to be allocated)
 */
void FFTW_plan_dim::allocate_plan(const int size_plan[DIM],double* data)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // allocate the plan
    //-------------------------------------------------------------------------
    if( _type == R2R || _type == MIX ){
        _allocate_plan_real(size_plan,data);
    }
    else if( _type == PERPER || _type == UNBUNB ){
        _allocate_plan_complex(size_plan,data);
    }
}

/**
 * @brief Allocate a plan that only treats real numbers
 * 
 * @note
 * The howmany is recomputed if we consider a Green function since it has been symmetrized
 * on the entire domain
 * 
 * @warning
 * If there is a r2c tranform we will perfom upto 2 spurious transform.
 * This is not important since the Green's transform is only performed once
 * 
 * @param size_ordered the size of the transposed data
 * @param offset the offset in memory computed by FFW_Solver in double indexing unit
 * @param data the pointer to the transposed data (has to be allocated)
 * 
 */
void FFTW_plan_dim::_allocate_plan_real(const int size_ordered[DIM], double *data)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    assert(data != NULL);

    //-------------------------------------------------------------------------
    /** - If is Green and #_type is R2R, exit */
    //-------------------------------------------------------------------------
    if(_isGreen && _type == R2R){
        _plan = NULL;
        
        INFOLOG("------------------------------------------\n");
        INFOLOG("## no real to real plan created for Green\n");
        INFOLOG ("------------------------------------------\n");
        return;
    }

    //-------------------------------------------------------------------------
    /** - Get compute the rank and the stides  */
    //-------------------------------------------------------------------------
    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    int rank = 1;

    // set the sizemult array usefull for strides and dists
    // int sizemult[DIM]; sizemult[0] = 1;
    // for(int id=1; id<DIM; id++) sizemult[id]=sizemult[id-1]*size_ordered[id-1];
    // if(_isDataComplex) for(int id=1; id<DIM; id++) sizemult[id]*= 2;

    // int istride = sizemult[(_orderID  )%DIM];
    // int idist   = sizemult[(_orderID+1)%DIM];
    // int ostride = sizemult[(_orderID  )%DIM];
    // int odist   = sizemult[(_orderID+1)%DIM];
    int istride = 1;
    int ostride = 1;
    int idist   = size_ordered[_dimID];
    int odist   = size_ordered[_dimID];
    

    //-------------------------------------------------------------------------
    /** - If is Green, compute #_howmany  */
    //-------------------------------------------------------------------------
    if(_isGreen){
        _howmany = 1;
        for(int id=0         ; id<_orderID; id++) _howmany *= size_ordered[id];
        for(int id=_orderID+1; id<DIM     ; id++) _howmany *= size_ordered[id];
        // if the data are complex we need to count twice the fastest rotating index
        // if(_isDataComplex && _orderID>0) _howmany *= 2;
    }

    //-------------------------------------------------------------------------
    /** - Create the plan  */
    //-------------------------------------------------------------------------
    _plan = fftw_plan_many_r2r( rank,(int*) (&_n_in),_howmany,
                                data,NULL,istride,idist,
                                data,NULL,ostride,odist,&_kind,FFTW_FLAG);
    
    INFOLOG ("------------------------------------------\n");
    if      (_type == R2R   ) {INFOLOG2("## R2R plan created for plan r2r (=%d)\n",_type);}
    else if (_type == MIX   ) {INFOLOG2("## R2R plan created for plan mix (=%d)\n",_type);}
    INFOLOG2("orderedID = %d\n",_orderID);
    if      (DIM == 2) {INFOLOG3("size = %d x %d\n",size_ordered[0],size_ordered[1]);}
    else if (DIM == 3) {INFOLOG4("size = %d x %d x %d\n",size_ordered[0],size_ordered[1],size_ordered[2]);}
    INFOLOG2("howmany   = %d\n",_howmany);
    INFOLOG2("size n    = %d\n",_n_in);
    INFOLOG3("istride (double) = %d - idist = %d\n",istride,idist);
    INFOLOG3("ostride (double) = %d - odist = %d\n",ostride,odist);
    if(_type == R2R) {INFOLOG2("starting offset  = %ld\n",offset);}
    INFOLOG ("------------------------------------------\n");

}

/**
 * @brief allocate a plan that treats complex numbers (r2c or c2c)
 * 
 * @note
 * The howmany is recomputed if we consider a Green function since it has been symmetrized
 * on the entire domain
 * 
 * @warning
 * If there is a r2c tranform we will perfom upto 2 spurious transform.
 * This is not important since the Green's transform is only performed once
 * 
 * 
 * @param size_ordered 
 * @param data 
 */
void FFTW_plan_dim::_allocate_plan_complex(const int size_ordered[DIM],double* data){
    BEGIN_FUNC

    assert(data != NULL);

    if(_isGreen && _type == PERPER){
        _plan = NULL;
        
        INFOLOG("------------------------------------------\n");
        INFOLOG("## no DFT plan created for Green\n");
        INFOLOG ("------------------------------------------\n");
        return;
    }

    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    int rank = 1;

    // if we are green we need to recompute howmany
    if(_isGreen){
        _howmany = 1;
        for(int id=0         ; id<_orderID; id++) _howmany *= size_ordered[id];
        for(int id=_orderID+1; id<DIM     ; id++) _howmany *= size_ordered[id];
    }

    // set the sizemult array usefull for strides and dists
    // int sizemult[DIM]; sizemult[0] = 1;
    // for(int id=1; id<DIM; id++) sizemult[id]=sizemult[id-1]*size_ordered[id-1];
    

    // ostrid = complex
    int istride = 1;
    int ostride = 1;
    int idist   = size_ordered[_dimID];
    int odist   = size_ordered[_dimID];
    
    
    // incomming arrays depends if we are a complex switcher or not
    if(_switch2Complex){
        // if the data are complex we need to count twice the fastest rotating index
        // if(_isDataComplex) for(int id=1; id<DIM; id++) sizemult[id]*= 2;
        // if(_isDataComplex && _orderID>0) _howmany *= 2;

        UP_CHECK0(odist%2==0,"odist has to be multiple of 2");
        odist = odist/2;
        

        INFOLOG ("------------------------------------------\n");
        if      (_type == PERPER) {INFO2("## R2C plan created for plan periodic-periodic (=%d)\n",_type);}
        else if (_type == UNBUNB) {INFO2("## R2C plan created for plan unbounded (=%d)\n",_type);}
        INFOLOG2("orderedID = %d\n",_orderID);
        if      (DIM == 2) {INFOLOG3("size = %d x %d\n",size_ordered[0],size_ordered[1]);}
        else if (DIM == 3) {INFOLOG4("size = %d x %d x %d\n",size_ordered[0],size_ordered[1],size_ordered[2]);}
        INFOLOG2("howmany   = %d\n",_howmany);
        INFOLOG2("size n    = %d\n",_n_in);
        INFOLOG3("istride (double)  = %d - idist = %d\n",istride,idist);
        INFOLOG3("ostride (complex) = %d - odist = %d\n",ostride,odist);
        INFOLOG ("------------------------------------------\n");

        // set the plan - there is no offset in r2c or c2c possible
        if(_sign == FFTW_FORWARD){
            _plan = fftw_plan_many_dft_r2c( rank,(int*) (&_n_in),_howmany,
                                            data,NULL,istride,idist,
                                            (fftw_complex*) data,NULL,ostride,odist,FFTW_FLAG);
        }
        else{
            _plan = fftw_plan_many_dft_c2r( rank,(int*) (&_n_in),_howmany,
                                            (fftw_complex*) data,NULL,ostride,odist,
                                            data,NULL,istride,idist,FFTW_FLAG);
        }

        
    }
    else{
        // set the plan
        _plan = fftw_plan_many_dft( rank,(int*) (&_n_in),_howmany,
                                    (fftw_complex*) data,NULL,istride,idist,
                                    (fftw_complex*) data,NULL,ostride,odist,_sign,FFTW_FLAG);

        INFOLOG ("------------------------------------------\n");
        if      (_type == PERPER) {INFO2("## C2C plan created for plan periodic-periodic (=%d)\n",_type);}
        else if (_type == UNBUNB) {INFO2("## C2C plan created for plan unbounded (=%d)\n",_type);}
        INFOLOG2("orderedID = %d\n",_orderID);
        if      (DIM == 2) {INFOLOG3("size = %d x %d\n",size_ordered[0],size_ordered[1]);}
        else if (DIM == 3) {INFOLOG4("size = %d x %d x %d\n",size_ordered[0],size_ordered[1],size_ordered[2]);}
        INFOLOG2("howmany = %d\n",_howmany);
        INFOLOG2("size n = %d\n",_n_in);
        INFOLOG3("istride (complex) = %d - idist = %d\n",istride,idist);
        INFOLOG3("ostride (complex) = %d - odist = %d\n",ostride,odist);
        INFOLOG ("------------------------------------------\n");
    }
}

/**
 * @brief Executes the plan
 * 
 */
void FFTW_plan_dim::execute_plan(){
    BEGIN_FUNC
    // run the plan
    if      (_type == R2R   ) {INFO2(">> Doing plan real2real for order %d\n",_orderID);}
    else if (_type == MIX   ) {INFO2(">> Doing plan mix for order %d\n",_orderID);}
    else if (_type == PERPER) {INFO2(">> Doing plan periodic-periodic for order %d\n",_orderID);}
    else if (_type == UNBUNB) {INFO2(">> Doing plan unbounded for order %d\n",_orderID);}
    fftw_execute(_plan);
}

// /**
//  * @brief 
//  * 
//  * @param id 
//  * @param size 
//  */
// void FFTW_plan_dim::get_outsize (const int id, int size[DIM]) const {
//     BEGIN_FUNC
//     size[id] = _n_out;
// }
// /**
//  * @brief returns the size that comes into the FFT plan
//  * 
//  * @param id 
//  * @param size 
//  */
// void FFTW_plan_dim::get_outsize_double (const int id, int size[DIM]) const {
//     BEGIN_FUNC
//     if(_switch2Complex) size[id] = _n_out*2;
//     else                size[id] = _n_out;
// }
// /**
//  * @brief 
//  * 
//  * @param id 
//  * @param start 
//  */
// void FFTW_plan_dim::get_fieldstart(const int id, int start[DIM]) const {
//     BEGIN_FUNC
//     start[id] = _fieldstart;
// }
// /**
//  * @brief 
//  * 
//  * @param id 
//  * @param dimID 
//  */
// void FFTW_plan_dim::get_dimID   (const int id, int dimID[DIM]) const {
//     BEGIN_FUNC
//     dimID[id] = _dimID;
// }
// /**
//  * @brief 
//  * 
//  * @param id 
//  */
// void FFTW_plan_dim::set_order(const int id){
//     BEGIN_FUNC
//     _orderID = id;
// }

/**
 * @brief display the FFTW_plan_dim object
 * 
 */
void FFTW_plan_dim::disp(){
    BEGIN_FUNC
    INFO ("------------------------------------------\n");
    INFO3("## Plan num %d created for dimension %d\n",_orderID,_dimID);
    if      (_type == R2R   ) {INFO2("- type = real2real (=%d)\n",_type);}
    else if (_type == MIX   ) {INFO2("- type = mix (=%d)\n",_type);}
    else if (_type == PERPER) {INFO2("- type = periodic-periodic (=%d)\n",_type);}
    else if (_type == UNBUNB) {INFO2("- type = unbounded (=%d)\n",_type);}
    if      (_bc[0] == EVEN) {INFO("- bc = { EVEN ,");}
    else if (_bc[0] == ODD)  {INFO("- bc = { ODD  ,");}
    else if (_bc[0] == UNB)  {INFO("- bc = { UNB  ,");}
    else if (_bc[0] == PER)  {INFO("- bc = { PER  ,");}
    if      (_bc[1] == EVEN) {INFO(" EVEN}\n");}
    else if (_bc[1] == ODD)  {INFO(" ODD}\n");}
    else if (_bc[1] == UNB)  {INFO(" UNB}\n");}
    else if (_bc[1] == PER)  {INFO(" PER}\n");}
    if (_type == R2R || _type == MIX ){
        if(_kind == FFTW_REDFT00) {INFO("- kind = REDFT00 = DCT type I\n")}
        if(_kind == FFTW_REDFT10) {INFO("- kind = REDFT10 = DCT type II\n")}
        if(_kind == FFTW_REDFT01) {INFO("- kind = REDFT01 = DCT type III\n")}
        if(_kind == FFTW_REDFT11) {INFO("- kind = REDFT11 = DCT type IV\n")}
        if(_kind == FFTW_RODFT00) {INFO("- kind = RODFT00 = DST type I\n")}
        if(_kind == FFTW_RODFT10) {INFO("- kind = RODFT10 = DST type II\n")}
        if(_kind == FFTW_RODFT01) {INFO("- kind = RODFT01 = DST type III\n")}
        if(_kind == FFTW_RODFT11) {INFO("- kind = RODFT11 = DST type IV\n")}
    }
    INFO2("- orderID    = %d\n",_orderID);
    INFO2("- is Green   ? %d\n",_isGreen);
    INFO2("- s2Complex  ? %d\n",_switch2Complex);
    INFO2("- n_in       = %d\n",_n_in);
    INFO2("- n_out      = %d\n",_n_out);
    INFO2("- howmany    = %d\n",_howmany);
    INFO2("- fieldstart = %d\n",_fieldstart);
    INFO2("- shiftgreen = %d\n",_shiftgreen);
    INFO2("- dospectral ? %d\n",_dospectral);
    if      (_sign == FFTW_FORWARD ) {INFO("- FORWARD plan\n");}
    else if (_sign == FFTW_BACKWARD) {INFO("- BACKWARD plan\n");}
    
    INFO ("------------------------------------------\n");
}