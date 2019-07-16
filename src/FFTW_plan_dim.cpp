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
        _fact   = 1.0; // no convolution so no multiplication by h
        _k_fact = C2PI/(2.0*L[_dimID]);
    }
    else if( mytype <= MIX ){
        _type   = MIX;
        _fact   = h[_dimID]; 
        _k_fact = C2PI/(4.0*L[_dimID]);
    }
    else if( mytype == PERPER ){
        _type   = PERPER;
        _fact   = 1.0; // no convolution so no multiplication by h
        _k_fact = C2PI/(L[_dimID]);
    }
    else if( mytype == UNBUNB ){
        _type   = UNBUNB;
        _fact   = h[_dimID];
        _k_fact = C2PI/(2.0*L[_dimID]);
    }
}
FFTW_plan_dim::~FFTW_plan_dim(){
    BEGIN_FUNC
    if(_plan != NULL) fftw_destroy_plan(_plan);
}
void FFTW_plan_dim::init(const size_t size[DIM],const bool isComplex){
    BEGIN_FUNC
    // !! has to be called IN THE CORRECT ORDER GIVEN BY PRIORITIES
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(size[_dimID] >= 0);

    //-------------------------------------------------------------------------
    // Get important informations
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
void FFTW_plan_dim::_init_real2real(const size_t size[DIM],const bool isComplex)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(isComplex == false);
    
    //-------------------------------------------------------------------------
    // Get the memory details (size in/out, _fieldstart, isComplex, howmany)
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _n_in  = size[_dimID];
        _n_out = size[_dimID];
    }
    else if( _isGreen){
        _n_in  = 1;
        _n_out = 1;
        INFOLOG("no need for DST/DCT for the Green's function\n");
    }

    _fieldstart  = 0; // no padding
    _howmany   = 1;
    _switch2Complex = false;
    for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
    for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];

    //-------------------------------------------------------------------------
    // update scaling factor
    //-------------------------------------------------------------------------
    _fact *= 1.0 /(2.0*size[_dimID]);

    //-------------------------------------------------------------------------
    // Get the type of Fourier transforms and imult
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

        _imult = true; // we DO have to multiply by i=sqrt(-1)

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
        ERROR("unable to init the solver required\n")
    }
}
void FFTW_plan_dim::_init_mixpoisson(const size_t size[DIM],const bool isComplex){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(isComplex == false);

    //-------------------------------------------------------------------------
    // Get the memory details (size in/out, _fieldstart, isComplex, howmany)
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

    _howmany = 1;
    for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
    for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];

    //-------------------------------------------------------------------------
    // update scaling factor
    //-------------------------------------------------------------------------
    _fact *= 1.0/(4.0*size[_dimID]);
    
    //-------------------------------------------------------------------------
    // Get the type of Fourier transforms
    //-------------------------------------------------------------------------
    if((_bc[0] == EVEN && _bc[1] == UNB)||(_bc[0] == UNB && _bc[1] == EVEN)){ // We have a DCT - we are EVEN / EVEN
        _imult = false; // we do NOT have to multiply by i=sqrt(-1)
        if(_isGreen){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT00; // DCT type I
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT00;
        }
        else{
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT10; // DCT type II
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT01;
        }
    }
    else if(_bc[0] == UNB && _bc[1] == ODD){ // We have a DCT - we are EVEN / ODD
        _imult = false; // we do NOT have to multiply by i=sqrt(-1)
        if(_isGreen){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT01; // DCT type III
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT10;
        }
        else{
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT11; // DCT type IV
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT11;
        }   
    }
    else if(_bc[0] == ODD && _bc[1] == UNB){ // we have a DST - we are ODD - EVEN
        _imult = true; // we DO have to multiply by i=sqrt(-1)
        if(_isGreen){
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT01; // DST type III
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT10;
        }
        else{
            if(_sign == FFTW_FORWARD ) _kind = FFTW_REDFT11; // DST type IV
            if(_sign == FFTW_BACKWARD) _kind = FFTW_REDFT11;
        }
    }
    else{
        ERROR("unable to init the solver required\n")     
    }
}
void FFTW_plan_dim::_init_periodic(const size_t size[DIM],const bool isComplex){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // Get the size and multiplication factor
    //-------------------------------------------------------------------------
    if(_isGreen){
        _n_in  = 1;
        _n_out = 1;
        INFOLOG("no need for DST/DCT for the Green's function\n");
    }
    else if(isComplex){
        _n_in  = size[_dimID]; // takes n complex, return n complex
        _n_out = size[_dimID];

        _switch2Complex = false;
    }
    else{
        _n_in  = size[_dimID]; // takes n real
        _n_out = _n_in/2+1; // return n_in/2 + 1 complex

        _switch2Complex = true;
    }

    _fieldstart  = 0;
    _howmany   = 1;
    for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
    for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];

    //-------------------------------------------------------------------------
    // udpate scaling factor
    //-------------------------------------------------------------------------
    _fact *= 1.0/(size[_dimID]);
}
void FFTW_plan_dim::_init_unbounded(const size_t size[DIM],const bool isComplex){
    BEGIN_FUNC

    printf(">> incomming size = %ld %ld\n",size[0],size[1]);
    //-------------------------------------------------------------------------
    // Get the size and multiplication factor
    //-------------------------------------------------------------------------
    if(isComplex){
        _n_in  = 2*size[_dimID]; // takes 2n complex, return 2n complex
        _n_out = 2*size[_dimID];

        _switch2Complex = false;
    }
    else{
        _n_in  = 2*size[_dimID]; // takes 2n real
        _n_out = _n_in/2+1; // return n_in/2 + 1 complex
        
        _switch2Complex = true;
    }

    _fieldstart  = 0;
    _howmany   = 1;
    for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
    for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];

    //-------------------------------------------------------------------------
    // udpate scaling factor
    //-------------------------------------------------------------------------
    _fact *= 1.0/(2*size[_dimID]);
}

void FFTW_plan_dim::allocate_plan(const size_t size_plan[DIM],const bool isComplex, void* data){
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // remember if the data allocated are complex or not
    //-------------------------------------------------------------------------
    _isDataComplex = isComplex;
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
void FFTW_plan_dim::_allocate_plan_real(const size_t size_ordered[DIM],void* data){
    BEGIN_FUNC

    assert(data != NULL);
    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    int rank = 1;

    // set the sizemult array usefull for strides and dists
    int sizemult[DIM]; sizemult[0] = 1;
    for(int id=1; id<DIM; id++) sizemult[id]=sizemult[id-1]*size_ordered[id-1];
    if(_isDataComplex) for(int id=1; id<DIM; id++) sizemult[id]*= 2;

    int istride = sizemult[(_orderID  )%DIM];
    int idist   = sizemult[(_orderID+1)%DIM];
    int ostride = sizemult[(_orderID  )%DIM];
    int odist   = sizemult[(_orderID+1)%DIM];
    
    // set the plan
    _plan = fftw_plan_many_r2r(  rank,(int*) (&_n_in),_howmany,
                                    (double*) data,NULL,istride,idist,
                                    (double*) data,NULL,ostride,odist,&_kind,FFTW_FLAG);

    INFOLOG ("------------------------------------------\n");
    if      (_type == R2R   ) {INFO2("## R2R plan created for plan r2r (=%d)\n",_type);}
    else if (_type == MIX   ) {INFO2("## R2R plan created for plan mix (=%d)\n",_type);}
    INFOLOG2("orderedID = %d\n",_orderID);
    INFOLOG2("howmany   = %d\n",_howmany);
    INFOLOG2("size n    = %ld\n",_n_in);
    INFOLOG3("istride (double) = %d - idist = %d\n",istride,idist);
    INFOLOG3("ostride (double) = %d - odist = %d\n",ostride,odist);
    INFOLOG ("------------------------------------------\n");
}
void FFTW_plan_dim::_allocate_plan_complex(const size_t size_ordered[DIM],void* data){
    BEGIN_FUNC

    assert(data != NULL);
    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    int rank = 1;

    // set the sizemult array usefull for strides and dists
    int sizemult[DIM]; sizemult[0] = 1;
    for(int id=1; id<DIM; id++) sizemult[id]=sizemult[id-1]*size_ordered[id-1];
    

    // ostrid = complex
    int ostride = sizemult[ _orderID   %DIM];
    int odist   = sizemult[(_orderID+1)%DIM];
    
    // incomming arrays depends if we are a complex switcher or not
    if(_switch2Complex){
        for(int id=1; id<DIM; id++) if(_isDataComplex) sizemult[id]*= 2;
        int istride = sizemult[ _orderID   %DIM];
        int idist   = sizemult[(_orderID+1)%DIM];

        INFOLOG ("------------------------------------------\n");
        if      (_type == PERPER) {INFO2("## R2C plan created for plan periodic-periodic (=%d)\n",_type);}
        else if (_type == UNBUNB) {INFO2("## R2C plan created for plan unbounded (=%d)\n",_type);}
        INFOLOG2("orderedID = %d\n",_orderID);
        INFOLOG2("howmany   = %d\n",_howmany);
        INFOLOG2("size n    = %ld\n",_n_in);
        INFOLOG3("istride (double)  = %d - idist = %d\n",istride,idist);
        INFOLOG3("ostride (complex) = %d - odist = %d\n",ostride,odist);
        INFOLOG ("------------------------------------------\n");

        // set the plan
        if(_sign == FFTW_FORWARD){
            _plan = fftw_plan_many_dft_r2c( rank,(int*) (&_n_in),_howmany,
                                            (double*) data,NULL,istride,idist,
                                            (fftw_complex*) data,NULL,ostride,odist,FFTW_FLAG);
        }
        else{
            _plan = fftw_plan_many_dft_c2r( rank,(int*) (&_n_in),_howmany,
                                            (fftw_complex*) data,NULL,ostride,odist,
                                            (double*) data,NULL,istride,idist,FFTW_FLAG);
        }

        
    }
    else{
        int istride = sizemult[ _orderID   %DIM];
        int idist   = sizemult[(_orderID+1)%DIM];
        // set the plan
        _plan = fftw_plan_many_dft( rank,(int*) (&_n_in),_howmany,
                                    (fftw_complex*) data,NULL,istride,idist,
                                    (fftw_complex*) data,NULL,ostride,odist,_sign,FFTW_FLAG);

        INFOLOG ("------------------------------------------\n");
        if      (_type == PERPER) {INFO2("## C2C plan created for plan periodic-periodic (=%d)\n",_type);}
        else if (_type == UNBUNB) {INFO2("## C2C plan created for plan unbounded (=%d)\n",_type);}
        INFOLOG2("orderedID = %d\n",_orderID);
        INFOLOG2("howmany = %d\n",_howmany);
        INFOLOG2("size n = %ld\n",_n_in);
        INFOLOG3("istride (complex) = %d - idist = %d\n",istride,idist);
        INFOLOG3("ostride (complex) = %d - odist = %d\n",ostride,odist);
        INFOLOG ("------------------------------------------\n");
    }
}


bool FFTW_plan_dim::get_isComplex() const {
    BEGIN_FUNC
    return _switch2Complex;
}
int FFTW_plan_dim::get_type() const {
    BEGIN_FUNC
    return _type;
}
void FFTW_plan_dim::get_outsize (size_t size[DIM]) const {
    BEGIN_FUNC
    size[_dimID] = _n_out;
}
void FFTW_plan_dim::get_fieldstart(size_t start[DIM]) const {
    BEGIN_FUNC
    start[_dimID] = _fieldstart;
}
void FFTW_plan_dim::get_isComplex(bool* isComplex) const {
    BEGIN_FUNC
    (*isComplex) = (*isComplex) || _switch2Complex;
}
void FFTW_plan_dim::get_outsize (const int id, size_t size[DIM]) const {
    BEGIN_FUNC
    size[id] = _n_out;
}
void FFTW_plan_dim::get_fieldstart(const int id, size_t start[DIM]) const {
    BEGIN_FUNC
    start[id] = _fieldstart;
}
void FFTW_plan_dim::get_dimID   (const int id, int dimID[DIM]) const {
    BEGIN_FUNC
    dimID[id] = _dimID;
}
void FFTW_plan_dim::set_order(const int id){
    BEGIN_FUNC
    // if the plan is called in ith position
    // its dimension is located in DIM-id-1th index
    _orderID = id;
}

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
        if(_kind == FFTW_RODFT00) {INFO("- kind = REDFT00 = DST type I\n")}
        if(_kind == FFTW_RODFT10) {INFO("- kind = REDFT10 = DST type II\n")}
        if(_kind == FFTW_RODFT01) {INFO("- kind = REDFT01 = DST type III\n")}
        if(_kind == FFTW_RODFT11) {INFO("- kind = REDFT11 = DST type IV\n")}
    }
    INFO2("- is Green   ? %d\n",_isGreen);
    INFO2("- s2Complex  ? %d\n",_switch2Complex);
    INFO2("- n_in       = %ld\n",_n_in);
    INFO2("- n_out      = %ld\n",_n_out);
    INFO2("- howmany    = %d\n",_howmany);
    INFO2("- fieldstart = %ld\n",_fieldstart);
    if      (_sign == FFTW_FORWARD ) {INFO("- FORWARD plan\n");}
    else if (_sign == FFTW_BACKWARD) {INFO("- BACKWARD plan\n");}
    
    INFO ("------------------------------------------\n");
}