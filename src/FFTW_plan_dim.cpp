
#include "FFTW_plan_dim.hpp"


FFTW_plan_dim::FFTW_plan_dim(const int dimID,const double h[DIM],const double L[DIM],const BoundaryType mybc[2],const int sign, const bool isGreen):
_dimID(dimID),
_sign(sign),
_isGreen(isGreen)
{
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
    if(_plan != NULL) fftw_destroy_plan(*_plan);
}
void FFTW_plan_dim::init(const size_t size[DIM],const bool isComplex){
    // !! has to be called IN THE CORRECT ORDER GIVEN BY PRIORITIES
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    #ifndef NDEBUG
    assert(size[_dimID] >= 0);
    #endif
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
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(isComplex == false);
    
    //-------------------------------------------------------------------------
    // Get the memory details (size in/out, padstart, isComplex, howmany)
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _n_in  = size[_dimID];
        _n_out = size[_dimID];
    }
    else if( _isGreen){
        _n_in  = 1;
        _n_out = 1;
    }

    _isComplex = false;
    _padstart  = 0; // no padding
    _howmany   = 1;
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
        INFO("no need for DST/DCT for the Green's function\n");
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
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(isComplex == false);

    //-------------------------------------------------------------------------
    // Get the memory details (size in/out, padstart, isComplex, howmany)
    //-------------------------------------------------------------------------
    if(!_isGreen){
        _n_in  =  2*size[_dimID];
        _n_out =  2*size[_dimID];
    }
    else if( _isGreen){
        _n_in  =  2*size[_dimID]+1;
        _n_out =  2*size[_dimID]+1;
    }

    _isComplex = false;

    if (_isGreen) _padstart = 0;
    else if(_bc[0] == UNB) _padstart = size[_dimID];   // padding to the left
    else if(_bc[1] == UNB) _padstart = 0;              // padding to the right

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
    //-------------------------------------------------------------------------
    // Get the size and multiplication factor
    //-------------------------------------------------------------------------
    if(_isGreen){
        _n_in  = 1;
        _n_out = 1;
    }
    else if(isComplex){
        _n_in  = size[_dimID]; // takes n complex, return n complex
        _n_out = size[_dimID];
    }
    else{
        _n_in  = size[_dimID]; // takes n real
        _n_out = _n_in/2+1; // return n_in/2 + 1 complex
    }

    _isComplex = true;
    _padstart  = 0;
    _howmany   = 1;
    for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
    for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];

    //-------------------------------------------------------------------------
    // udpate scaling factor
    //-------------------------------------------------------------------------
    _fact *= 1.0/(size[_dimID]);
}
void FFTW_plan_dim::_init_unbounded(const size_t size[DIM],const bool isComplex){
    //-------------------------------------------------------------------------
    // Get the size and multiplication factor
    //-------------------------------------------------------------------------
    if(isComplex){
        _n_in  = 2*size[_dimID]; // takes 2n complex, return 2n complex
        _n_out = 2*size[_dimID];
    }
    else{
        _n_in  = 2*size[_dimID]; // takes 2n real
        _n_out = _n_in/2+1; // return n_in/2 + 1 complex
    }

    _isComplex = true; // we keep going into complex
    _padstart  = 0;
    _howmany   = 1;
    for(int id=0       ; id<_dimID; id++) _howmany *= size[id];
    for(int id=_dimID+1; id<DIM   ; id++) _howmany *= size[id];

    //-------------------------------------------------------------------------
    // udpate scaling factor
    //-------------------------------------------------------------------------
    _fact *= 1.0/(2*size[_dimID]);
}

void FFTW_plan_dim::allocate_plan(const size_t size_plan[DIM],const bool isComplex, void* data){
    //-------------------------------------------------------------------------
    // remember if the data allocated are complex or not
    //-------------------------------------------------------------------------
    _isDataComplex = isComplex;
    //-------------------------------------------------------------------------
    // allocate the plan
    //-------------------------------------------------------------------------
    if( _type == R2R ){
        _allocate_plan_real2real(size_plan,data);
    }
    else if( _type == MIX ){
        // _n = 2*size[dimID]; // we have to double the size
        // _allocate_plan_mixpoisson(size_plan,data);
    }
    else if( _type == PERPER ){
        // _allocate_plan_periodic(size_plan,data);
    }
    else if( _type == UNBUNB ){
        // _allocate_plan_unbounded(size_plan,data);
    }
}
void FFTW_plan_dim::_allocate_plan_real2real(const size_t size_ordered[DIM],void* data){

    assert(data != NULL);
    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    INFOLOG("computation of stides\n");
    printf("data complex ? %d\n",_isDataComplex);
    int rank = 1;
    // IN - have to be in complex or double
    int sizemult[DIM];
    sizemult[0] = 1;
    for(int id=1; id<DIM; id++) sizemult[id]=sizemult[id-1]*size_ordered[id-1];
    for(int id=1; id<DIM; id++) if(_isDataComplex) sizemult[id]*= 2;

    printf("size mult = ");
    for(int id=0; id<DIM; id++) printf(" %d ",sizemult[id]);
    printf("\n");

    // if the field is complex, we have to double the sizes since we speek double here!!
    int ifact = 1;
    int ofact = 1;
    if(_isDataComplex){
        ifact = 2;
        ofact = 2;
    }
    
    printf("istride @ %d\n",_orderID   %(DIM));
    int istride = sizemult[ _orderID   %DIM];
    int idist   = sizemult[(_orderID+1)%DIM];

    int ostride = sizemult[ _orderID   %DIM];
    int odist   = sizemult[(_orderID+1)%DIM];
    
    printf("howmany = %d\n",_howmany);
    printf("size n = %ld\n",_n_in);
    printf("istrid = %d - idist = %d\n",istride,idist);
    printf("ostrid = %d - odist = %d\n",ostride,odist);
    // set the plan
    (*_plan) = fftw_plan_many_r2r(rank,(int*) &_n_in,_howmany,
                                    (double*) data,NULL,istride,idist,
                                    (double*) data,NULL,ostride,odist,&_kind,FFTW_FLAG);

    // fftw_print_plan(*_plan);

}
// void FFTW_plan_dim::_allocate_plan_mixpoisson(const size_t size_res[DIM],double* data){
//     // the array has to be (n[3] x n[2] x n[1])
//     // the jth element of transform k is at k*idist+j*istride
//     int rank = 1;
//     // IN
//     int istride = 1;
//     for(int id=0; id<_dimID; id++) istride *= size_res[id];
//     int idist   = 1;
//     for(int id=_dimID+1; id<DIM; id++) idist *= size_res[id];
//     // OUT
//     int ostride = 1;
//     for(int id=0; id<_dimID; id++) istride *= size_res[id];
//     int odist   = 1;
//     for(int id=_dimID+1; id<DIM; id++) idist *= size_res[id];
    
//     (*_plan) = fftw_plan_many_r2r(rank,(int*) &_n_plan,_howmany,data,NULL,istride,idist,data,NULL,ostride,odist,&_kind,FFTW_FLAG);

// }

bool FFTW_plan_dim::get_isComplex() const {
    return _isComplex;
}
int FFTW_plan_dim::get_type() const {
    return _type;
}
void FFTW_plan_dim::get_outsize (size_t size[DIM]) const {
    size[_dimID] = _n_out;
}
void FFTW_plan_dim::get_padstart(size_t start[DIM]) const {
    start[_dimID] = _padstart;
}
void FFTW_plan_dim::get_isComplex(bool* isComplex) const {
    (*isComplex) = (*isComplex) || _isComplex;
}
void FFTW_plan_dim::get_outsize (const int id, size_t size[DIM]) const {
    size[id] = _n_out;
}
void FFTW_plan_dim::get_padstart(const int id, size_t start[DIM]) const {
    start[id] = _padstart;
}
void FFTW_plan_dim::get_dimID   (const int id, int dimID[DIM]) const {
    dimID[id] = _dimID;
}
void FFTW_plan_dim::set_order(const int id){
    _orderID = id;
}

void FFTW_plan_dim::disp(){
    INFO ("------------------------------------------\n");
    INFO3("Plan num %d created for dimension %d\n",_orderID,_dimID);
    if      (_type == R2R   ) {INFO2("- type = real2real (=%d)\n",_type);}
    else if (_type == MIX   ) {INFO2("- type = mix (=%d)\n",_type);}
    else if (_type == PERPER) {INFO2("- type = periodic-periodic (=%d)\n",_type);}
    else if (_type == UNBUNB) {INFO2("- type = unbounded (=%d)\n",_type);}
    if      (_bc[0] == EVEN) {INFO("- { EVEN ,");}
    else if (_bc[0] == ODD)  {INFO("- { ODD  ,");}
    else if (_bc[0] == UNB)  {INFO("- { UNB  ,");}
    else if (_bc[0] == PER)  {INFO("- { PER  ,");}
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
    INFO2("- is Green  ? %d\n",_isGreen);
    INFO2("- isComplex ? %d\n",_isComplex);
    INFO2("- n_in      = %ld\n",_n_in);
    INFO2("- n_out     = %ld\n",_n_out);
    INFO2("- howmany   = %d\n",_howmany);
    INFO2("- padstart  = %ld\n",_padstart);
    if      (_sign == FFTW_FORWARD ) {INFO("- FORWARD plan\n");}
    else if (_sign == FFTW_BACKWARD) {INFO("- BACKWARD plan\n");}
    
    INFO ("------------------------------------------\n");
}