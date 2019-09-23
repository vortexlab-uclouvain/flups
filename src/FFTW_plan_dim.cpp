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

using namespace FLUPS;

/**
 * @brief Construct a new FFTW_plan_dim object
 * 
 * @param dimID the dimension id in the non-transpose reference = the field reference
 * @param h the grid spacing
 * @param L the lenght of the computational domain
 * @param mybc the boundary condition to use for this plan
 * @param sign the sign of the plan (FLUPS_FORWARD or FLUPS_BACKWARD)
 * @param isGreen boolean to indicate if the plan is intended for Green's function
 */
FFTW_plan_dim::FFTW_plan_dim(const int dimID, const double h[3], const double L[3], const BoundaryType mybc[2], const int sign, const bool isGreen) : _dimID(dimID),
                                                                                                                                                          _sign(sign),
                                                                                                                                                          _isGreen(isGreen) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    FLUPS_CHECK(dimID >= 0 && dimID < 3,"we are only creating plans on dim from 0 to 2", LOCATION);

    //-------------------------------------------------------------------------
    // Initialisation of the sizes and types
    //-------------------------------------------------------------------------
    _bc[0] = mybc[0];
    _bc[1] = mybc[1];
    //determine the type of the solver
    int mytype = _bc[0] + _bc[1];
    //-------------------------------------------------------------------------
    // Get type and mult factors
    //-------------------------------------------------------------------------
    if (mytype <= SYMSYM) {
        _type     = SYMSYM;
        _normfact = 1.0;
        _volfact  = 1.0;  // no convolution so no multiplication by h
        _kfact    = c_2pi / (2.0 * L[_dimID]);
        _koffset  = 0.0;
        if (_isGreen) _isSpectral = true;
    } else if (mytype <= MIXUNB) {
        _type     = MIXUNB;
        _normfact = 1.0;
        _volfact  = h[_dimID];
        _kfact    = c_2pi / (4.0 * L[_dimID]);
        _koffset  = 0.0;
    } else if (mytype == PERPER) {
        _type     = PERPER;
        _normfact = 1.0;
        _volfact  = 1.0;  // no convolution so no multiplication by h
        _kfact    = c_2pi / (L[_dimID]);
        _koffset  = 0.0;
        if (_isGreen) _isSpectral = true;
    } else if (mytype == UNBUNB) {
        _type     = UNBUNB;
        _normfact = 1.0;
        _volfact  = h[_dimID];
        _kfact    = c_2pi / (2.0 * L[_dimID]);
        _koffset  = 0.0;
    } else {
        FLUPS_ERROR("Invalid combination of BCs", LOCATION);
    }
}
FFTW_plan_dim::~FFTW_plan_dim() {
    BEGIN_FUNC;
    if (_plan != NULL) fftw_destroy_plan(_plan);
}

/**
 * @brief Initialize the FFTW_plan_dim by performing a 'dry run'
 * 
 * The function redirects to one of the init functions depending on the type:
 * - _init_real2real()
 * - _init_mixunbounded()
 * - _init_periodic()
 * - _init_unbounded()
 * 
 * Each of the sub-function initializes the following variables
 * - #_n_in the size of data provided as input to the FFTW (i.e. the number of real or complex numbers)
 * - #_n_out the size of data that comes out of the FFTW
 * - #_fieldstart the index to start the FFTW (non zero for mixunbounded solvers)
 * - #_isr2c is true if this plan switches to the complex numbers
 * - #_imult is true if we used a DST
 * - #_kind the kind of FFTW plan to execute (for SYMSYM and MIXUNB plans only)
 * - #_symstart the symmetry start = id of symmetry, for the Green's function only
 * 
 * @param size the current size of data in during dry run (hence already partially transformed)
 * @param isComplex the current complex state of the data
 */
void FFTW_plan_dim::init(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(size[_dimID] >= 0);

    //-------------------------------------------------------------------------
    // redirect to the corresponding subfunction
    //-------------------------------------------------------------------------
    if (_type == SYMSYM) {
        _init_real2real(size, isComplex);
    } else if (_type == MIXUNB) {
        _init_mixunbounded(size, isComplex);
    } else if (_type == PERPER) {
        //this is the only transform that could give a R2C on data and being spectral for green
        _init_periodic(size, isComplex); 
    } else if (_type == UNBUNB) {
        _init_unbounded(size, isComplex);
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
void FFTW_plan_dim::_init_real2real(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    FLUPS_CHECK(isComplex == false,"the data cannot be complex", LOCATION);

    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__isr2c)  */
    //-------------------------------------------------------------------------
    if (!_isGreen) {
        _n_in  = size[_dimID];
        _n_out = size[_dimID];
    } else {
        // if ODD-EVEN or EVEN-ODD the FFT for the Green function would have been DCT Type III
        // hence on a size of n
        if (_bc[0] != _bc[1]) {
            _n_in  = size[_dimID];
            _n_out = size[_dimID];
        } else {
            _n_in       = size[_dimID] + 1;
            _n_out      = size[_dimID] + 1;
            _ignoreMode = true;
        }
    }

    _fieldstart = 0;

    // no switch to complex
    _isr2c      = false;

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0;  // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (2.0 * size[_dimID]);

    //-------------------------------------------------------------------------
    /** - Get the #_kind of Fourier transforms and #_imult */
    //-------------------------------------------------------------------------
    if (_isGreen) {
        // if we are doing odd-even we have to use shifted FFTW plans
        if (_bc[0] != _bc[1]) {
            _koffset = 0.5;
        }
        if (_bc[0] == ODD && _bc[1] == ODD) {
            // if we do a ODD ODD, we have to shift the Green's function
            _shiftgreen = 1;
        }
        return;
    } else if (_bc[0] == EVEN) {  // We have a DCT

        _imult = false;  // we do NOT have to multiply by i=sqrt(-1)

        if (_bc[1] == EVEN) {
            if (_sign == FLUPS_FORWARD) _kind = FFTW_REDFT10;  // DCT type II
            if (_sign == FLUPS_BACKWARD) _kind = FFTW_REDFT01; // DCT type III
        } else if (_bc[1] == ODD) {
            if (_sign == FLUPS_FORWARD) _kind = FFTW_REDFT11;  // DCT type IV
            if (_sign == FLUPS_BACKWARD) _kind = FFTW_REDFT11; // DCT type IV
            _koffset = 0.5;
        }
    } else if (_bc[0] == ODD) {  // We have a DST

        _imult = true;  // we DO have to multiply by -i=-sqrt(-1)

        if (_bc[1] == ODD) {
            if (_sign == FLUPS_FORWARD) _kind = FFTW_RODFT10;  // DST type II
            if (_sign == FLUPS_BACKWARD) _kind = FFTW_RODFT01; // DST type III
        } else if (_bc[1] == EVEN) {
            if (_sign == FLUPS_FORWARD) _kind = FFTW_RODFT11;  // DST type IV
            if (_sign == FLUPS_BACKWARD) _kind = FFTW_RODFT11; // DST type IV
            _koffset = 0.5;
        }
    } else {
        FLUPS_ERROR("unable to init the solver required", LOCATION);
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
void FFTW_plan_dim::_init_mixunbounded(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    FLUPS_CHECK(isComplex == false,"the data cannot be complex", LOCATION);

    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart and #__isr2c)  */
    //-------------------------------------------------------------------------
    if (!_isGreen) {
        _n_in  = 2 * size[_dimID];
        _n_out = 2 * size[_dimID];
    } else if (_isGreen) {
        //Different because the Green's function is to be seen as "vertex centered",
        //as opposed to data which are "cell cenetered".
        _n_in       = 2 * size[_dimID] + 1;
        _n_out      = 2 * size[_dimID] + 1;
        _ignoreMode = true;
    }

    _isr2c = false;

    if (_isGreen)
        _fieldstart = 0;
    else if (_bc[0] == UNB)
        _fieldstart = size[_dimID];  // padding to the left
    else if (_bc[1] == UNB)
        _fieldstart = 0;  // padding to the right

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0;  // if no symmetry is needed, set to 0
    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (4.0 * size[_dimID]);
    //-------------------------------------------------------------------------
    /** - Get the #_kind of Fourier transforms, #_imult and #_shiftgreen */
    //-------------------------------------------------------------------------
    if (_isGreen) {
        _imult = false;
        // set the shiftg Green to 1 if we do ODD-ODD bc
        if ((_bc[0] == UNB && _bc[1] == ODD) || (_bc[0] == ODD && _bc[1] == UNB)) {
            _shiftgreen = 1;
        }
        // The Green function is ALWAYS EVEN - EVEN
        if (_sign == FLUPS_FORWARD) _kind = FFTW_REDFT00;  // DCT type I
        if (_sign == FLUPS_BACKWARD) _kind = FFTW_REDFT00;

    } else {
        if ((_bc[0] == EVEN && _bc[1] == UNB) || (_bc[0] == UNB && _bc[1] == EVEN)) {  // We have a DCT - we are EVEN - EVEN over 2L
            _imult      = false;
            if (_sign == FLUPS_FORWARD) _kind = FFTW_REDFT10;  // DCT type II
            if (_sign == FLUPS_BACKWARD) _kind = FFTW_REDFT01; // DCT type III
        } else if ((_bc[0] == UNB && _bc[1] == ODD) || (_bc[0] == ODD && _bc[1] == UNB)) {  // We have a DST - we are ODD - ODD over 2L
            _imult      = true;
            if (_sign == FLUPS_FORWARD) _kind = FFTW_RODFT10;  // DST type II
            if (_sign == FLUPS_BACKWARD) _kind = FFTW_RODFT01; // DST type III
        } else {
            FLUPS_ERROR("unable to init the solver required", LOCATION);
        }
    }
}

/**
 * @brief Initialize for a periodic plan
 * 
 * @param size 
 * @param isComplex if the data is already complex
 * 
 * ----------------------------------------
 * We do the following operations
 */
void FFTW_plan_dim::_init_periodic(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__isr2c)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        _n_in  = size[_dimID];  // takes n complex, return n complex
        _n_out = size[_dimID];

        _isr2c = false;
    } else {
        _n_in  = size[_dimID];   // takes n real
        _n_out = _n_in / 2 + 1;  // return n_in/2 + 1 complex

        _isr2c = true;
    }
    
    _fieldstart = 0;
    _shiftgreen = 0;

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0;  // if no symmetry is needed, set to 0
    if(isComplex ){
        _symstart = size[_dimID]/2.0;
    }
    
    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (size[_dimID]);
    //-------------------------------------------------------------------------
    /** - Get the #_imult factor */
    //-------------------------------------------------------------------------
    _imult = false;
}

/**
 * @brief initialize the plan for unbounded solvers
 * 
 * @param size the size
 * @param isComplex if the data is already complex
 * 
 * 
 *--------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim::_init_unbounded(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__isr2c)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        _n_in  = 2 * size[_dimID];  // takes 2n complex, return 2n complex
        _n_out = 2 * size[_dimID];

        _isr2c = false;
    } else {
        _n_in  = 2 * size[_dimID];  // takes 2n real (because of the padding)
        _n_out = _n_in / 2 + 1;     // return n_in/2 + 1 complex

        _isr2c = true;
    }

    _fieldstart = 0;
    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = size[_dimID];
    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (2.0 * size[_dimID]);
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
 * @param size_plan the size of the data BEFORE THE PLAN is executed
 * @param isComplex if the transpoed data is complex or real
 * @param data the pointer to the transposed data (has to be allocated)
 */
void FFTW_plan_dim::allocate_plan(const Topology *topo, double* data) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // allocate the plan
    //-------------------------------------------------------------------------
    if (_type == SYMSYM || _type == MIXUNB) {
        _allocate_plan_real(topo, data);
    } else if (_type == PERPER || _type == UNBUNB) {
        _allocate_plan_complex(topo, data);
    }
}

/**
 * @brief Allocate a plan that only treats real numbers
 * 
 * @note
 * The howmany is computed using the memory information provided
 * 
 * @warning
 * If there is a r2c tranform we will perfom upto 2 spurious transform.
 * This is not important since the Green's transform is only performed once
 * 
 * @param memsize the size of the data BEFORE THE PLAN is executed
 * @param data the pointer to the transposed data (has to be allocated)
 * 
 */
void FFTW_plan_dim::_allocate_plan_real(const Topology *topo, double* data) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    assert(data != NULL);

    const int memsize[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)}; //the "current" size, corresponding to size_tmp during the dry run, see _init_plansAndTopos

    //-------------------------------------------------------------------------
    /** - If is Green and #_type is SYMSYM, exit */
    //-------------------------------------------------------------------------
    if (_isGreen && _type == SYMSYM) {
        _plan = NULL;

        FLUPS_INFO("------------------------------------------");
        FLUPS_INFO("## no real to real plan created for Green");
        FLUPS_INFO("------------------------------------------");
        return;
    }

    //-------------------------------------------------------------------------
    /** - Get compute the rank and the stides  */
    //-------------------------------------------------------------------------
    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    int rank    = 1;
    int istride = topo->nf() ;
    int ostride = topo->nf() ;
    int idist   = memsize[_dimID]*topo->nf();
    int odist   = memsize[_dimID]*topo->nf();

    //-------------------------------------------------------------------------
    /** - Compute howmany  */
    //-------------------------------------------------------------------------
    int howmany = 1;
    for (int id = 0; id < _dimID; id++) howmany *= memsize[id];
    for (int id = _dimID + 1; id < 3; id++) howmany *= memsize[id];

    //-------------------------------------------------------------------------
    /** - Create the plan  */
    //-------------------------------------------------------------------------
    _plan = fftw_plan_many_r2r(rank, (int*)(&_n_in), howmany,
                               data, NULL, istride, idist,
                               data, NULL, ostride, odist, &_kind, FFTW_FLAG);

    FLUPS_INFO("------------------------------------------");
    if (_type == SYMSYM) {
        FLUPS_INFO("## SYMSYM plan created for plan r2r (=%d)", _type);
    } else if (_type == MIXUNB) {
        FLUPS_INFO("## MIXUNB plan created for plan mix (=%d)", _type);
    }
    FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
    FLUPS_INFO("howmany   = %d", howmany);
    FLUPS_INFO("size n    = %d", _n_in);
    FLUPS_INFO("istride (double) = %d - idist = %d", istride, idist);
    FLUPS_INFO("ostride (double) = %d - odist = %d", ostride, odist);
    FLUPS_INFO("------------------------------------------");
}

/**
 * @brief allocate a plan that treats complex numbers (r2c or c2c)
 * 
 * @note
 * The howmany is computed using the memory information provided
 * 
 * @warning
 * If there is a r2c tranform we will perfom upto 2 spurious transform.
 * This is not important since the Green's transform is only performed once
 * 
 * 
 * @param memsize the size of the data BEFORE THE PLAN is executed
 * @param data memory
 */
void FFTW_plan_dim::_allocate_plan_complex(const Topology *topo, double* data) {
    BEGIN_FUNC;

    assert(data != NULL);

    const int memsize[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)}; //the "current" size, corresponding to size_tmp during the dry run, see _init_plansAndTopos

    if (_isGreen && _type == PERPER) {
        _plan = NULL;

        FLUPS_INFO("------------------------------------------");
        FLUPS_INFO("## no DFT plan created for Green");
        FLUPS_INFO("------------------------------------------");
        return;
    }

    // the jth element of transform k is at k*idist+j*istride
    int rank = 1;
    // Compute howmany
    int howmany = 1;
    for (int id = 0; id < _dimID; id++) howmany *= memsize[id];
    for (int id = _dimID + 1; id < 3; id++) howmany *= memsize[id];

    // strides
    int istride = 1;
    int ostride = 1;
    int idist   = memsize[_dimID];
    int odist   = memsize[_dimID];

    // incomming arrays depends if we are a complex switcher or not
    if (_isr2c) {
        // idist has been obtained from a the incomming size
        odist /= 2;

        FLUPS_INFO("------------------------------------------");
        if (_type == PERPER) {
            FLUPS_INFO("## R2C plan created for plan periodic-periodic (=%d)", _type);
        } else if (_type == UNBUNB) {
            FLUPS_INFO("## R2C plan created for plan unbounded (=%d)", _type);
        }
        // FLUPS_INFO("orderedID = %d",_orderID);
        if (_sign == FLUPS_FORWARD) {
            FLUPS_INFO("FORWARD transfrom");
        } else if (_sign == FLUPS_BACKWARD) {
            FLUPS_INFO("BACKWARD transfrom");
        }
        FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
        FLUPS_INFO("dimID     = %d", _dimID);
        FLUPS_INFO("howmany   = %d", howmany);
        FLUPS_INFO("size n    = %d", _n_in);
        FLUPS_INFO("istride (double)  = %d - idist = %d", istride, idist);
        FLUPS_INFO("ostride (complex) = %d - odist = %d", ostride, odist);
        FLUPS_INFO("------------------------------------------");

        // set the plan - there is no offset in r2c or c2c possible
        if (_sign == FLUPS_FORWARD) {
            _plan = fftw_plan_many_dft_r2c(rank, (int*)(&_n_in), howmany,
                                           data, NULL, istride, idist,
                                           (fftw_complex*)data, NULL, ostride, odist, FFTW_FLAG);
        } else {
            _plan = fftw_plan_many_dft_c2r(rank, (int*)(&_n_in), howmany,
                                           (fftw_complex*)data, NULL, ostride, odist,
                                           data, NULL, istride, idist, FFTW_FLAG);
        }

    } else {
        // set the plan
        _plan = fftw_plan_many_dft(rank, (int*)(&_n_in), howmany,
                                   (fftw_complex*)data, NULL, istride, idist,
                                   (fftw_complex*)data, NULL, ostride, odist, _sign, FFTW_FLAG);

    }
}

/**
 * @brief Executes the plan
 * 
 */
void FFTW_plan_dim::execute_plan() {
    BEGIN_FUNC;

    FLUPS_CHECK(!_isSpectral,"Trying to execute a plan for data which is already spectral", LOCATION);

    // run the plan
    if (_type == SYMSYM) {
        FLUPS_INFO(">> Doing plan real2real for dim %d", _dimID);
    } else if (_type == MIXUNB) {
        FLUPS_INFO(">> Doing plan mix for dim %d", _dimID);
    } else if (_type == PERPER) {
        FLUPS_INFO(">> Doing plan periodic-periodic for dim %d", _dimID);
    } else if (_type == UNBUNB) {
        FLUPS_INFO(">> Doing plan unbounded for dim %d", _dimID);
    }
    fftw_execute(_plan);
}

/**
 * @brief display the FFTW_plan_dim object
 * 
 */
void FFTW_plan_dim::disp() {
    BEGIN_FUNC;
    FLUPS_INFO("------------------------------------------");
    FLUPS_INFO("## Plan num created for dimension %d", _dimID);
    if (_type == SYMSYM) {
        FLUPS_INFO("- type = real2real (=%d)", _type);
    } else if (_type == MIXUNB) {
        FLUPS_INFO("- type = mix (=%d)", _type);
    } else if (_type == PERPER) {
        FLUPS_INFO("- type = periodic-periodic (=%d)", _type);
    } else if (_type == UNBUNB) {
        FLUPS_INFO("- type = unbounded (=%d)", _type);
    }
    if (_bc[0] == EVEN) {
        FLUPS_INFO("- bc = { EVEN ,");
    } else if (_bc[0] == ODD) {
        FLUPS_INFO("- bc = { ODD  ,");
    } else if (_bc[0] == UNB) {
        FLUPS_INFO("- bc = { UNB  ,");
    } else if (_bc[0] == PER) {
        FLUPS_INFO("- bc = { PER  ,");
    }
    if (_bc[1] == EVEN) {
        FLUPS_INFO(" EVEN}");
    } else if (_bc[1] == ODD) {
        FLUPS_INFO(" ODD}");
    } else if (_bc[1] == UNB) {
        FLUPS_INFO(" UNB}");
    } else if (_bc[1] == PER) {
        FLUPS_INFO(" PER}");
    }
    if ((_type == SYMSYM && !_isGreen) || _type == MIXUNB) {
        if (_kind == FFTW_REDFT00) {
            FLUPS_INFO("- kind = REDFT00 = DCT type I");
        }
        if (_kind == FFTW_REDFT10) {
            FLUPS_INFO("- kind = REDFT10 = DCT type II");
        }
        if (_kind == FFTW_REDFT01) {
            FLUPS_INFO("- kind = REDFT01 = DCT type III");
        }
        if (_kind == FFTW_REDFT11) {
            FLUPS_INFO("- kind = REDFT11 = DCT type IV");
        }
        if (_kind == FFTW_RODFT00) {
            FLUPS_INFO("- kind = RODFT00 = DST type I");
        }
        if (_kind == FFTW_RODFT10) {
            FLUPS_INFO("- kind = RODFT10 = DST type II");
        }
        if (_kind == FFTW_RODFT01) {
            FLUPS_INFO("- kind = RODFT01 = DST type III");
        }
        if (_kind == FFTW_RODFT11) {
            FLUPS_INFO("- kind = RODFT11 = DST type IV");
        }
    }
    FLUPS_INFO("- dimID      = %d", _dimID);
    FLUPS_INFO("- is Green   ? %d", _isGreen);
    FLUPS_INFO("- s2Complex  ? %d", _isr2c);
    FLUPS_INFO("- n_in       = %d", _n_in);
    FLUPS_INFO("- n_out      = %d", _n_out);
    FLUPS_INFO("- fieldstart = %d", _fieldstart);
    FLUPS_INFO("- shiftgreen = %d", _shiftgreen);
    FLUPS_INFO("- isSpectral ? %d", _isSpectral);
    if (_sign == FLUPS_FORWARD) {
        FLUPS_INFO("- FORWARD plan");
    } else if (_sign == FLUPS_BACKWARD) {
        FLUPS_INFO("- BACKWARD plan");
    }

    FLUPS_INFO("------------------------------------------");
}