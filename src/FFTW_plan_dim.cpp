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
 * @param sign the sign of the plan (UP_FORWARD or UP_BACKWARD)
 * @param isGreen boolean to indicate if the plan is intended for Green's function
 */
FFTW_plan_dim::FFTW_plan_dim(const int dimID, const double h[DIM], const double L[DIM], const BoundaryType mybc[2], const int sign, const bool isGreen) : _dimID(dimID),
                                                                                                                                                          _sign(sign),
                                                                                                                                                          _isGreen(isGreen) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(dimID < DIM);
    assert(dimID >= 0);

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
        if (_isGreen) _isSpectral = true;
    } else if (mytype <= MIXUNB) {
        _type     = MIXUNB;
        _normfact = 1.0;
        _volfact  = h[_dimID];
        _kfact    = c_2pi / (4.0 * L[_dimID]);
    } else if (mytype == PERPER) {
        _type     = PERPER;
        _normfact = 1.0;
        _volfact  = 1.0;  // no convolution so no multiplication by h
        _kfact    = c_2pi / (L[_dimID]);
        if (_isGreen) _isSpectral = true;
    } else if (mytype == UNBUNB) {
        _type     = UNBUNB;
        _normfact = 1.0;
        _volfact  = h[_dimID];
        _kfact    = c_2pi / (2.0 * L[_dimID]);
    } else {
        UP_ERROR("Invalid combination of BCs")
    }
}
FFTW_plan_dim::~FFTW_plan_dim() {
    BEGIN_FUNC
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
 * - #_symstart the symmetry start, for the Green's function only
 * 
 * @param size the current size of data in during dry run (hence already partially transformed)
 * @param isComplex the current complex state of the data
 */
void FFTW_plan_dim::init(const int size[DIM], const bool isComplex) {
    BEGIN_FUNC
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
        // _n = 2*size[dimID]; // we have to double the size
        _init_mixunbounded(size, isComplex);
    } else if (_type == PERPER) {
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
void FFTW_plan_dim::_init_real2real(const int size[DIM], const bool isComplex) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0(isComplex == false,"the data cannot be complex");

    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__isr2c)  */
    //-------------------------------------------------------------------------
    _n_in       = size[_dimID];
    _n_out      = size[_dimID];
    _fieldstart = 0;

    // no switch to complex
    _isr2c      = false;
    _shiftgreen = 0;

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
        return;
    } else if (_bc[0] == EVEN) {  // We have a DCT

        _imult = false;  // we do NOT have to multiply by i=sqrt(-1)

        if (_bc[1] == EVEN) {
            if (_sign == UP_FORWARD) _kind = FFTW_REDFT10;  // DCT type II
            if (_sign == UP_BACKWARD) _kind = FFTW_REDFT01;
        } else if (_bc[1] == ODD) {
            if (_sign == UP_FORWARD) _kind = FFTW_REDFT01;  // DCT type III
            if (_sign == UP_BACKWARD) _kind = FFTW_REDFT10;
        }
    } else if (_bc[0] == ODD) {  // We have a DST

        _imult = true;  // we DO have to multiply by -i=-sqrt(-1)

        if (_bc[1] == ODD) {
            if (_sign == UP_FORWARD) _kind = FFTW_RODFT10;  // DST type II
            if (_sign == UP_BACKWARD) _kind = FFTW_RODFT01;
        } else if (_bc[1] == EVEN) {
            if (_sign == UP_FORWARD) _kind = FFTW_REDFT11;  // DST type IV
            if (_sign == UP_BACKWARD) _kind = FFTW_REDFT11;
        }
    } else {
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
void FFTW_plan_dim::_init_mixunbounded(const int size[DIM], const bool isComplex) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0(isComplex == false,"the data cannot be complex");

    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart and #__isr2c)  */
    //-------------------------------------------------------------------------
    if (!_isGreen) {
        _n_in  = 2 * size[_dimID];
        _n_out = 2 * size[_dimID];
    } else if (_isGreen) {
        //Different because the Green's function is to be seen as "vertex centered",
        //as opposed to data which are "cell cenetered".
        _n_in  = 2 * size[_dimID] + 1;
        _n_out = 2 * size[_dimID] + 1;
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
        _imult      = false;
        _shiftgreen = 0;
        // The Green function is ALWAYS EVEN - EVEN
        if (_sign == UP_FORWARD) _kind = FFTW_REDFT00;  // DCT type I
        if (_sign == UP_BACKWARD) _kind = FFTW_REDFT00;
    } else {
        if ((_bc[0] == EVEN && _bc[1] == UNB) || (_bc[0] == UNB && _bc[1] == EVEN)) {  // We have a DCT - we are EVEN - EVEN
            _imult      = false;
            _shiftgreen = 0;
            if (_sign == UP_FORWARD) _kind = FFTW_REDFT10;  // DCT type II
            if (_sign == UP_BACKWARD) _kind = FFTW_REDFT01;
        } else if ((_bc[0] == UNB && _bc[1] == ODD) || (_bc[0] == ODD && _bc[1] == UNB)) {  // We have a DCT - we are EVEN - ODD
            _imult      = true;
            _shiftgreen = 1;
            if (_sign == UP_FORWARD) _kind = FFTW_RODFT10;  // DST type II
            if (_sign == UP_BACKWARD) _kind = FFTW_RODFT01;
        } else {
            UP_ERROR("unable to init the solver required\n")
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
void FFTW_plan_dim::_init_periodic(const int size[DIM], const bool isComplex) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - get the memory details (#_n_in, #_n_out, #_fieldstart, #_shiftgreen and #__isr2c)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        _n_in  = 2 * size[_dimID];  // takes n complex, return n complex
        _n_out = 2 * size[_dimID];

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
void FFTW_plan_dim::_init_unbounded(const int size[DIM], const bool isComplex) {
    BEGIN_FUNC
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
    _shiftgreen = 0;
    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = size[_dimID];
    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (2 * size[_dimID]);
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
void FFTW_plan_dim::allocate_plan(const int size_plan[DIM], double* data) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    // allocate the plan
    //-------------------------------------------------------------------------
    if (_type == SYMSYM || _type == MIXUNB) {
        _allocate_plan_real(size_plan, data);
    } else if (_type == PERPER || _type == UNBUNB) {
        _allocate_plan_complex(size_plan, data);
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
 * @param memsize the size of the data BEFORE THE PLAN is executed
 * @param data the pointer to the transposed data (has to be allocated)
 * 
 */
void FFTW_plan_dim::_allocate_plan_real(const int memsize[DIM], double* data) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    assert(data != NULL);

    //-------------------------------------------------------------------------
    /** - If is Green and #_type is SYMSYM, exit */
    //-------------------------------------------------------------------------
    if (_isGreen && _type == SYMSYM) {
        _plan = NULL;

        INFOLOG("------------------------------------------\n");
        INFOLOG("## no real to real plan created for Green\n");
        INFOLOG("------------------------------------------\n");
        return;
    }

    //-------------------------------------------------------------------------
    /** - Get compute the rank and the stides  */
    //-------------------------------------------------------------------------
    // the array has to be (n[3] x n[2] x n[1])
    // the jth element of transform k is at k*idist+j*istride
    int rank    = 1;
    int istride = 1;
    int ostride = 1;
    int idist   = memsize[_dimID];
    int odist   = memsize[_dimID];

    //-------------------------------------------------------------------------
    /** - If is Green, compute #_howmany  */
    //-------------------------------------------------------------------------
    _howmany = 1;
    for (int id = 0; id < _dimID; id++) _howmany *= memsize[id];
    for (int id = _dimID + 1; id < DIM; id++) _howmany *= memsize[id];

    //-------------------------------------------------------------------------
    /** - Create the plan  */
    //-------------------------------------------------------------------------
    _plan = fftw_plan_many_r2r(rank, (int*)(&_n_in), _howmany,
                               data, NULL, istride, idist,
                               data, NULL, ostride, odist, &_kind, FFTW_FLAG);

    INFOLOG("------------------------------------------\n");
    if (_type == SYMSYM) {
        INFOLOG2("## SYMSYM plan created for plan r2r (=%d)\n", _type);
    } else if (_type == MIXUNB) {
        INFOLOG2("## SYMSYM plan created for plan mix (=%d)\n", _type);
    }
    INFOLOG4("memsize = %d x %d x %d\n", memsize[0], memsize[1], memsize[2]);
    INFOLOG2("howmany   = %d\n", _howmany);
    INFOLOG2("size n    = %d\n", _n_in);
    INFOLOG3("istride (double) = %d - idist = %d\n", istride, idist);
    INFOLOG3("ostride (double) = %d - odist = %d\n", ostride, odist);
    INFOLOG("------------------------------------------\n");
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
 * @param memsize the size of the data BEFORE THE PLAN is executed
 * @param data memory
 */
void FFTW_plan_dim::_allocate_plan_complex(const int memsize[DIM], double* data) {
    BEGIN_FUNC

    assert(data != NULL);

    if (_isGreen && _type == PERPER) {
        _plan = NULL;

        INFOLOG("------------------------------------------\n");
        INFOLOG("## no DFT plan created for Green\n");
        INFOLOG("------------------------------------------\n");
        return;
    }

    // the jth element of transform k is at k*idist+j*istride
    int rank = 1;
    // if we are green we need to recompute howmany
    _howmany = 1;
    for (int id = 0; id < _dimID; id++) _howmany *= memsize[id];
    for (int id = _dimID + 1; id < DIM; id++) _howmany *= memsize[id];

    // strides
    int istride = 1;
    int ostride = 1;
    int idist   = memsize[_dimID];
    int odist   = memsize[_dimID];

    // incomming arrays depends if we are a complex switcher or not
    if (_isr2c) {
        // idist has been obtained from a the incomming size
        odist /= 2;

        INFOLOG("------------------------------------------\n");
        if (_type == PERPER) {
            INFO2("## R2C plan created for plan periodic-periodic (=%d)\n", _type);
        } else if (_type == UNBUNB) {
            INFO2("## R2C plan created for plan unbounded (=%d)\n", _type);
        }
        // INFOLOG2("orderedID = %d\n",_orderID);
        if (_sign == UP_FORWARD) {
            INFOLOG("FORWARD transfrom\n");
        } else if (_sign == UP_BACKWARD) {
            INFOLOG("BACKWARD transfrom\n");
        }
        INFOLOG4("memsize = %d x %d x %d\n", memsize[0], memsize[1], memsize[2]);
        INFOLOG2("dimID     = %d\n", _dimID);
        INFOLOG2("howmany   = %d\n", _howmany);
        INFOLOG2("size n    = %d\n", _n_in);
        INFOLOG3("istride (double)  = %d - idist = %d\n", istride, idist);
        INFOLOG3("ostride (complex) = %d - odist = %d\n", ostride, odist);
        INFOLOG("------------------------------------------\n");

        // set the plan - there is no offset in r2c or c2c possible
        if (_sign == UP_FORWARD) {
            _plan = fftw_plan_many_dft_r2c(rank, (int*)(&_n_in), _howmany,
                                           data, NULL, istride, idist,
                                           (fftw_complex*)data, NULL, ostride, odist, FFTW_FLAG);
        } else {
            _plan = fftw_plan_many_dft_c2r(rank, (int*)(&_n_in), _howmany,
                                           (fftw_complex*)data, NULL, ostride, odist,
                                           data, NULL, istride, idist, FFTW_FLAG);
        }

    } else {
        // set the plan
        _plan = fftw_plan_many_dft(rank, (int*)(&_n_in), _howmany,
                                   (fftw_complex*)data, NULL, istride, idist,
                                   (fftw_complex*)data, NULL, ostride, odist, _sign, FFTW_FLAG);

        // INFOLOG ("------------------------------------------\n");
        // if      (_type == PERPER) {INFO2("## C2C plan created for plan periodic-periodic (=%d)\n",_type);}
        // else if (_type == UNBUNB) {INFO2("## C2C plan created for plan unbounded (=%d)\n",_type);}
        // // INFOLOG2("orderedID = %d\n",_orderID);
        // if      (_sign == UP_FORWARD){ INFOLOG("FORWARD transfrom\n");}
        // else if (_sign == UP_BACKWARD){ INFOLOG("BACKWARD transfrom\n");}
        // INFOLOG4("size = %d x %d x %d\n",memsize[0],memsize[1],memsize[2]);
        // INFOLOG2("dimID = %d\n",_dimID);
        // INFOLOG2("howmany = %d\n",_howmany);
        // INFOLOG2("size n = %d\n",_n_in);
        // INFOLOG3("istride (complex) = %d - idist = %d\n",istride,idist);
        // INFOLOG3("ostride (complex) = %d - odist = %d\n",ostride,odist);
        // INFOLOG ("------------------------------------------\n");
    }
}

/**
 * @brief Executes the plan
 * 
 */
void FFTW_plan_dim::execute_plan() {
    BEGIN_FUNC
    // run the plan
    if (_type == SYMSYM) {
        INFO2(">> Doing plan real2real for dim %d\n", _dimID);
    } else if (_type == MIXUNB) {
        INFO2(">> Doing plan mix for dim %d\n", _dimID);
    } else if (_type == PERPER) {
        INFO2(">> Doing plan periodic-periodic for dim %d\n", _dimID);
    } else if (_type == UNBUNB) {
        INFO2(">> Doing plan unbounded for dim %d\n", _dimID);
    }
    fftw_execute(_plan);
}

/**
 * @brief display the FFTW_plan_dim object
 * 
 */
void FFTW_plan_dim::disp() {
    BEGIN_FUNC
    INFO("------------------------------------------\n");
    INFO2("## Plan num created for dimension %d\n", _dimID);
    if (_type == SYMSYM) {
        INFO2("- type = real2real (=%d)\n", _type);
    } else if (_type == MIXUNB) {
        INFO2("- type = mix (=%d)\n", _type);
    } else if (_type == PERPER) {
        INFO2("- type = periodic-periodic (=%d)\n", _type);
    } else if (_type == UNBUNB) {
        INFO2("- type = unbounded (=%d)\n", _type);
    }
    if (_bc[0] == EVEN) {
        INFO("- bc = { EVEN ,");
    } else if (_bc[0] == ODD) {
        INFO("- bc = { ODD  ,");
    } else if (_bc[0] == UNB) {
        INFO("- bc = { UNB  ,");
    } else if (_bc[0] == PER) {
        INFO("- bc = { PER  ,");
    }
    if (_bc[1] == EVEN) {
        INFO(" EVEN}\n");
    } else if (_bc[1] == ODD) {
        INFO(" ODD}\n");
    } else if (_bc[1] == UNB) {
        INFO(" UNB}\n");
    } else if (_bc[1] == PER) {
        INFO(" PER}\n");
    }
    if (_type == SYMSYM || _type == MIXUNB) {
        if (_kind == FFTW_REDFT00) {
            INFO("- kind = REDFT00 = DCT type I\n")
        }
        if (_kind == FFTW_REDFT10) {
            INFO("- kind = REDFT10 = DCT type II\n")
        }
        if (_kind == FFTW_REDFT01) {
            INFO("- kind = REDFT01 = DCT type III\n")
        }
        if (_kind == FFTW_REDFT11) {
            INFO("- kind = REDFT11 = DCT type IV\n")
        }
        if (_kind == FFTW_RODFT00) {
            INFO("- kind = RODFT00 = DST type I\n")
        }
        if (_kind == FFTW_RODFT10) {
            INFO("- kind = RODFT10 = DST type II\n")
        }
        if (_kind == FFTW_RODFT01) {
            INFO("- kind = RODFT01 = DST type III\n")
        }
        if (_kind == FFTW_RODFT11) {
            INFO("- kind = RODFT11 = DST type IV\n")
        }
    }
    INFO2("- dimID      = %d\n", _dimID);
    INFO2("- is Green   ? %d\n", _isGreen);
    INFO2("- s2Complex  ? %d\n", _isr2c);
    INFO2("- n_in       = %d\n", _n_in);
    INFO2("- n_out      = %d\n", _n_out);
    INFO2("- howmany    = %d\n", _howmany);
    INFO2("- fieldstart = %d\n", _fieldstart);
    INFO2("- shiftgreen = %d\n", _shiftgreen);
    INFO2("- isSpectral ? %d\n", _isSpectral);
    if (_sign == UP_FORWARD) {
        INFO("- FORWARD plan\n");
    } else if (_sign == UP_BACKWARD) {
        INFO("- BACKWARD plan\n");
    }

    INFO("------------------------------------------\n");
}