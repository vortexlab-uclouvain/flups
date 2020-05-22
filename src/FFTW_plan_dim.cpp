/**
 * @file FFTW_plan_dim.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
 * @param sign the sign of the plan (FLUPS_FORWARD or FLUPS_BACKWARD)
 * @param isGreen boolean to indicate if the plan is intended for Green's function
 */
FFTW_plan_dim::FFTW_plan_dim(const int lda, const int dimID, const double h[3], const double L[3], BoundaryType* mybc[2], const int sign, const bool isGreen) : _lda(lda),
                                                                                                                                                                    _dimID(dimID),
                                                                                                                                                                    _sign(sign),
                                                                                                                                                                    _isGreen(isGreen) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    FLUPS_CHECK(dimID >= 0 && dimID < 3,"we are only creating plans on dim from 0 to 2", LOCATION);

    //-------------------------------------------------------------------------
    // get the boundary conditions for each dimnension
    //-------------------------------------------------------------------------
    // allocate the bc space
    _bc[0] =(BoundaryType*) flups_malloc(sizeof(int)*_lda);
    _bc[1] =(BoundaryType*) flups_malloc(sizeof(int)*_lda);

    //store the other dimension and check if the type is correct
    for (int lia = 0; lia < _lda; lia++) {
        _bc[0][lia] = mybc[0][lia];
        _bc[1][lia] = mybc[1][lia];
    }

    // setup the type of solver, given by the first dimension
    int mytype = _bc[0][0] + _bc[1][0];

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
        _type       = MIXUNB;
        _normfact   = 1.0;
        _volfact    = h[_dimID];
        _kfact      = c_2pi / (4.0 * L[_dimID]);
        _isSpectral = false;
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
    } else if (mytype == EMPTY) {
        _type = EMPTY;
        // chosen to have no influence
        _normfact   = 1.0;
        _volfact    = 1.0;
        _kfact      = 0.0;
        _isSpectral = false;
    } else {
        FLUPS_ERROR("Invalid combination of BCs", LOCATION);
    }

    //-------------------------------------------------------------------------
    // Get type and mult factors
    //-------------------------------------------------------------------------
    for(int lia= 0 ; lia<_lda; lia++){
        FLUPS_CHECK(_bc[0][lia] + _bc[1][lia] <= _type, "dimension %d's bc = %d %d is not compatible with the plan choosen = %d", lia, _bc[0][lia], _bc[1][lia], _type,LOCATION);
    }
    END_FUNC;
}

/**
 * @brief Destroy the fftw plan
 * 
 */
FFTW_plan_dim::~FFTW_plan_dim() {
    BEGIN_FUNC;

    if (_type == SYMSYM || _type == MIXUNB) {
        // if the solver is SYMSYM or MIXUNB, each dimension has its own plan
        for (int lia = 0; lia < _lda; lia++) {
            if (_plan != NULL) fftw_destroy_plan(_plan[lia]);
        }
    } else {
        // else, the first plan is the same as all the other ones
        if (_plan != NULL) fftw_destroy_plan(_plan[0]);
    }
    // free the allocated arrays
    if (_bc[0] != NULL) flups_free(_bc[0]);
    if (_bc[1] != NULL) flups_free(_bc[1]);
    if (_imult != NULL) flups_free(_imult);
    if (_kind != NULL) flups_free(_kind);
    if (_corrtype != NULL) flups_free(_corrtype);
    if (_plan != NULL) flups_free(_plan);
    END_FUNC;
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
    } else if (_type == EMPTY) {
        FLUPS_INFO_1("No plan required for this direction");
    }
    END_FUNC;
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
    /** - check that the BC given for the different component are compatible,
     *    i.e. the size of the transform is the same for every compoment */
    //-------------------------------------------------------------------------
    // get if the first dimension asks for a type 4 transform?
    bool istype4 = _bc[0][0] != _bc[1][0];
    for(int lia=1; lia<_lda; lia++){
        // a boundary condition of type4 = left != right has to be the case for EVERY component
        bool type4 = _bc[0][lia] != _bc[1][lia];
        if((type4 && !istype4) || (!type4 && istype4) ){
            FLUPS_ERROR("one component has an EVEN-ODD condition, while one of the other uses EVEN-EVEN or ODD-ODD, which is not supported",LOCATION);
        }
    }

    //-------------------------------------------------------------------------
    /** - get #_fieldstart, #_isr2c and #_symstart for Green */
    //-------------------------------------------------------------------------
    _symstart   = 0;  // if no symmetry is needed, set to 0
    _fieldstart = 0;
    _isr2c      = false;

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (2.0 * size[_dimID]);

    //-------------------------------------------------------------------------
    /** - Get the #_kind of Fourier transforms, the #_koffset for each dimension */
    //-------------------------------------------------------------------------
    _imult    = (bool*)flups_malloc(sizeof(bool) * _lda);
    _kind     = (fftw_r2r_kind*)flups_malloc(sizeof(fftw_r2r_kind) * _lda);
    _corrtype = (PlanCorrectionType*)flups_malloc(sizeof(int) * _lda);

    // because of the constrain on the BC, we only the kind argument is linked to the lia
    // while the other values (n_in, n_out and koffset) will remain unchanged accross the lda
    // yet its easier to read if we set them lda times...
    for (int lia = 0; lia < _lda; lia++) {
        if (_isGreen) {
            _corrtype[lia] = CORRECTION_NONE;
            _imult[lia]    = false;
            // if we are doing odd-even we have to use shifted FFTW plans
            if (_bc[0][lia] != _bc[1][lia]) {
                // we would go for a DCT/DST type III
                // -> the size of unknows: DST missing first point, DCT missing last one
                _n_in  = size[_dimID];
                _n_out = size[_dimID];
                // -> the modes are shifted by 1/2
                _koffset = 0.5;
            } else {
                // we go for DST/DCT of type I or III
                // -> we have to add one information because of the vertex-centered
                _n_in  = size[_dimID] + 1;
                _n_out = size[_dimID] + 1;
                // no shift in the mode is required
                _koffset = 0.0;
            }
            return;
        } else if (_bc[0][lia] == EVEN) {  // We have a DCT
            // the information coming in does not change
            _n_in = size[_dimID];
            // we do a DCT, so no imult
            _imult[lia] = false;
            if (_bc[1][lia] == EVEN) {
                // -> we add the flip-flop mode by hand
                _n_out = size[_dimID] + 1;
                // the correction is the one of the DCT = put 0 in the flip-flop mode
                _corrtype[lia] = CORRECTION_DCT;
                _koffset       = 0.0;
                // choose the correct type
                if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_REDFT10;   // DCT type II
                if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_REDFT01;  // DCT type III
            } else if (_bc[1][lia] == ODD) {
                // no additional mode is required
                _n_out = size[_dimID];
                // no correction is needed for the types 4 but an offset of 1/2 in fourier
                _corrtype[lia] = CORRECTION_NONE;
                _koffset       = 0.5;
                // always the samed DCT
                if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_REDFT11;   // DCT type IV
                if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_REDFT11;  // DCT type IV
            }
        } else if (_bc[0][lia] == ODD) {  // We have a DST
                                          // the information coming in does not change
            _n_in = size[_dimID];
            // we do a DST, so no imult
            _imult[lia] = true;
            if (_bc[1][lia] == ODD) {
                // -> we add the 0 mode by hand
                _n_out = size[_dimID] + 1;
                // the correction is the one of the DST = put 0 in the 0 mode
                _corrtype[lia] = CORRECTION_DST;
                _koffset       = 0.0;
                // always the correct DST
                if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_RODFT10;   // DST type II
                if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_RODFT01;  // DST type III
            } else if (_bc[1][lia] == EVEN) {
                // no additional mode is required
                _n_out = size[_dimID];
                // no correction is needed for the types 4 but an offset of 1/2 in fourier
                _corrtype[lia] = CORRECTION_NONE;
                _koffset       = 0.5;
                // always the samed DST
                if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_RODFT11;   // DST type IV
                if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_RODFT11;  // DST type IV
            }
        } else {
            FLUPS_ERROR("unable to init the solver required", LOCATION);
        }
    }
    END_FUNC;
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
    FLUPS_CHECK(isComplex == false, "the data cannot be complex", LOCATION);

    //-------------------------------------------------------------------------
    /** - get the memory details: #_fieldstart and #_isr2c */
    //-------------------------------------------------------------------------
    _isr2c = false;

    if (_isGreen)
        _fieldstart = 0;
    else if (_bc[0][0] == UNB)
        _fieldstart = size[_dimID];  // padding to the left - only the first dim is enough
    else if (_bc[1][0] == UNB)
        _fieldstart = 0;  // padding to the right - only the first dim is enough

    //-------------------------------------------------------------------------
    /** - get the #_symstart if is Green */
    //-------------------------------------------------------------------------
    _symstart = 0;  // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - update #_normfact factor */
    //-------------------------------------------------------------------------
    _normfact *= 1.0 / (4.0 * size[_dimID]);

    //-------------------------------------------------------------------------
    /** - Get the #_kind of Fourier transforms */
    //-------------------------------------------------------------------------
    _imult    = (bool*)flups_malloc(sizeof(bool) * _lda);
    _kind     = (fftw_r2r_kind*)flups_malloc(sizeof(fftw_r2r_kind) * _lda);
    _corrtype = (PlanCorrectionType*)flups_malloc(sizeof(int) * _lda);

    for (int lia = 0; lia < _lda; lia++) {
        if (_isGreen) {
            // the sizes have to be augmented by 1 compared to the cell-centered approach
            _n_in  = 2 * size[_dimID] + 1;
            _n_out = 2 * size[_dimID] + 1;
            // since we do a pure DCT/DST, no offset
            _koffset = 0.0;
            // no correction is needed
            _corrtype[lia] = CORRECTION_NONE;
            // we do a DCT, so no imult
            _imult[lia] = false;
            // The Green function is ALWAYS EVEN - EVEN
            if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_REDFT00;  // DCT type I
            if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_REDFT00;

        } else {
            // we double the size of the data
            _n_in = 2 * size[_dimID];
            // we add a mode for the outgoing dct/dst
            _n_out = 2 * size[_dimID] + 1;
            // no offset after the correction
            _koffset = 0.0;

            if ((_bc[0][lia] == EVEN && _bc[1][lia] == UNB) || (_bc[0][lia] == UNB && _bc[1][lia] == EVEN)) {  // We have a DCT - we are EVEN - EVEN over 2L
                // we need a DCT correction
                _corrtype[lia] = CORRECTION_DCT;
                // we do a DCT, so no imult
                _imult[lia] = false;
                if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_REDFT10;   // DCT type II
                if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_REDFT01;  // DCT type III

            } else if ((_bc[0][lia] == UNB && _bc[1][lia] == ODD) || (_bc[0][lia] == ODD && _bc[1][lia] == UNB)) {  // We have a DST - we are ODD - ODD over 2L
                                                                                                                    // we need a DST correction
                _corrtype[lia] = CORRECTION_DST;
                // we do a DCT, so no imult
                _imult[lia] = true;
                if (_sign == FLUPS_FORWARD) _kind[lia] = FFTW_RODFT10;   // DST type II
                if (_sign == FLUPS_BACKWARD) _kind[lia] = FFTW_RODFT01;  // DST type III
                _koffset = 0.0;
            } else {
                FLUPS_ERROR("unable to init the solver required", LOCATION);
            }
        }
    }
    END_FUNC;
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
    /** - Get the #_koffset factor */
    //-------------------------------------------------------------------------
    _corrtype = (PlanCorrectionType*)flups_malloc(sizeof(int) * _lda);
    _imult    = (bool*)flups_malloc(sizeof(bool) * _lda);
    for (int lia = 0; lia < _lda; lia++) {
        _corrtype[lia] = CORRECTION_NONE;
        // we do a DFT, so no imult
        _imult[lia] = false;
    }
    END_FUNC;
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
    /** - Get the #_koffset factor */
    //-------------------------------------------------------------------------
    _corrtype = (PlanCorrectionType*)flups_malloc(sizeof(int) * _lda);
    _imult    = (bool*)flups_malloc(sizeof(bool) * _lda);
    for (int lia = 0; lia < _lda; lia++) {
        _corrtype[lia] = CORRECTION_NONE;
        // we do a DFT, so no imult
        _imult[lia] = false;
    }
    END_FUNC;
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
    END_FUNC;
}

/**
 * @brief Allocate a plan that only treats real numbers
 * 
 * @note
 * The howmany is computed using the memory information provided
 * 
 * @warning
 * If this plan creation is called with a complex topo, we will do
 * a transform using a stride of 1 on the real part of the array
 * This may happen in the Green's transform if there is a combination of real-real 
 * directions and periodic direction. Hence, only a DCT/DST on the real part is needed
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
    FLUPS_CHECK(data != NULL,"data cannot be null",LOCATION);
    const int memsize[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)}; //the "current" size, corresponding to size_tmp during the dry run, see _init_plansAndTopos

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
    /** - Compute howmany and the stride to run the loop in the execute  */
    //-------------------------------------------------------------------------
    _howmany = 1;
    for (int id = 0; id < _dimID; id++) _howmany *= topo->nloc(id);
    for (int id = _dimID + 1; id < 3; id++) _howmany *= topo->nloc(id);

    //-------------------------------------------------------------------------
    /** - Create the plan  */
    //-------------------------------------------------------------------------
    // we make sure to use only 1 thread, the multi-threading is used in the solver, not inside a plan
    fftw_plan_with_nthreads(1);

    // allocate the plan
    _plan =(fftw_plan*) flups_malloc(sizeof(fftw_plan) * _lda);

    // we initiate the plan with the size #_n_in, because this is the real number of data needed
    for (int lia = 0; lia < _lda; lia++) {
        if (topo->nf() == 1) {
            _fftw_stride = memsize[_dimID];
            _plan[lia]   = fftw_plan_r2r_1d(_n_in, data, data, _kind[lia], FFTW_FLAG);

        } else if (topo->nf() == 2) {
            _fftw_stride = memsize[_dimID] * topo->nf();
            _plan[lia]   = fftw_plan_many_r2r(1, (int*)(&_n_in), 1,
                                            data, NULL, topo->nf(), memsize[_dimID] * topo->nf(),
                                            data, NULL, topo->nf(), memsize[_dimID] * topo->nf(), _kind + lia, FFTW_FLAG);
        }
    }

    FLUPS_INFO("------------------------------------------");
    if (_type == SYMSYM) {
        FLUPS_INFO("## SYMSYM plan created for plan r2r (=%d)", _type);
    } else if (_type == MIXUNB) {
        FLUPS_INFO("## MIXUNB plan created for plan mix (=%d)", _type);
    }
    FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
    FLUPS_INFO("dimID     = %d", _dimID);
    FLUPS_INFO("howmany   = %d", _howmany);
    FLUPS_INFO("fftw stride   = %d", _fftw_stride);
    FLUPS_INFO("size n    = %d", _n_in);
    if (topo->nf() == 1) {
        FLUPS_INFO("plan created with the simple interface");
    } else if (topo->nf() == 2) {
        FLUPS_INFO("plan created with the many interface for non-unit stride)");
    }
    FLUPS_INFO("------------------------------------------");
    END_FUNC;
}

/**
 * @brief allocate a plan that treats complex numbers (r2c or c2c)
 * 
 * @note
 * The howmany is computed using the memory information provided
 * 
 * 
 * @param memsize the size of the data BEFORE THE PLAN is executed
 * @param data memory
 */
void FFTW_plan_dim::_allocate_plan_complex(const Topology *topo, double* data) {
    BEGIN_FUNC;

    assert(data != NULL);

    const int memsize[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)}; //the "current" size, corresponding to size_tmp during the dry run, see _init_plansAndTopos

    if (_isGreen && _type == PERPER) {
        _plan = NULL;

        FLUPS_INFO("------------------------------------------");
        FLUPS_INFO("## no DFT plan created for Green");
        FLUPS_INFO("------------------------------------------");
        return;
    }

    //-------------------------------------------------------------------------
    /** - Compute howmany and the stride to run the loop in the execute  */
    //-------------------------------------------------------------------------
    // Compute howmany
    _howmany = 1;
    for (int id = 0; id < _dimID; id++) _howmany *= topo->nloc(id);
    for (int id = _dimID + 1; id < 3; id++) _howmany *= topo->nloc(id);
    
    // compute the stride
    _fftw_stride = memsize[_dimID];

    // allocate the plan
    _plan =(fftw_plan*) flups_malloc(sizeof(fftw_plan) * _lda);
       
    if (_isr2c) {
        FLUPS_CHECK(topo->nf() == 1, "the nf of the input topology has to be 1 = real topo",LOCATION);

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
        FLUPS_INFO("howmany   = %d", _howmany);
        FLUPS_INFO("fftw stride   = %d", _fftw_stride);
        FLUPS_INFO("size n    = %d", _n_in);
        FLUPS_INFO("------------------------------------------");

        if (_sign == FLUPS_FORWARD) {
            _plan[0] = fftw_plan_dft_r2c_1d(_n_in, data, (fftw_complex*)data, FFTW_FLAG);
        } else {
            _plan[0] = fftw_plan_dft_c2r_1d(_n_in, (fftw_complex*)data, data, FFTW_FLAG);
        }

    } else {
        FLUPS_CHECK(topo->nf() == 2, "the nf of the input topology has to be 1 = real topo",LOCATION);
        FLUPS_INFO("------------------------------------------");
        if (_type == PERPER) {
            FLUPS_INFO("## C2C plan created for plan periodic-periodic (=%d)", _type);
        } else if (_type == UNBUNB) {
            FLUPS_INFO("## C2C plan created for plan unbounded (=%d)", _type);
        }
        if (_sign == FLUPS_FORWARD) {
            FLUPS_INFO("FORWARD transfrom");
        } else if (_sign == FLUPS_BACKWARD) {
            FLUPS_INFO("BACKWARD transfrom");
        }
        FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
        FLUPS_INFO("dimID     = %d", _dimID);
        FLUPS_INFO("howmany   = %d", _howmany);
        FLUPS_INFO("fftw stride   = %d", _fftw_stride);
        FLUPS_INFO("size n    = %d", _n_in);
        FLUPS_INFO("------------------------------------------");

        _plan[0] = fftw_plan_dft_1d(_n_in, (fftw_complex*)data, (fftw_complex*)data, _sign, FFTW_FLAG);
    }

    // the plan is the same in every other direction
    for(int lia=1; lia<_lda;lia++){
        _plan[lia] = _plan[lia-1];
    }

    END_FUNC;
}

/**
 * @brief check that every starting pointer in a direction is well-aligned for the FFTW requirement
 * 
 * @warning to access the memory, we cannot use #_howmany since it is based on the local size of the topo on the input.
 * Then, we have to use the memdim() function of the Topology
 * 
 * @param topo 
 * @param data 
 */
void FFTW_plan_dim::_check_dataAlign(const Topology* topo, double* data) const {
#ifndef NDEBUG
    const size_t howmany = _howmany;
    const size_t onmax   = _howmany * _lda;
    const size_t memdim = topo->memdim();

    for (size_t id = 0; id < onmax; id++) {
        // get the current index
        size_t io = id%_howmany;
        size_t lia = id/_howmany;
        // get the memory
        double* mydata;
        if (_type == SYMSYM || _type == MIXUNB) {
            mydata = data + lia* memdim + io * _fftw_stride;
        } else if (_type == PERPER || _type == UNBUNB) {
            if (_isr2c) {
                mydata = data + lia* memdim + io * _fftw_stride;
            } else {
                mydata = data + lia* memdim + io * _fftw_stride * 2;
            }
        }
        // check the alignment
        FLUPS_CHECK(fftw_alignment_of(mydata) == 0, "data for FFTW have to be aligned on the FFTW alignement! Alignment is %d with id = %d and fftw_stride = %d", fftw_alignment_of(mydata), id, _fftw_stride, LOCATION);
    }
#endif
}

/**
 * @brief corrects the plan executed depending on #_corrtype and #_sign.
 * 
 * This function resets the correct mode at the correct place in the Topology
 * If going forward:
 * - the DST correction sets the 0-mode to 0 and shifts the modes by 1 on the right  (i-> i+1)
 * - the DCT correction sets the flip-flop mode to 0
 * 
 * If going backward:
 * - the DCT correction is not needed
 * - the DST correction shifts the mode to the left (i-> i-1)
 * 
 * @param data 
 */
void FFTW_plan_dim::correct_plan(const Topology* topo, double* data) {
    BEGIN_FUNC;
    // check the data alignment
    _check_dataAlign(topo,data);

    const int    nloc        = topo->nloc(topo->axis());
    const size_t howmany     = _howmany;
    const size_t memdim      = topo->memdim();
    const size_t fftw_stride = (size_t)_fftw_stride;

    for (int lia = 0; lia < _lda; lia++) {
        // get the starting point of the
        opt_double_ptr mydata = data + lia * memdim;
        // if we need a DCT correction and that we are doing forward (backward doesn't matter)
        if (_corrtype[lia] == CORRECTION_DCT && _sign == FLUPS_FORWARD) {
            // we need to enforce the flip-flop mode to be zero and that's it
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(mydata, fftw_stride, howmany, nloc)
            for (size_t io = 0; io < howmany; io++) {
                // get the memory
                opt_double_ptr dataloc = mydata + io * fftw_stride;
                // reset the flip-flop mode
                dataloc[nloc - 1]       = 0.0;
            }
        }
        else if(_corrtype[lia] == CORRECTION_DST && _sign == FLUPS_FORWARD){
            // we need to enforce the the mode 0 + shift everything from i to i+1
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(mydata, fftw_stride, howmany, nloc)
            for (size_t io = 0; io < howmany; io++) {
                // get the memory
                opt_double_ptr dataloc = mydata + io * fftw_stride;
                // shift everything i -> i+1
                for(int ii=nloc-2; ii >= 0; ii--){
                    dataloc[ii+1] = dataloc[ii];
                }
                dataloc[0] = 0.0;
            }

        }
        else if(_corrtype[lia] == CORRECTION_DST && _sign == FLUPS_BACKWARD){
            // we need to enforce the the mode 0 + shift everything from i to i-1
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(mydata, fftw_stride, howmany, nloc)
            for (size_t io = 0; io < howmany; io++) {
                // get the memory
                opt_double_ptr dataloc = mydata + io * fftw_stride;
                // shift everything i -> i+1
                for(int ii=1; ii < nloc; ii++){
                    dataloc[ii-1] = dataloc[ii];
                }
            }

        }
    }
    END_FUNC;
}

/**
 * @brief Executes the plan for a given Topology on a given data
 * 
 * The transform is done in-place on the data array
 * Every transform is done as a 1 thread 1d transform.
 * The multi-threading is used to perfom several FFT's at once
 * 
 * @warning to access the memory, we cannot use #_howmany since it is based on the local size of the topo on the input.
 * Then, we have to use the memdim() function of the Topology
 * 
 */
void FFTW_plan_dim::execute_plan(const Topology* topo, double* data) const {
    BEGIN_FUNC;

    FLUPS_CHECK(!_isSpectral, "Trying to execute a plan for data which is already spectral", LOCATION);
    FLUPS_CHECK(topo->lda() == _lda, "The given topology's lda does not match with the initialisation one", LOCATION);

    if (_type == SYMSYM) {
        FLUPS_INFO(">> Doing plan real2real for dim %d", _dimID);
    } else if (_type == MIXUNB) {
        FLUPS_INFO(">> Doing plan mix for dim %d", _dimID);
    } else if (_type == PERPER) {
        FLUPS_INFO(">> Doing plan periodic-periodic for dim %d", _dimID);
    } else if (_type == UNBUNB) {
        FLUPS_INFO(">> Doing plan unbounded for dim %d", _dimID);
    } else if (_type == EMPTY) {
        FLUPS_INFO(">> Doing no plan for dim %d", _dimID);
        return;
    }

    // copy the variable to avoid issues while compiling using openMP and gcc
    const size_t howmany     = _howmany;
    const size_t onmax       = _howmany * _lda;
    const size_t fftw_stride = (size_t)_fftw_stride;
    const size_t memdim      = topo->memdim();
    // get the plan pointer
    const fftw_plan* plan = _plan;

    //-------------------------------------------------------------------------
    /** - check the alignment if needed. Cannot be done inside the loop when compiling with GCC and default(none) */
    //-------------------------------------------------------------------------
    _check_dataAlign(topo,data);

    //-------------------------------------------------------------------------
    /** - run the plan on each FFT  */
    //-------------------------------------------------------------------------
    // incomming arrays depends if we are a complex switcher or not
    if (_type == SYMSYM || _type == MIXUNB) {  // R2R
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
        for (size_t id = 0; id < onmax; id++) {
            size_t lia = id / howmany;
            size_t io  = id % howmany;
            // get the memory
            double* mydata = (double*)data + lia * memdim + io * fftw_stride;
            // execute the plan on it
            fftw_execute_r2r(plan[lia], (double*)mydata, (double*)mydata);
        }
    } else if (_type == PERPER || _type == UNBUNB) {
        if (_isr2c) {
            if (_sign == FLUPS_FORWARD) {  // DFT - R2C
                FLUPS_CHECK(topo->nf() == 1, "nf should be 1 at this stage", LOCATION);
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
                for (size_t id = 0; id < onmax; id++) {
                    size_t lia = id / howmany;
                    size_t io  = id % howmany;
                    // get the memory
                    double* mydata = (double*)data + lia * memdim + io * fftw_stride;
                    // execute the plan on it
                    fftw_execute_dft_r2c(plan[lia], (double*)mydata, (fftw_complex*)mydata);
                }
            } else {  // DFT - C2R
                FLUPS_CHECK(topo->nf() == 2, "nf should be 2 at this stage", LOCATION);
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
                for (size_t id = 0; id < onmax; id++) {
                    size_t lia = id / howmany;
                    size_t io  = id % howmany;
                    // WARNING the stride is given in the input size =  REAL => id * _fftw_stride/2 * nf = id * _fftw_stride
                    double* mydata = (double*)data + lia * memdim + io * fftw_stride;
                    // execute the plan on it
                    fftw_execute_dft_c2r(plan[lia], (fftw_complex*)mydata, (double*)mydata);
                }
            }

        } else {  // DFT
            FLUPS_CHECK(topo->nf() == 2, "nf should be 2 at this stage", LOCATION);
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
            for (size_t id = 0; id < onmax; id++) {
                size_t lia = id / howmany;
                size_t io  = id % howmany;
                // we access complex info with a fftw_stride real
                double* mydata = (double*)data + lia * memdim + io * fftw_stride * 2;
                // execute the plan on it
                fftw_execute_dft(plan[lia], (fftw_complex*)mydata, (fftw_complex*)mydata);
            }
        }
    }
    END_FUNC;
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
    for (int lia = 0; lia < _lda; lia++) {
        char msg[512];
        sprintf(msg,"- bc = {");
        if (_bc[0][lia] == EVEN) {
            // FLUPS_INFO(EVEN ,");
            sprintf(msg,"%s EVEN",msg);
        } else if (_bc[0][lia] == ODD) {
            // FLUPS_INFO("- bc = { ODD  ,");
            sprintf(msg,"%s ODD",msg);
        } else if (_bc[0][lia] == UNB) {
            // FLUPS_INFO("- bc = { UNB  ,");
            sprintf(msg,"%s UNB",msg);
        } else if (_bc[0][lia] == PER) {
            // FLUPS_INFO("- bc = { PER  ,");
            sprintf(msg,"%s PER",msg);
        }
        if (_bc[1][lia] == EVEN) {
            // FLUPS_INFO(" EVEN}");
            sprintf(msg,"%s , EVEN}",msg);
        } else if (_bc[1][lia] == ODD) {
            // FLUPS_INFO(" ODD}");
            sprintf(msg,"%s , ODD}",msg);
        } else if (_bc[1][lia] == UNB) {
            // FLUPS_INFO(" UNB}");
            sprintf(msg,"%s , UNB}",msg);
        } else if (_bc[1][lia] == PER) {
            sprintf(msg,"%s , PER}",msg);
            // FLUPS_INFO(" PER}");
        }
        FLUPS_INFO(msg);
        if ((_type == SYMSYM && !_isGreen) || _type == MIXUNB) {
            if (_kind[lia] == FFTW_REDFT00) {
                FLUPS_INFO("- kind = REDFT00 = DCT type I");
            }
            if (_kind[lia] == FFTW_REDFT10) {
                FLUPS_INFO("- kind = REDFT10 = DCT type II");
            }
            if (_kind[lia] == FFTW_REDFT01) {
                FLUPS_INFO("- kind = REDFT01 = DCT type III");
            }
            if (_kind[lia] == FFTW_REDFT11) {
                FLUPS_INFO("- kind = REDFT11 = DCT type IV");
            }
            if (_kind[lia] == FFTW_RODFT00) {
                FLUPS_INFO("- kind = RODFT00 = DST type I");
            }
            if (_kind[lia] == FFTW_RODFT10) {
                FLUPS_INFO("- kind = RODFT10 = DST type II");
            }
            if (_kind[lia] == FFTW_RODFT01) {
                FLUPS_INFO("- kind = RODFT01 = DST type III");
            }
            if (_kind[lia] == FFTW_RODFT11) {
                FLUPS_INFO("- kind = RODFT11 = DST type IV");
            }
        }
    }
    FLUPS_INFO("- dimID      = %d", _dimID);
    FLUPS_INFO("- is Green   ? %d", _isGreen);
    FLUPS_INFO("- s2Complex  ? %d", _isr2c);
    FLUPS_INFO("- n_in       = %d", _n_in);
    FLUPS_INFO("- n_out      = %d", _n_out);
    FLUPS_INFO("- fieldstart = %d", _fieldstart);
    FLUPS_INFO("- isSpectral ? %d", _isSpectral);
    if (_sign == FLUPS_FORWARD) {
        FLUPS_INFO("- FORWARD plan");
    } else if (_sign == FLUPS_BACKWARD) {
        FLUPS_INFO("- BACKWARD plan");
    }

    FLUPS_INFO("------------------------------------------");
    END_FUNC;
}