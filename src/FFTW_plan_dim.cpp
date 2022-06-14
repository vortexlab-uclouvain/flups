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
FFTW_plan_dim::FFTW_plan_dim(const int lda, const int dimID, const double h[3], const double L[3], BoundaryType* mybc[2], const int sign, const bool isGreen) : lda_(lda),
                                                                                                                                                                isGreen_(isGreen),
                                                                                                                                                                dimID_(dimID),
                                                                                                                                                                sign_(sign) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    FLUPS_CHECK(dimID >= 0 && dimID < 3,"we are only creating plans on dim from 0 to 2");

    //-------------------------------------------------------------------------
    // get the boundary conditions for each dimnension
    //-------------------------------------------------------------------------
    // allocate the bc space
    bc_[0] =(BoundaryType*) m_calloc(sizeof(int)*lda_);
    bc_[1] =(BoundaryType*) m_calloc(sizeof(int)*lda_);

    //store the other dimension and check if the type is correct
    for (int lia = 0; lia < lda_; lia++) {
        bc_[0][lia] = mybc[0][lia];
        bc_[1][lia] = mybc[1][lia];
    }

    // setup the type of solver, given by the first dimension
    int mytype = bc_[0][0] + bc_[1][0];

    //-------------------------------------------------------------------------
    // Get type and mult factors
    //-------------------------------------------------------------------------
    if (mytype <= SYMSYM) {
        type_     = SYMSYM;
        volfact_  = 1.0;  // no convolution so no multiplication by h
        kfact_    = c_2pi / (2.0 * L[dimID_]);
        if (isGreen_) isSpectral_ = true;
    } else if (mytype <= MIXUNB) {
        type_       = MIXUNB;
        volfact_    = h[dimID_];
        kfact_      = c_2pi / (4.0 * L[dimID_]);
        isSpectral_ = false;
    } else if (mytype == PERPER) {
        type_     = PERPER;
        volfact_  = 1.0;  // no convolution so no multiplication by h
        kfact_    = c_2pi / (L[dimID_]);
        if (isGreen_) isSpectral_ = true;
    } else if (mytype == UNBUNB) {
        type_     = UNBUNB;
        volfact_  = h[dimID_];
        kfact_    = c_2pi / (2.0 * L[dimID_]);
        isSpectral_ = false;
    } else if (mytype == EMPTY) {
        type_ = EMPTY;
        // chosen to have no influence
        volfact_    = 1.0;
        kfact_      = 0.0;
        isSpectral_ = false;
    } else {
        FLUPS_CHECK(false, "Invalid combination of BCs");
    }

    //-------------------------------------------------------------------------
    // Allocate the component dependent stuffs
    //-------------------------------------------------------------------------
    n_in_     = (int*)m_calloc(sizeof(int) * lda_);
    fftwstart_= (int*)m_calloc(sizeof(int) * lda_);
    // normfact_ = (double*)m_calloc(sizeof(double) * lda_);
    for(int lia = 0; lia < lda_; lia ++){
        n_in_[lia]      = 1;
        fftwstart_[lia] = 0;
    }
    normfact_ = 1.0;
    corrtype_ = (PlanCorrectionType*)m_calloc(sizeof(int) * lda_);
    imult_    = (bool*)m_calloc(sizeof(bool) * lda_);

    for(int lia= 0 ; lia<lda_; lia++){
        FLUPS_CHECK(bc_[0][lia] + bc_[1][lia] <= type_, "dimension %d's bc = %d %d is not compatible with the plan choosen = %d", lia, bc_[0][lia], bc_[1][lia], type_);
    }
    END_FUNC;
}

/**
 * @brief Destroy the fftw plan
 * 
 */
FFTW_plan_dim::~FFTW_plan_dim() {
    BEGIN_FUNC;
    //-------------------------------------------------------------------
    
    if (type_ == SYMSYM || type_ == MIXUNB) {
        // if the solver is SYMSYM or MIXUNB, each dimension has its own plan
        for (int lia = 0; lia < lda_; lia++) {
            if (plan_ != NULL) fftw_destroy_plan(plan_[lia]);
        }
    } else {
        // else, the first plan is the same as all the other ones
        if (plan_ != NULL) fftw_destroy_plan(plan_[0]);
    }
    
    // free the allocated arrays
    if (bc_[0] != NULL) m_free(bc_[0]);
    if (bc_[1] != NULL) m_free(bc_[1]);

    if (n_in_ != NULL) m_free(n_in_);
    if (fftwstart_ != NULL) m_free(fftwstart_);
    if (imult_ != NULL) m_free(imult_);
    if (kind_ != NULL) m_free(kind_);
    if (corrtype_ != NULL) m_free(corrtype_);
    if (plan_ != NULL) m_free(plan_);
    //-------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Initialize the FFTW_plan_dim by performing a 'dry run'
 * 
 * The function redirects to one of the init functions depending on the type:
 * - init_real2real_()
 * - init_mixunbounded_()
 * - init_periodic_()
 * - init_unbounded_()
 * 
 * Each of the sub-function initializes the following variables
 * - #n_in_ the size of data provided as input to the FFTW (i.e. the number of real or complex numbers)
 * - #n_out_ the size of data that comes out of the FFTW
 * - #fieldstart_ the index to start the FFTW (non zero for mixunbounded solvers)
 * - #isr2c_ is true if this plan switches to the complex numbers
 * - #kind_ the kind of FFTW plan to execute (for SYMSYM and MIXUNB plans only)
 * - #symstart_ the symmetry start = id of symmetry, for the Green's function only
 * 
 * @param size the current size of data in during dry run (hence already partially transformed)
 * @param isComplex the current complex state of the data
 */
void FFTW_plan_dim::init(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // sanity checks
    //-------------------------------------------------------------------------
    assert(size[dimID_] >= 0);

    //-------------------------------------------------------------------------
    // redirect to the corresponding subfunction
    //-------------------------------------------------------------------------
    if (type_ == SYMSYM) {
        init_real2real_(size, isComplex);
    } else if (type_ == MIXUNB) {
        init_mixunbounded_(size, isComplex);
    } else if (type_ == PERPER) {
        //this is the only transform that could give a R2C on data and being spectral for green
        init_periodic_(size, isComplex);
    } else if (type_ == UNBUNB) {
        init_unbounded_(size, isComplex);
    } else if (type_ == EMPTY) {
        FLUPS_INFO_1("No plan required for this direction");
    }
    END_FUNC;
}



/**
 * @brief allocate the plan based on the information computed by init_()
 * 
 * The function redirects to one of the init functions depending on the type:
 * - allocate_plan_real_()
 * - allocate_plan_complex_()
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
    if (type_ == SYMSYM || type_ == MIXUNB) {
        allocate_plan_real_(topo, data);
    } else if (type_ == PERPER || type_ == UNBUNB) {
        allocate_plan_complex_(topo, data);
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
void FFTW_plan_dim::allocate_plan_real_(const Topology *topo, double* data) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    FLUPS_CHECK(data != NULL,"data cannot be null");
    const int memsize[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)}; //the "current" size, corresponding to size_tmp during the dry run, see init_plansAndTopos_

    //-------------------------------------------------------------------------
    /** - If is Green and #type_ is SYMSYM, exit */
    //-------------------------------------------------------------------------
    if (isGreen_ && type_ == SYMSYM) {
        plan_ = NULL;

        FLUPS_INFO("------------------------------------------");
        FLUPS_INFO("## no real to real plan created for Green");
        FLUPS_INFO("------------------------------------------");
        return;
    }

    //-------------------------------------------------------------------------
    /** - Compute howmany and the stride to run the loop in the execute  */
    //-------------------------------------------------------------------------
    howmany_ = 1;
    for (int id = 0; id < dimID_; id++) howmany_ *= topo->nloc(id);
    for (int id = dimID_ + 1; id < 3; id++) howmany_ *= topo->nloc(id);

    //-------------------------------------------------------------------------
    /** - Create the plan  */
    //-------------------------------------------------------------------------
    // we make sure to use only 1 thread, the multi-threading is used in the solver, not inside a plan
    fftw_plan_with_nthreads(1);

    // allocate the plan
    plan_ =(fftw_plan*) m_calloc(sizeof(fftw_plan) * lda_);

    // we initiate the plan with the size #n_in_, because this is the real number of data needed
    for (int lia = 0; lia < lda_; lia++) {
        if (topo->nf() == 1) {
            fftw_stride_ = memsize[dimID_];
            plan_[lia]   = fftw_plan_r2r_1d(n_in_[lia], data, data, kind_[lia], FLUPS_FFTW_FLAG);

        } else if (topo->nf() == 2) {
            fftw_stride_ = memsize[dimID_] * topo->nf();
            plan_[lia]   = fftw_plan_many_r2r(1, (int*)(&n_in_[lia]), 1,
                                            data, NULL, topo->nf(), memsize[dimID_] * topo->nf(),
                                            data, NULL, topo->nf(), memsize[dimID_] * topo->nf(), kind_ + lia, FLUPS_FFTW_FLAG);
        }
    }

    FLUPS_INFO("------------------------------------------");
    if (type_ == SYMSYM) {
        FLUPS_INFO("## SYMSYM plan created for plan r2r (=%d)", type_);
    } else if (type_ == MIXUNB) {
        FLUPS_INFO("## MIXUNB plan created for plan mix (=%d)", type_);
    }
    FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
    FLUPS_INFO("dimID     = %d", dimID_);
    FLUPS_INFO("howmany   = %d", howmany_);
    FLUPS_INFO("fftw stride   = %d", fftw_stride_);
    FLUPS_INFO("lda       = %d", lda_);
    FLUPS_INFO(" size n   = %d", n_in_[0]);
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
void FFTW_plan_dim::allocate_plan_complex_(const Topology *topo, double* data) {
    BEGIN_FUNC;

    assert(data != NULL);

    const int memsize[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)}; //the "current" size, corresponding to size_tmp during the dry run, see init_plansAndTopos_

    if (isGreen_ && type_ == PERPER) {
        plan_ = NULL;

        FLUPS_INFO("------------------------------------------");
        FLUPS_INFO("## no DFT plan created for Green");
        FLUPS_INFO("------------------------------------------");
        return;
    }

    //-------------------------------------------------------------------------
    /** - Compute howmany and the stride to run the loop in the execute  */
    //-------------------------------------------------------------------------
    // Compute howmany
    howmany_ = 1;
    for (int id = 0; id < dimID_; id++) howmany_ *= topo->nloc(id);
    for (int id = dimID_ + 1; id < 3; id++) howmany_ *= topo->nloc(id);
    
    // compute the stride
    fftw_stride_ = memsize[dimID_];

    // allocate the plan
    plan_ =(fftw_plan*) m_calloc(sizeof(fftw_plan) * lda_);
       
    if (isr2c_) {
        FLUPS_CHECK(topo->nf() == 1, "the nf of the input topology has to be 1 = real topo");

        FLUPS_INFO("------------------------------------------");
        if (type_ == PERPER) {
            FLUPS_INFO("## R2C plan created for plan periodic-periodic (=%d)", type_);
        } else if (type_ == UNBUNB) {
            FLUPS_INFO("## R2C plan created for plan unbounded (=%d)", type_);
        }
        // FLUPS_INFO("orderedID = %d",orderID_);
        if (sign_ == FLUPS_FORWARD) {
            FLUPS_INFO("FORWARD transfrom");
        } else if (sign_ == FLUPS_BACKWARD) {
            FLUPS_INFO("BACKWARD transfrom");
        }
        FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
        FLUPS_INFO("dimID     = %d", dimID_);
        FLUPS_INFO("howmany   = %d", howmany_);
        FLUPS_INFO("fftw stride   = %d", fftw_stride_);
        FLUPS_INFO("size n    = %d", n_in_[0]);
        FLUPS_INFO("------------------------------------------");

        if (sign_ == FLUPS_FORWARD) {
            plan_[0] = fftw_plan_dft_r2c_1d(n_in_[0], data, (fftw_complex*)data, FLUPS_FFTW_FLAG);
        } else {
            plan_[0] = fftw_plan_dft_c2r_1d(n_in_[0], (fftw_complex*)data, data, FLUPS_FFTW_FLAG);
        }

    } else {
        FLUPS_CHECK(topo->nf() == 2, "the nf of the input topology has to be 1 = real topo");
        FLUPS_INFO("------------------------------------------");
        if (type_ == PERPER) {
            FLUPS_INFO("## C2C plan created for plan periodic-periodic (=%d)", type_);
        } else if (type_ == UNBUNB) {
            FLUPS_INFO("## C2C plan created for plan unbounded (=%d)", type_);
        }
        if (sign_ == FLUPS_FORWARD) {
            FLUPS_INFO("FORWARD transfrom");
        } else if (sign_ == FLUPS_BACKWARD) {
            FLUPS_INFO("BACKWARD transfrom");
        }
        FLUPS_INFO("memsize = %d x %d x %d", memsize[0], memsize[1], memsize[2]);
        FLUPS_INFO("dimID     = %d", dimID_);
        FLUPS_INFO("howmany   = %d", howmany_);
        FLUPS_INFO("fftw stride   = %d", fftw_stride_);
        FLUPS_INFO("size n    = %d", n_in_[0]);
        FLUPS_INFO("------------------------------------------");

        plan_[0] = fftw_plan_dft_1d(n_in_[0], (fftw_complex*)data, (fftw_complex*)data, sign_, FLUPS_FFTW_FLAG);
    }

    // the plan is the same in every other direction
    for(int lia=1; lia<lda_;lia++){
        plan_[lia] = plan_[lia-1];
    }

    END_FUNC;
}

/**
 * @brief check that every starting pointer in a direction is well-aligned for the FFTW requirement
 *
 * @warning to access the memory, we cannot use #howmany_ since it is based on the local size of the topo on the input.
 * Then, we have to use the memdim() function of the Topology
 *
 * @param topo
 * @param data
 */
void FFTW_plan_dim::check_dataAlign_(const Topology* topo, double* data) const {
#ifndef NDEBUG
    const size_t howmany = howmany_;
    const size_t onmax   = howmany_ * lda_;
    const size_t memdim  = topo->memdim();

    for (size_t id = 0; id < onmax; id++) {
        // get the current index
        size_t io  = id % howmany;
        size_t lia = id / howmany;
        // get the memory
        double* mydata = nullptr;
        if (type_ == SYMSYM || type_ == MIXUNB) {
            mydata = data + lia * memdim + io * fftw_stride_;
        } else if (type_ == PERPER || type_ == UNBUNB) {
            if (isr2c_) {
                mydata = data + lia * memdim + io * fftw_stride_;
            } else {
                mydata = data + lia * memdim + io * fftw_stride_ * 2;
            }
        }
        // check the alignment
        FLUPS_CHECK(fftw_alignment_of(mydata) == 0, "data for FFTW have to be aligned on the FFTW alignement! Alignment is %d with id = %zu and fftw_stride = %d", fftw_alignment_of(mydata), id, fftw_stride_);
    }
#endif
}

/**
 * @brief corrects the plan executed depending on #corrtype_ and #sign_.
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
    check_dataAlign_(topo,data);

    const int    nloc        = topo->nloc(topo->axis());
    const size_t howmany     = howmany_;
    const size_t memdim      = topo->memdim();
    const size_t fftw_stride = (size_t)fftw_stride_;

    for (int lia = 0; lia < lda_; lia++) {
        // get the starting point of the
        opt_double_ptr mydata = data + lia * memdim;
        // if we need a DCT correction and that we are doing forward (backward doesn't matter)
        if (corrtype_[lia] == CORRECTION_DCT && sign_ == FLUPS_FORWARD) {
            // we need to enforce the flip-flop mode to be zero and that's it
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(mydata, fftw_stride, howmany, nloc)
            for (size_t io = 0; io < howmany; io++) {
                // get the memory
                opt_double_ptr dataloc = mydata + io * fftw_stride;
                // reset the flip-flop mode
                dataloc[nloc - 1]       = 0.0;
            }
        }
        else if(corrtype_[lia] == CORRECTION_DST && sign_ == FLUPS_FORWARD){
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
        else if (corrtype_[lia] == CORRECTION_NDST && sign_ == FLUPS_FORWARD) {
            // we need to shift everything from i+1 to i
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(mydata, fftw_stride, howmany, nloc)
            for (size_t io = 0; io < howmany; io++) {
                // get the memory
                opt_double_ptr dataloc = mydata + io * fftw_stride;
                // shift everything i+1 -> i
                for (int ii = 1; ii < nloc; ii++) {
                    dataloc[ii - 1] = dataloc[ii];
                }
            }
        } 

        else if (corrtype_[lia] == CORRECTION_DST && sign_ == FLUPS_BACKWARD) {
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
        else if (corrtype_[lia] == CORRECTION_NDST && sign_ == FLUPS_BACKWARD) {
            // we need to shift everything from i to i+1
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(mydata, fftw_stride, howmany, nloc)
            for (size_t io = 0; io < howmany; io++) {
                // get the memory
                opt_double_ptr dataloc = mydata + io * fftw_stride;
                // shift everything i+1 -> i
                for (int ii = nloc-2; ii >= 0; ii--) {
                    dataloc[ii + 1] = dataloc[ii];
                }
                // enforce the first point in the physical domain to be 0
                dataloc[0] = 0.0;
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
 * @warning to access the memory, we cannot use #howmany_ since it is based on the local size of the topo on the input.
 * Then, we have to use the memdim() function of the Topology
 * 
 */
void FFTW_plan_dim::execute_plan(const Topology* topo, double* data) const {
    BEGIN_FUNC;
    FLUPS_CHECK(!isSpectral_, "Trying to execute a plan for data which has already been setup spectraly");
    FLUPS_CHECK(topo->lda() == lda_, "The given topology's lda does not match with the initialisation one");

    if (type_ == SYMSYM) {
        FLUPS_INFO(">> Doing plan real2real for dim %d", dimID_);
    } else if (type_ == MIXUNB) {
        FLUPS_INFO(">> Doing plan mix for dim %d", dimID_);
    } else if (type_ == PERPER) {
        FLUPS_INFO(">> Doing plan periodic-periodic for dim %d", dimID_);
    } else if (type_ == UNBUNB) {
        FLUPS_INFO(">> Doing plan unbounded for dim %d", dimID_);
    } else if (type_ == EMPTY) {
        FLUPS_INFO(">> Doing no plan for dim %d", dimID_);
        return;
    }

    // copy the variable to avoid issues while compiling using openMP and gcc
    const size_t howmany     = howmany_;
    const size_t onmax       = howmany_ * lda_;
    const size_t fftw_stride = (size_t)fftw_stride_;
    const size_t memdim      = topo->memdim();
    // get the plan pointer
    const fftw_plan* plan = plan_;

    //-------------------------------------------------------------------------
    /** - check the alignment if needed. Cannot be done inside the loop when compiling with GCC and default(none) */
    //-------------------------------------------------------------------------
    check_dataAlign_(topo,data);

    //-------------------------------------------------------------------------
    /** - run the plan on each FFT  */
    //-------------------------------------------------------------------------
    // incomming arrays depends if we are a complex switcher or not
    if (type_ == SYMSYM || type_ == MIXUNB) {  // R2R
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
        for (size_t id = 0; id < onmax; id++) {
            size_t lia = id / howmany;
            size_t io  = id % howmany;
            // get the memory
            double* mydata = (double*)data + lia * memdim + io * fftw_stride;
            // execute the plan on it
            fftw_execute_r2r(plan[lia], (double*)mydata + fftwstart_[lia], (double*)mydata + fftwstart_[lia]);
        }
    } else if (type_ == PERPER || type_ == UNBUNB) {
        if (isr2c_) {
            if (sign_ == FLUPS_FORWARD) {  // DFT - R2C
                FLUPS_CHECK(topo->nf() == 1, "nf should be 1 at this stage");
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
                for (size_t id = 0; id < onmax; id++) {
                    size_t lia = id / howmany;
                    size_t io  = id % howmany;
                    // get the memory
                    double* mydata = (double*)data + lia * memdim + io * fftw_stride;
                    // execute the plan on it
                    fftw_execute_dft_r2c(plan[lia], (double*)mydata + fftwstart_[lia], (fftw_complex*)mydata  + fftwstart_[lia]);
                }
            } else {  // DFT - C2R
                FLUPS_CHECK(topo->nf() == 2, "nf should be 2 at this stage");
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
                for (size_t id = 0; id < onmax; id++) {
                    size_t lia = id / howmany;
                    size_t io  = id % howmany;
                    // WARNING the stride is given in the input size =  REAL => id * fftw_stride_/2 * nf = id * fftw_stride_
                    double* mydata = (double*)data + lia * memdim + io * fftw_stride;
                    // execute the plan on it
                    fftw_execute_dft_c2r(plan[lia], (fftw_complex*)mydata + fftwstart_[lia], (double*)mydata + fftwstart_[lia]);
                }
            }

        } else {  // DFT
            FLUPS_CHECK(topo->nf() == 2, "nf should be 2 at this stage");
#pragma omp parallel for proc_bind(close) schedule(static) default(none) firstprivate(plan, data, fftw_stride, onmax, howmany, memdim)
            for (size_t id = 0; id < onmax; id++) {
                size_t lia = id / howmany;
                size_t io  = id % howmany;
                // we access complex info with a fftw_stride real
                double* mydata = (double*)data + lia * memdim + io * fftw_stride * 2;
                // execute the plan on it
                fftw_execute_dft(plan[lia], (fftw_complex*)mydata + fftwstart_[lia], (fftw_complex*)mydata  + fftwstart_[lia]);
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
    FLUPS_INFO("## Plan num created for dimension %d", dimID_);
    if (type_ == SYMSYM) {
        FLUPS_INFO("- type = real2real (=%d)", type_);
    } else if (type_ == MIXUNB) {
        FLUPS_INFO("- type = mix (=%d)", type_);
    } else if (type_ == PERPER) {
        FLUPS_INFO("- type = periodic-periodic (=%d)", type_);
    } else if (type_ == UNBUNB) {
        FLUPS_INFO("- type = unbounded (=%d)", type_);
    }
    for (int lia = 0; lia < lda_; lia++) {
        char msg[512];
        sprintf(msg,"- bc = {");
        if (bc_[0][lia] == EVEN) {
            // FLUPS_INFO(EVEN ,");
            sprintf(msg,"%s EVEN",msg);
        } else if (bc_[0][lia] == ODD) {
            // FLUPS_INFO("- bc = { ODD  ,");
            sprintf(msg,"%s ODD",msg);
        } else if (bc_[0][lia] == UNB) {
            // FLUPS_INFO("- bc = { UNB  ,");
            sprintf(msg,"%s UNB",msg);
        } else if (bc_[0][lia] == PER) {
            // FLUPS_INFO("- bc = { PER  ,");
            sprintf(msg,"%s PER",msg);
        }
        if (bc_[1][lia] == EVEN) {
            // FLUPS_INFO(" EVEN}");
            sprintf(msg,"%s , EVEN}",msg);
        } else if (bc_[1][lia] == ODD) {
            // FLUPS_INFO(" ODD}");
            sprintf(msg,"%s , ODD}",msg);
        } else if (bc_[1][lia] == UNB) {
            // FLUPS_INFO(" UNB}");
            sprintf(msg,"%s , UNB}",msg);
        } else if (bc_[1][lia] == PER) {
            sprintf(msg,"%s , PER}",msg);
            // FLUPS_INFO(" PER}");
        }
        FLUPS_INFO("%s", msg);
        if ((type_ == SYMSYM && !isGreen_) || type_ == MIXUNB) {
            if (kind_[lia] == FFTW_REDFT00) {
                FLUPS_INFO("- kind = REDFT00 = DCT type I");
            }
            if (kind_[lia] == FFTW_REDFT10) {
                FLUPS_INFO("- kind = REDFT10 = DCT type II");
            }
            if (kind_[lia] == FFTW_REDFT01) {
                FLUPS_INFO("- kind = REDFT01 = DCT type III");
            }
            if (kind_[lia] == FFTW_REDFT11) {
                FLUPS_INFO("- kind = REDFT11 = DCT type IV");
            }
            if (kind_[lia] == FFTW_RODFT00) {
                FLUPS_INFO("- kind = RODFT00 = DST type I");
            }
            if (kind_[lia] == FFTW_RODFT10) {
                FLUPS_INFO("- kind = RODFT10 = DST type II");
            }
            if (kind_[lia] == FFTW_RODFT01) {
                FLUPS_INFO("- kind = RODFT01 = DST type III");
            }
            if (kind_[lia] == FFTW_RODFT11) {
                FLUPS_INFO("- kind = RODFT11 = DST type IV");
            }
        }
    }
    FLUPS_INFO("- dimID      = %d", dimID_);
    FLUPS_INFO("- is Green   ? %d", isGreen_);
    FLUPS_INFO("- s2Complex  ? %d", isr2c_);
    FLUPS_INFO("- n_in       = %d", n_in_[0]);
    FLUPS_INFO("- n_out      = %d", n_out_);
    FLUPS_INFO("- fieldstart = %d", fieldstart_);
    FLUPS_INFO("- isSpectral ? %d", isSpectral_);
    if (sign_ == FLUPS_FORWARD) {
        FLUPS_INFO("- FORWARD plan");
    } else if (sign_ == FLUPS_BACKWARD) {
        FLUPS_INFO("- BACKWARD plan");
    }

    FLUPS_INFO("------------------------------------------");
    END_FUNC;
}