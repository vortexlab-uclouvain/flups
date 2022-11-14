/**
 * @file FFTW_plan_dim_cell.cpp
 * @copyright Copyright © UCLouvain 2020
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include "FFTW_plan_dim_cell.hpp"

/**
 * @brief Initialize for a real to real plan
 * 
 * @param size 
 * @param isComplex 
 * 
 * -------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim_cell::init_real2real_(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    FLUPS_CHECK(isComplex == false,"the data cannot be complex");

    //-------------------------------------------------------------------------
    /** - check that the BC given for the different component are compatible,
     *    i.e. the size of the transform is the same for every compoment */
    //-------------------------------------------------------------------------
    // get if the first dimension asks for a type 4 transform?
    bool istype4 = bc_[0][0] != bc_[1][0];
    for(int lia=1; lia<lda_; lia++){
        // a boundary condition of type4 = left != right has to be the case for EVERY component
        bool type4 = bc_[0][lia] != bc_[1][lia];
        bool isOk = !((type4 && !istype4) || (!type4 && istype4));
        FLUPS_CHECK(isOk, "one component has an EVEN-ODD condition, while one of the other uses EVEN-EVEN or ODD-ODD, which is not supported");
    }

    //-------------------------------------------------------------------------
    /** - get #fieldstart_, #isr2c_ and #symstart_ for Green */
    //-------------------------------------------------------------------------
    symstart_   = 0;  // if no symmetry is needed, set to 0
    fieldstart_ = 0;
    isr2c_      = false;

    //-------------------------------------------------------------------------
    /** - update #normfact_ factor */
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (2.0 * (size[dimID_]));

    //-------------------------------------------------------------------------
    /** - Get the #kind_ of Fourier transforms, the #koffset_ for each dimension */
    //-------------------------------------------------------------------------
    kind_     = (fftw_r2r_kind*)m_calloc(sizeof(fftw_r2r_kind) * lda_);

    // because of the constrain on the BC, we only the kind argument is linked to the lia
    // while the other values (n_in, n_out and koffset) will remain unchanged accross the lda
    // yet its easier to read if we set them lda times...
    for (int lia = 0; lia < lda_; lia++) {
        //-------------------------------------------------------------------------
        /** - Take care of the Green function                                    */
        //-------------------------------------------------------------------------
        if (isGreen_) {
            postpro_type_[lia]  = POSTPRO_NONE;
            fftwstart_in_[lia]  = 0;
            fftwstart_out_[lia] = 0;
            imult_[lia]         = false;

            // if we are doing odd-even we have to use shifted FFTW plans
            if (bc_[0][lia] != bc_[1][lia]) {
                // we would go for a DCT/DST type IV
                // -> the size of unknows: DST missing first point, DCT missing last one
                n_in_[lia]  = size[dimID_];
                n_out_ = size[dimID_];
                // -> the modes are shifted by 1/2
                koffset_ = 0.5;
            } else {
                // we go for DST/DCT of type II or III
                // -> we have to add one information because of the vertex-centered
                n_in_[lia]  = size[dimID_];
                n_out_ = size[dimID_] + 1;
                // no shift in the mode is required
                koffset_ = 0.0;
            }
            return;
        //-------------------------------------------------------------------------
        /** - Take care of the DCTs                                              */
        //-------------------------------------------------------------------------
        } else if (bc_[0][lia] == EVEN) {  // We have a DCT
            // the information coming in does not change
            n_in_[lia] = size[dimID_];
            // we do a DCT, so no imult
            imult_[lia]         = false;
            
            // -------- Even - Even --------
            if (bc_[1][lia] == EVEN) {
                // -> we add the flip-flop mode by hand
                n_out_ = size[dimID_] + 1;
                koffset_       = 0.0;
                // choose the correct type
                if (sign_ == FLUPS_FORWARD){
                    kind_[lia] = FFTW_REDFT10;   // DCT type II
                    // the correction is the one of the DCT = put 0 in the flip-flop mode
                    fftwstart_in_[lia]  = 0;
                    fftwstart_out_[lia] = 0;
                    postpro_type_[lia]  = NULL_LAST_POINT;
                }
                if (sign_ == FLUPS_BACKWARD){
                    kind_[lia] = FFTW_REDFT01;  // DCT type III
                    // the correction is the one of the DCT = put 0 in the flip-flop mode
                    fftwstart_in_[lia]  = 0;
                    fftwstart_out_[lia] = 0;
                    postpro_type_[lia] = POSTPRO_NONE ;
                }
            // -------- Even - Odd --------
            } else if (bc_[1][lia] == ODD) {
                // no additional mode is required
                n_out_ = size[dimID_];
                // no correction is needed for the types 4 but an offset of 1/2 in fourier
                koffset_       = 0.5;
                // always the samed DCT
                if (sign_ == FLUPS_FORWARD){
                    // -- DCT type IV --
                    kind_[lia]          = FFTW_REDFT11;
                    fftwstart_in_[lia]  = 0;
                    fftwstart_out_[lia] = 0;
                    postpro_type_[lia]  = POSTPRO_NONE;
                }
                if (sign_ == FLUPS_BACKWARD) {
                    kind_[lia]          = FFTW_REDFT11;  // DCT type IV
                    fftwstart_in_[lia]  = 0;
                    fftwstart_out_[lia] = 0;
                    postpro_type_[lia]  = POSTPRO_NONE;
                }
            }
        //-------------------------------------------------------------------------
        /** - Take care of the DSTs                                              */
        //-------------------------------------------------------------------------
        } else if (bc_[0][lia] == ODD) {  // We have a DST
            // the information coming in does not change
            n_in_[lia] = size[dimID_];
            // we do a DST, so no imult
            imult_[lia] = true;
            // -------- Odd - Odd --------
            if (bc_[1][lia] == ODD) {
                // -> we add the 0 mode by hand
                n_out_ = size[dimID_] + 1;

                koffset_ = 0.0;
                // always the correct DST
                if (sign_ == FLUPS_FORWARD) {
                    //--- DST type-II ---
                    kind_[lia]     = FFTW_RODFT10;  
                    fftwstart_in_[lia] = 0;
                    fftwstart_out_[lia] = 1;
                    postpro_type_[lia] = NULL_FIRST_POINT;
                }
                if (sign_ == FLUPS_BACKWARD) {
                    //--- DST type-III ---
                    kind_[lia]          = FFTW_RODFT01;
                    fftwstart_in_[lia]  = 1;
                    fftwstart_out_[lia] = 0;
                    postpro_type_[lia] = NULL_LAST_POINT;
                }
            // -------- Odd - Even --------
            } else if (bc_[1][lia] == EVEN) {
                // no additional mode is required
                n_out_ = size[dimID_];
                // no correction is needed for the types 4 but an offset of 1/2 in fourier
                fftwstart_in_[lia]  = 0;
                fftwstart_out_[lia] = 0;
                postpro_type_[lia]  = POSTPRO_NONE;
                koffset_       = 0.5;
                // always the samed DST
                if (sign_ == FLUPS_FORWARD) kind_[lia] = FFTW_RODFT11;   // DST type IV
                if (sign_ == FLUPS_BACKWARD) kind_[lia] = FFTW_RODFT11;  // DST type IV
            }
        } else {
            FLUPS_CHECK(false, "unable to init the solver required");
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
void FFTW_plan_dim_cell::init_mixunbounded_(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    FLUPS_CHECK(isComplex == false, "the data cannot be complex");

    //-------------------------------------------------------------------------
    /** - get the memory details: #fieldstart_ and #isr2c_ */
    //-------------------------------------------------------------------------
    isr2c_ = false;

    if (isGreen_)
        fieldstart_ = 0;
    else if (bc_[0][0] == UNB)
        fieldstart_ = size[dimID_];  // padding to the left - only the first dim is enough
    else if (bc_[1][0] == UNB)
        fieldstart_ = 0;  // padding to the right - only the first dim is enough

    //-------------------------------------------------------------------------
    /** - get the #symstart_ if is Green */
    //-------------------------------------------------------------------------
    symstart_ = 0;  // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - update #normfact_ factor */
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (4.0 * size[dimID_]);

    //-------------------------------------------------------------------------
    /** - Get the #kind_ of Fourier transforms */
    //-------------------------------------------------------------------------
    kind_ = (fftw_r2r_kind*)m_calloc(sizeof(fftw_r2r_kind) * lda_);

    for (int lia = 0; lia < lda_; lia++) {
        if (isGreen_) {
            // the sizes have to be augmented by 1 compared to the cell-centered approach
            n_in_[lia] = 2 * size[dimID_] + 1;
            n_out_     = 2 * size[dimID_] + 1;
            // since we do a pure DCT/DST, no offset
            koffset_ = 0.0;
            // no correction is needed
            fftwstart_in_[lia] = 0;
            fftwstart_out_[lia] = 0;
            postpro_type_[lia] = POSTPRO_NONE;
            // we do a DCT, so no imult
            imult_[lia] = false;
            // The Green function is ALWAYS EVEN - EVEN
            if (sign_ == FLUPS_FORWARD) kind_[lia] = FFTW_REDFT00;  // DCT type I
            if (sign_ == FLUPS_BACKWARD) kind_[lia] = FFTW_REDFT00;

        } else {
            // we double the size of the data
            n_in_[lia] = 2 * size[dimID_];
            // we add a mode for the outgoing dct/dst
            n_out_ = 2 * size[dimID_] + 1;
            // no offset after the correction
            koffset_ = 0.0;

            if ((bc_[0][lia] == EVEN && bc_[1][lia] == UNB) || (bc_[0][lia] == UNB && bc_[1][lia] == EVEN)) {  // We have a DCT - we are EVEN - EVEN over 2L
                // we do a DCT, so no imult
                imult_[lia] = false;
                fftwstart_in_[lia] = 0;
                fftwstart_out_[lia] = 0;
                if (sign_ == FLUPS_FORWARD) {
                    kind_[lia]         = FFTW_REDFT10;  // DCT type II
                    postpro_type_[lia] = NULL_LAST_POINT;
                }
                if (sign_ == FLUPS_BACKWARD) {
                    kind_[lia]     = FFTW_REDFT01;  // DCT type III
                    postpro_type_[lia] = POSTPRO_NONE;
                }
            } else if ((bc_[0][lia] == UNB && bc_[1][lia] == ODD) || (bc_[0][lia] == ODD && bc_[1][lia] == UNB)) {  // We have a DST - we are ODD - ODD over 2L
                // we do a DCT, so no imult
                imult_[lia] = true;
                if (sign_ == FLUPS_FORWARD) {
                    kind_[lia]     = FFTW_RODFT10;  // DST type II
                    fftwstart_in_[lia] = 0;
                    fftwstart_out_[lia] = 1;
                    // postpro_type_[lia] = CORRECTION_SHIFTRIGHT + CORRECTION_ZEROMODE;
                    postpro_type_[lia] = NULL_FIRST_POINT;
                }
                if (sign_ == FLUPS_BACKWARD) {
                    kind_[lia]          = FFTW_RODFT01;  // DST type III
                    fftwstart_in_[lia]  = 1;
                    fftwstart_out_[lia] = 0;
                    postpro_type_[lia] = NULL_LAST_POINT;
                }
                koffset_ = 0.0;
            } else {
                FLUPS_CHECK(false, "unable to init the solver required");
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
void FFTW_plan_dim_cell::init_periodic_(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the memory details (#n_in_, #n_out_, #fieldstart_, #shiftgreen_ and #_isr2c_)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        n_in_[0]  = size[dimID_];  // takes n complex, return n complex
        n_out_ = size[dimID_];

        isr2c_ = false;
    } else {
        n_in_[0]  = size[dimID_];   // takes n real
        n_out_ = n_in_[0] / 2 + 1;  // return n_in/2 + 1 complex

        isr2c_ = true;
    }
    
    fieldstart_ = 0;

    //-------------------------------------------------------------------------
    /** - get the #symstart_ if is Green */
    //-------------------------------------------------------------------------
    symstart_ = 0;  // if no symmetry is needed, set to 0
    if(isComplex ){
        symstart_ = size[dimID_]/2.0;
    }

    //-------------------------------------------------------------------------
    /** - update #normfact_ factor */
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (size[dimID_]);

    //-------------------------------------------------------------------------
    /** - Get the #koffset_ factor */
    //-------------------------------------------------------------------------
    for (int lia = 0; lia < lda_; lia++) {
        postpro_type_[lia]  = POSTPRO_NONE;
        fftwstart_in_[lia]  = 0;
        fftwstart_out_[lia] = 0;
        // we do a DFT, so no imult
        imult_[lia] = false;
    }
    END_FUNC;
}

/**
 * @brief initialize the plan for unbounded solvers
 * 
 * @param size the size
 * @param isComplex if the data is already complex
 * 
 *--------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim_cell::init_unbounded_(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the memory details (#n_in_, #n_out_, #fieldstart_, #shiftgreen_ and #_isr2c_)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        n_in_[0]  = 2 * size[dimID_];  // takes 2n complex, return 2n complex
        n_out_ = 2 * size[dimID_];

        isr2c_ = false;
    } else {
        n_in_[0]  = 2 * size[dimID_];  // takes 2n real (because of the padding)
        n_out_ = n_in_[0] / 2 + 1;     // return n_in/2 + 1 complex

        isr2c_ = true;
    }

    fieldstart_ = 0;
    //-------------------------------------------------------------------------
    /** - get the #symstart_ if is Green */
    //-------------------------------------------------------------------------
    symstart_ = size[dimID_];
    //-------------------------------------------------------------------------
    /** - update #normfact_ factor */
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (2.0 * size[dimID_]);
    //-------------------------------------------------------------------------
    /** - Get the #koffset_ factor */
    //-------------------------------------------------------------------------
    for (int lia = 0; lia < lda_; lia++) {
        fftwstart_in_[lia] = 0;
        fftwstart_out_[lia] = 0;
        postpro_type_[lia] = POSTPRO_NONE;
        // we do a DFT, so no imult
        imult_[lia] = false;
    }
    END_FUNC;
}