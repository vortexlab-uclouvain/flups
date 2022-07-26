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


#include "FFTW_plan_dim_node.hpp"

/**
 * @brief Initialize for a real to real plan
 * 
 * @param size 
 * @param isComplex 
 * 
 * -------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim_node::init_real2real_(const int size[3], const bool isComplex) {
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
    for (int lia = 1; lia < lda_; lia++) {
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
    /** - Get the #normfact_  The normfactor is independant of the component but depend on the number of point we give to the fft*/
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (2.0 * ( size[dimID_] - 1));

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
            //In this case, the green functions are spectral. There is no need for correction
            corrtype_[lia] = CORRECTION_NONE;
            imult_[lia]    = false;

            if (bc_[0][lia] == EVEN) {
                if (bc_[1][lia] == EVEN) {
                    n_in_[lia]  = size[dimID_];
                    n_out_ = size[dimID_];
                    koffset_   = 0.0;

                } else if (bc_[1][lia] == ODD) {
                    // We need to remove the 0 at the end of the input
                    n_in_[lia]  = size[dimID_] - 1;
                    n_out_ = size[dimID_];
                    koffset_   = 0.5;
                }

            }else if (bc_[0][lia] == ODD) {  // We have a DST
                if (bc_[1][lia] == ODD) {
                    n_in_[lia]  = size[dimID_] - 2;
                    n_out_ = size[dimID_];
                    koffset_   = 0.0;

                } else if (bc_[1][lia] == EVEN) {
                    // no additional mode is required
                    n_in_[lia]  = size[dimID_] - 1;
                    n_out_ = size[dimID_];
                    koffset_       = 0.5;
                    // always the samed DST
                }

            }

            return;
        
        //-------------------------------------------------------------------------
        /** - Take care of the DCTs                                              */
        //-------------------------------------------------------------------------
            } else if (bc_[0][lia] == EVEN) {
                // we do a DCT, so no imult
                imult_[lia] = false;

                if (bc_[1][lia] == EVEN) {
                    // the information coming in does not change
                    n_in_[lia]  = size[dimID_];
                    n_out_ = size[dimID_];

                    // No correction or offset needed
                    corrtype_[lia] = CORRECTION_NONE;
                    koffset_       = 0.0;

                    // Both the forward and the backward tranform uses a REDFT00 transform
                    kind_[lia] = FFTW_REDFT00;  // DCT type I

                } else if (bc_[1][lia] == ODD) {
                    // We need to remove the 0 at the end of the input
                    n_in_[lia] = size[dimID_] - 1;
                    n_out_     = size[dimID_];

                    koffset_ = 0.0;
                    // always the samed DCT
                    if (sign_ == FLUPS_FORWARD) {
                        // DCT type III
                        kind_[lia] = FFTW_REDFT01;
                        // DCT type III: must set the flip-flop to 0
                        corrtype_[lia] = CORRECTION_FLIPFLOP;
                    } else if (sign_ == FLUPS_BACKWARD) {
                        // DCT type II (inverse of the type III)
                        kind_[lia]     = FFTW_REDFT10;
                        // no correction is needed
                        corrtype_[lia] = CORRECTION_NONE;
                    }
                }

                //-------------------------------------------------------------------------
                /** - Take care of the DSTs                                              */
                //-------------------------------------------------------------------------
            } else if (bc_[0][lia] == ODD) {  // We have a DST

                // we do a DST, so imult
                imult_[lia] = true;
                if (bc_[1][lia] == ODD) {
                    // -> we remove the first and the last data, as FFTW don't need them
                    n_in_[lia]  = size[dimID_] - 2;
                    n_out_ = size[dimID_];

                    // The first data of the memory is not given to fftw
                    fftwstart_[lia] = 1;
                    koffset_        = 0.0;
                    fieldstart_     = 0;

                    // always the correct DST
                    if (sign_ == FLUPS_FORWARD) {
                        // DST type I
                        kind_[lia]     = FFTW_RODFT00;
                        // sets the zero mode to 0 + the flipflop mode
                        corrtype_[lia] = CORRECTION_ZEROMODE + CORRECTION_FLIPFLOP;
                    }
                    if (sign_ == FLUPS_BACKWARD) {
                        kind_[lia]     = FFTW_RODFT00;  // DST type I
                        // no correction
                        corrtype_[lia] = CORRECTION_NONE;
                    }

                } else if (bc_[1][lia] == EVEN) {
                    // The first data is not need by FFTW, we have to manually remove it
                    n_in_[lia] = size[dimID_] - 1;
                    n_out_     = size[dimID_];

                    fftwstart_[lia] = 1;
                    koffset_        = 0.0;
                    fieldstart_     = 0;
                    // always the samed DST
                    if (sign_ == FLUPS_FORWARD) {
                        // DST type III
                        kind_[lia]     = FFTW_RODFT01;
                        // shifts all the mode to the left (i -> i-1) + sets flipflop
                        corrtype_[lia] = CORRECTION_FLIPFLOP + CORRECTION_SHIFTLEFT;
                    }
                    if (sign_ == FLUPS_BACKWARD) {
                        // DST type II
                        kind_[lia]     = FFTW_RODFT10;
                        // shifts all the mode to the right (i -> i+1) + sets zeromode
                        corrtype_[lia] = CORRECTION_ZEROMODE + CORRECTION_SHIFTRIGHT;
                    }
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
void FFTW_plan_dim_node::init_mixunbounded_(const int size[3], const bool isComplex) {
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
        fieldstart_ = size[dimID_] - 1;  // padding to the left - only the first dim is enough
    else if (bc_[1][0] == UNB)
        fieldstart_ = 0;  // padding to the right - only the first dim is enough

    //-------------------------------------------------------------------------
    /** - get the #symstart_ if is Green */
    //-------------------------------------------------------------------------
    symstart_ = 0;  // if no symmetry is needed, set to 0

    //-------------------------------------------------------------------------
    /** - Get the #kind_ of Fourier transforms */
    //-------------------------------------------------------------------------
    kind_     = (fftw_r2r_kind*)m_calloc(sizeof(fftw_r2r_kind) * lda_);

    //-------------------------------------------------------------------------
    /** - Get the #normfact_  The normfactor is independant of the component but depend on the number of point we give to the fft*/
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (2.0 *((2 * size[dimID_] - 1) - 1.0));

    for (int lia = 0; lia < lda_; lia++) {
        if (isGreen_) {
            // We have a DCT - we are EVEN - EVEN over 2L
            if ((bc_[0][lia] == EVEN && bc_[1][lia] == UNB) || (bc_[0][lia] == UNB && bc_[1][lia] == EVEN)) {
                // In node centered, we need to remove one point when doubling the domain
                n_in_[lia]  = 2 * size[dimID_] - 1;
                n_out_ = n_in_[lia];
                // since we do a pure DCT/DST, no offset
                koffset_ = 0.0;
            } else if ((bc_[0][lia] == UNB && bc_[1][lia] == ODD) || (bc_[0][lia] == ODD && bc_[1][lia] == UNB)) {
                //In node centered, we need to remove one point when doubling the domain
                n_in_[lia]  = (2 * size[dimID_] - 1);
                n_out_ = n_in_[lia];
                // since we do a pure DCT/DST, no offset
                koffset_ = 0.0;
            }
            // no correction is needed
            corrtype_[lia] = CORRECTION_NONE;
            // we do a DCT, so no imult
            imult_[lia] = false;
            // The Green function is ALWAYS EVEN - EVEN
            if (sign_ == FLUPS_FORWARD) kind_[lia] = FFTW_REDFT00;  // DCT type I
            if (sign_ == FLUPS_BACKWARD) kind_[lia] = FFTW_REDFT00;

        } else {
            if ((bc_[0][lia] == EVEN && bc_[1][lia] == UNB) || (bc_[0][lia] == UNB && bc_[1][lia] == EVEN)) {  // We have a DCT - we are EVEN - EVEN over 2L
                // In node centered, we need to remove one point when doubling the domain
                n_in_[lia] = 2 * size[dimID_] - 1 ;
                // we add a mode for the outgoing dct/dst
                n_out_ = n_in_[lia];
                // no offset after the correction
                koffset_ = 0.0;
                // we need a DCT correction
                corrtype_[lia] = CORRECTION_NONE;
                // we do a DCT, so no imult
                imult_[lia] = false;
                if (sign_ == FLUPS_FORWARD) kind_[lia] = FFTW_REDFT00;   // DCT type I
                if (sign_ == FLUPS_BACKWARD) kind_[lia] = FFTW_REDFT00;  // DCT type I

            } else if ((bc_[0][lia] == UNB && bc_[1][lia] == ODD) || (bc_[0][lia] == ODD && bc_[1][lia] == UNB)) {
                // We have a DST - we are ODD - ODD over 2L
                // we double the size of the data
                n_in_[lia] = (2 * size[dimID_] - 1) - 2;
                // n_out_ must be equal to the n_out_ of the UNB EVEN transform
                n_out_ = n_in_[lia] + 2;
                // no offset after the correction
                koffset_ = 0.0;
                // we do a DCT, so no imult
                imult_[lia]     = true;
                fftwstart_[lia] = 1;
                if (sign_ == FLUPS_FORWARD) {
                    kind_[lia]     = FFTW_RODFT00;  // DST type I
                    corrtype_[lia] = CORRECTION_ZEROMODE + CORRECTION_FLIPFLOP;
                }
                if (sign_ == FLUPS_BACKWARD) {
                    kind_[lia]     = FFTW_RODFT00;  // DST type I
                    corrtype_[lia] = CORRECTION_NONE;
                }
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
void FFTW_plan_dim_node::init_periodic_(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the memory details (#n_in_[lia], #n_out_, #fieldstart_, #shiftgreen_ and #_isr2c_)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        n_in_[0]  = size[dimID_] - 1;   // takes n complex, return n complex
        n_out_ = size[dimID_] - 1;
        isr2c_ = false;

    } else {
        n_in_[0]  = size[dimID_] - 1;   // takes n-1 real
        n_out_ = n_in_[0] / 2 + 1;      // return n_in/2 + 1 complex

        isr2c_ = true;
    }
    
    fieldstart_ = 0;

    //-------------------------------------------------------------------------
    /** - get the #symstart_ if is Green */
    //-------------------------------------------------------------------------
    symstart_ = 0;  // if no symmetry is needed, set to 0
    if(isComplex ){
        symstart_ = (size[dimID_] - 1)/2.0;
    }

    //-------------------------------------------------------------------------
    /** - update #normfact_ factor */
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (size[dimID_] - 1);

    //-------------------------------------------------------------------------
    /** - Get the #koffset_ factor */
    //-------------------------------------------------------------------------
    for (int lia = 0; lia < lda_; lia++) {
        corrtype_[lia] = CORRECTION_NONE;
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
 * 
 *--------------------------------------
 * We do the following operations:
 */
void FFTW_plan_dim_node::init_unbounded_(const int size[3], const bool isComplex) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the memory details (#n_in_[lia], #n_out_, #fieldstart_, #shiftgreen_ and #_isr2c_)  */
    //-------------------------------------------------------------------------
    if (isComplex) {
        n_in_[0]  = 2 * (size[dimID_] - 1);  // takes 2n complex, return 2n complex
        n_out_ = 2 * (size[dimID_] - 1);

        isr2c_ = false;
    } else {
        n_in_[0]  = 2 * (size[dimID_] - 1);  // takes 2n real (because of the padding)
        n_out_ = n_in_[0] / 2 + 1;     // return n_in/2 + 1 complex

        isr2c_ = true;
    }

    fieldstart_ = 0;
    //------------------------------------------------------------------------- .///
    /** - get the #symstart_ if is Green */
    //-------------------------------------------------------------------------
    symstart_ = size[dimID_] - 1;
    //-------------------------------------------------------------------------
    /** - update #normfact_ factor */
    //-------------------------------------------------------------------------
    normfact_ *= 1.0 / (2.0 * (size[dimID_] - 1 ));
    //-------------------------------------------------------------------------
    /** - Get the #koffset_ factor */
    //-------------------------------------------------------------------------
    for (int lia = 0; lia < lda_; lia++) {
        corrtype_[lia] = CORRECTION_NONE;
        // we do a DFT, so no imult
        imult_[lia] = false;
    }
    END_FUNC;
}