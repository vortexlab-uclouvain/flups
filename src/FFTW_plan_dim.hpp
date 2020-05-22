/**
 * @file FFTW_plan_dim.hpp
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

#ifndef FFTW_PLAN_DIM_HPP
#define FFTW_PLAN_DIM_HPP

#include "Topology.hpp"
#include "defines.hpp"
#include "fftw3.h"

/**
 * @brief A FFTW plan in one dimension
 * 
 */
class FFTW_plan_dim {
   public:
    /**
     * @brief PlanType is the type of plan considered and is computed as the sum of both BoundaryType variables
     * 
     * The integer value associated gives is the priority of processing.
     * We first have to do the real to real transforms, then the padded real to real (mix direction = unbounded + boundary condition),
     * then the periodic (DFT) directions and finally the padded periodic boundary condition.
     * This order is chosen in order to reduce the computational cost.
     * 
     * If a multi-dimension FFT is asked, one plan is created for each dimension as it may be different.
     * If the plans are the same, we keep the first plan issued
     */
    enum PlanType {
        SYMSYM = 2, /**< type real 2 real (DCT / DST) : EE (0) , EO/OE (1) , OO (2) */
        MIXUNB = 5, /**< type unbounded and a symetry condition: UE/EU (4) , UO/OU (5) */
        PERPER = 6, /**< type periodic - periodic: PERPER (6) */
        UNBUNB = 8, /**< type fully unbounded UU (8) */
        EMPTY  = 18 /**< type empty, i.e. this direction is not used */
    };

    /**
     * @brief Type of real plan, this will drive the correction step in the execute function
     * - CORRECTION_NONE: no correction is needed
     * - CORRECTION_DCT: the correction of a DCT is needed (while going forward, put 0 in the flip-flop mode)
     * - CORRECTION_DST: the correction of a DST is needed (forward: shift the modes FORWARD and put 0, backward: shift the mode backward)
     * 
     */
    enum PlanCorrectionType{
        CORRECTION_NONE = 0,
        CORRECTION_DCT = 1,
        CORRECTION_DST = 2
    };

   protected:
    const int  _lda;     /**<@brief the dimension of the solver */
    const bool _isGreen; /**< @brief boolean is true if this plan is for a Green's function */
    const int  _dimID;   /**< @brief the dimension of the plan in the field reference */
    const int  _sign;    /**< @brief FFT_FORWARD (-1) or FFT_BACKWARD(+1) */

    bool   _isr2c       = false; /**< @brief is this plan the one that changes to complex?*/
    bool   _isSpectral  = false; /**< @brief indicate if the Green's function has to be done spectrally (leading to a helmolz problem) */
    int    _fftw_stride = 0;     /**<@brief the memory space between two ffts */
    int    _howmany     = 0;     /**<@brief the number of FFT's to do */
    int    _fieldstart  = 0;     /**< @brief the starting index for the field copy in the direction of the plan*/
    int    _n_in        = 1;     /**< @brief the number of element in the transform*/
    int    _n_out       = 1;     /**< @brief the number of element coming out of the transform*/
    double _symstart    = 0.0;   /**< @brief the first index to be copied for the symmetry done on the Green's function, set to 0 if no symmetry is needed*/
    double _normfact    = 1.0;   /**< @brief factor you need to multiply to get the transform on the right scaling*/
    double _volfact     = 1.0;   /**< @brief volume factor*/
    double _kfact       = 0.0;   /**< @brief multiplication factor to have the correct k numbers*/
    double _koffset     = 0.0;   /**< @brief additive factor to have the correct k numbers*/

    PlanType            _type;                    /**< @brief type of this plan, see #PlanType*/
    BoundaryType*       _bc[2]    = {NULL, NULL}; /**< @brief boundary condition for the ith component [0][i]=LEFT/MIN - [1][i]=RIGHT/MAX*/
    PlanCorrectionType* _corrtype = NULL;         /**< @brief correction type of this plan, see #PlanCorrectionType*/
    bool*               _imult    = NULL;        /**< @brief boolean indicating that we have to multiply by (-i) in forward and (i) in backward*/
    fftw_r2r_kind*      _kind     = NULL;         /**< @brief kind of transfrom to perform (used by r2r and mix plan only)*/
    fftw_plan*          _plan     = NULL;         /**< @brief the array of FFTW plan*/

   public:
    FFTW_plan_dim(const int lda, const int dimID, const double h[3], const double L[3], BoundaryType* mybc[2], const int sign, const bool isGreen);
    ~FFTW_plan_dim();

    void init(const int size[3], const bool isComplex);

    void allocate_plan(const Topology* topo, double* data);
    void correct_plan(const Topology*, double* data);
    void execute_plan(const Topology* topo, double* data) const;

    /**
     * @name Getters - return the value
     * 
     */
    /**@{ */
    inline bool   isSpectral() const { return _isSpectral; }
    inline bool   isr2c() const { return _isr2c; }
    inline bool   imult(const int lia) const { return _imult[lia]; }
    inline bool   isr2c_doneByFFT() const { return _isr2c && (!_isSpectral); }
    inline int    dimID() const { return _dimID; }
    inline int    type() const { return _type; }
    inline double symstart() const { return _symstart; }
    inline double normfact() const { return _normfact; }
    inline double volfact() const { return _volfact; }
    inline double kfact() const { return _kfact; }
    inline double koffset() const { return _koffset; }
    inline void   get_outsize(int* size) const { size[_dimID] = _n_out; };
    inline void   get_fieldstart(int* start) const { start[_dimID] = _fieldstart; };
    inline void   get_isNowComplex(bool* isComplex) const { (*isComplex) = (*isComplex) || _isr2c; };
    /**@} */

    void disp();

   protected:
    void _check_dataAlign(const Topology* topo, double* data) const;
    /**
     * @name Initialization
     */
    /**@{ */
    void _init_real2real(const int size[3], bool isComplex);
    void _init_mixunbounded(const int size[3], bool isComplex);
    void _init_periodic(const int size[3], bool isComplex);
    void _init_unbounded(const int size[3], bool isComplex);
    void _init_empty(const int size[3], bool isComplex);
    /**@} */

    /**
     * @name Plan allocation
     */
    /**@{ */
    void _allocate_plan_real(const Topology* topo, double* data);
    void _allocate_plan_complex(const Topology* topo, double* data);
    /**@} */
};

#endif