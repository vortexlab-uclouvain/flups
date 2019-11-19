/**
 * @file FFTW_plan_dim.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
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
    /**
     * @brief PlanType is the type of plan considered and is computed as the sum of both BoundaryType variables
     * 
     * The integer value associated gives is the priority of processing.
     * We first have to do the real to real transforms, then the padded real to real (mix direction = unbounded + boundary condition),
     * then the periodic (DFT) directions and finally the padded periodic boundary condition.
     * This order is chosen in order to reduce the computational cost.
     */
    enum PlanType {
        SYMSYM = 2, /**< type real 2 real (DCT / DST) : EE (0) , EO/OE (1) , OO (2) */
        MIXUNB = 5, /**< type unbounded and a symetry condition: UE/EU (4) , UO/OU (5) */
        PERPER = 6, /**< type periodic - periodic: PERPER (6) */
        UNBUNB = 8, /**< type fully unbounded UU (8) */
        EMPTY  = 18 /**< type empty, i.e. this direction is not used */
    };

   protected:
    const bool _isGreen; /**< @brief boolean is true if this plan is for a Green's function */
    const int  _dimID;   /**< @brief the dimension of the plan in the field reference */
    const int  _sign;    /**< @brief FFT_FORWARD (-1) or FFT_BACKWARD(+1) */

    bool   _ignoreMode = false; /**< @brief do we have to ignore a mode in the output? k=0 if _shiftgreen=1 or k=end if _shiftgreen = 0*/
    bool   _isr2c      = false; /**< @brief is this plan the one that changes to complex?*/
    bool   _imult      = false; /**< @brief boolean to determine if we have to multiply by (i=sqrt(-1)) or not*/
    bool   _isSpectral = false; /**< @brief indicate if the Green's function has to be done spectrally (leading to a helmolz problem) */
    int    _fftw_stride = 0;
    int    _howmany         = 0;
    int    _fieldstart = 0;     /**< @brief the starting index for the field copy in the direction of the plan*/
    int    _n_in       = 0;     /**< @brief the number of element in the transform*/
    int    _n_out      = 0;     /**< @brief the number of element coming out of the transform*/
    int    _shiftgreen = 0;     /**< @brief the shift to set in the Green's function when doing the convolution*/
    double _symstart   = 0.0;   /**< @brief the first index to be copied for the symmetry done on the Green's function, set to 0 if no symmetry is needed*/
    double _normfact   = 0.0;   /**< @brief factor you need to multiply to get the transform on the right scaling*/
    double _volfact    = 0.0;   /**< @brief volume factor*/
    double _kfact      = 0.0;   /**< @brief multiplication factor to have the correct k numbers*/
    double _koffset    = 0.0;   /**< @brief additive factor to have the correct k numbers*/

    PlanType     _type;  /**< @brief type of this plan, see #PlanType*/
    BoundaryType _bc[2]; /**< @brief boundary condition [0]=LEFT/MIN - [1]=RIGHT/MAX*/

    fftw_r2r_kind _kind;        /**< @brief kind of transfrom to perform (used by r2r and mix plan only)*/
    fftw_plan     _plan = NULL; /**< @brief the actual FFTW plan*/

   public:
    FFTW_plan_dim(const int dimID, const double h[3], const double L[3], const BoundaryType mybc[2], const int sign, const bool isGreen);
    ~FFTW_plan_dim();

    void init(const int size[3], const bool isComplex);

    void allocate_plan(const Topology* topo, double* data);
    void execute_plan(const Topology *topo, double* data) const;

    /**
     * @name Getters - return the value
     * 
     */
    /**@{ */
    inline bool ignoreMode() const { return _ignoreMode; }
    inline bool isSpectral() const { return _isSpectral; }
    inline bool isr2c() const { return _isr2c; }
    // inline bool   isr2c_doneByFFT() const { return _isr2c && (!_isSpectral || !_isGreen) ; } //bug
    inline bool   isr2c_doneByFFT() const { return _isr2c && (!_isSpectral); }
    inline int    dimID() const { return _dimID; }
    inline int    imult() const { return _imult; }
    inline int    shiftgreen() const { return _shiftgreen; }
    inline int    type() const { return _type; }
    inline double symstart() const { return _symstart; }
    inline double normfact() const { return _normfact; }
    inline double volfact() const { return _volfact; }
    inline double kfact() const { return _kfact; }
    inline double koffset() const { return _koffset; }

    inline void get_outsize(int* size) const { size[_dimID] = _n_out; };
    inline void get_fieldstart(int* start) const { start[_dimID] = _fieldstart; };
    inline void get_isNowComplex(bool* isComplex) const { (*isComplex) = (*isComplex) || _isr2c; };
    /**@} */

    void disp();

   protected:
    /**
     * @name Initialization
     */
    /**@{ */
    void _init_real2real(const int size[3], bool isComplex);
    void _init_mixunbounded(const int size[3], bool isComplex);
    void _init_periodic(const int size[3], bool isComplex);
    void _init_unbounded(const int size[3], bool isComplex);
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