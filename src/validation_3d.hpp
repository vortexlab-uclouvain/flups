/**
 * @file validation_3d.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
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

#ifndef VALIDATION_3D_HPP
#define VALIDATION_3D_HPP

#include "Solver.hpp"
#include "defines.hpp"
#include "expint.hpp"

#include "hdf5_io.hpp"
#include "mpi.h"

#define MANUFACTURED_SOLUTION


struct DomainDescr {
    int          nglob[3]   = {64, 64, 64};
    int          nproc[3]   = {1, 2, 1};
    double       L[3]       = {1.0, 1.0, 1.0};
    FLUPS::BoundaryType mybc[3][2] = {{FLUPS::UNB, FLUPS::UNB}, {FLUPS::UNB, FLUPS::UNB}, {FLUPS::UNB, FLUPS::UNB}};
    // BoundaryType mybc[3][2] = {{UNB, UNB}, {UNB, UNB}, {UNB, UNB}};
};

/**
 * @name validation of the solver using a gaussian blob
 * 
 */
/**@{ */
void validation_3d(const DomainDescr myCase, const FLUPS::SolverType type, const FLUPS::GreenType typeGreen);
void validation_3d(const DomainDescr myCase, const FLUPS::SolverType type, const FLUPS::GreenType typeGreen, const int nSolve);
/**@} */


#ifdef MANUFACTURED_SOLUTION


struct manuParams {
    double       freq    = 1; //an integer or 0.5
    double       sign[2] = {0., 0.};
    double       sigma   = 0.15;
    double       center  = 0.5;
};

typedef double (*manuF)(const double, const double, const manuParams);


static inline double fOddOdd(const double x, const double L, const manuParams params) {
    return sin((c_2pi / L * params.freq) * x);
}
static inline double d2dx2_fOddOdd(const double x, const double L, const manuParams params) {
    return -(c_2pi / L * params.freq) * (c_2pi / L * params.freq) * sin((c_2pi / L * params.freq) * x);
}

static inline double fEvenEven(const double x, const double L, const manuParams params) {
    return cos((M_PI / L * params.freq) * x);
}
static inline double d2dx2_fEvenEven(const double x, const double L, const manuParams params) {
    return -(M_PI / L * params.freq) * (M_PI / L * params.freq) * cos((M_PI / L * params.freq) * x);
}

static inline double fOddEven(const double x, const double L, const manuParams params) {
    return sin((M_PI / L * (params.freq+.5)) * x);
}
static inline double d2dx2_fOddEven(const double x, const double L, const manuParams params) {
    return -(M_PI / L * (params.freq+.5)) * (M_PI / L * (params.freq+.5)) * sin((M_PI / L * (params.freq+.5)) * x);
}

static inline double fEvenOdd(const double x, const double L, const manuParams params) {
    return cos((M_PI / L * (params.freq+.5)) * x);
}
static inline double d2dx2_fEvenOdd(const double x, const double L, const manuParams params) {
    return -(M_PI / L * (params.freq+.5)) * (M_PI / L * (params.freq+.5)) * cos((M_PI / L * (params.freq+.5)) * x);
}

static inline double fUnb(const double x, const double L, const manuParams params) {
    const double x0 = (x -       params.center  * L) / params.sigma;
    const double x1 = (x +       params.center  * L) / params.sigma;
    const double x2 = (x - (2. - params.center) * L) / params.sigma;
    return                   exp(-x0*x0) + \
            params.sign[0] * exp(-x1*x1) + \
            params.sign[1] * exp(-x2*x2) ;
}
static inline double d2dx2_fUnb(const double x, const double L, const manuParams params) {
    const double x0 = (x -       params.center  * L) / params.sigma;
    const double x1 = (x +       params.center  * L) / params.sigma;
    const double x2 = (x - (2. - params.center) * L) / params.sigma;
    return                   -2. / (params.sigma * params.sigma) * exp(-x0*x0 ) * (1. - 2. * ( x0 * x0)) + \
            params.sign[0] * -2. / (params.sigma * params.sigma) * exp(-x1*x1 ) * (1. - 2. * ( x1 * x1)) + \
            params.sign[1] * -2. / (params.sigma * params.sigma) * exp(-x2*x2 ) * (1. - 2. * ( x2 * x2)) ;
}

#endif


#endif