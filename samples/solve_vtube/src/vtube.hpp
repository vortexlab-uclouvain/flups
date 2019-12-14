/**
 * @file vtube.hpp
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

#ifndef VTUBE_HPP
#define VTUBE_HPP

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <cstring>

#include "mpi.h"
#include "flups.h"


static const double c_1opi     = 1.0 / (1.0 * M_PI);
static const double c_1o2pi    = 1.0 / (2.0 * M_PI);
static const double c_1o4pi    = 1.0 / (4.0 * M_PI);
static const double c_1osqrtpi = 1.0 / sqrt(M_PI);
static const double c_1o2      = 1.0 / 2.0;
static const double c_1o4      = 1.0 / 4.0;
static const double c_7o4      = 7. / 4.;
static const double c_19o8     = 19. / 8;
static const double c_2o3      = 2. / 3;
static const double c_1o24     = 1. / 24;
static const double c_1osqrt2 = 1.0 / M_SQRT2;
static const double c_2pi = 2.0 * M_PI;

struct DomainDescr {
    double             xcntr          = 0.5;
    double             ycntr          = 0.5;
    double             zcntr          = 0.5;
    double             xsign          = 0.0;
    double             ysign          = 0.0;
    double             zsign          = 0.0;
    int                nglob[3]       = {64, 64, 64};
    int                nproc[3]       = {1, 2, 2};
    double             L[3]           = {1.0, 1.0, 1.0};
    FLUPS_BoundaryType mybcv[3][2][3] = {{{UNB, UNB, UNB}, {UNB, UNB, UNB}}, {{UNB, UNB, UNB}, {UNB, UNB, UNB}}, {{UNB, UNB, UNB}, {UNB, UNB, UNB}}};
};

/**
 * @name validation of the solver using a gaussian blob
 * 
 */
/**@{ */
void vtube(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int nSolve, int type);
/**@} */

#endif