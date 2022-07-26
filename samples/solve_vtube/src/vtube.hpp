/**
 * @file vtube.hpp
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
#include "h3lpr/profiler.hpp"
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
    FLUPS_CenterType   center[3]      = {CELL_CENTER, CELL_CENTER, CELL_CENTER};
};

/**
 * @name validation of the solver using a gaussian blob
 * 
 */
/**@{ */
void vtube(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int nSolve, int type, FLUPS_DiffType order, int vdir);
/**@} */

#endif