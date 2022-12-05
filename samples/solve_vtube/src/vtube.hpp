/**
 * @file vtube.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
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

/**********************************************************************/
/*                      Needed constant                               */
/**********************************************************************/
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

/**********************************************************************/
/*                      double expint()                              */
/**********************************************************************/
template <int P>
double gexpint(const double z) {
    // DLMF 8.19.12 (https://dlmf.nist.gov/8.19#E12)
    // or Abramowitz and stegun 5.1.14
    return (std::exp(-z) - z * gexpint<P - 1>(z)) / (P - 1.0);
}
template <> 
double gexpint<1>(const double z);
//{
//    // for real values E1(x) = -Ei(-x):
//    // according to
//    // https://en.cppreference.com/w/cpp/numeric/special_functions/expint and
//    // https://en.wikipedia.org/wiki/Exponential_integral
//    return (-std::expint(-z));
//}

/**********************************************************************/
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
    double             rc             = 0.1;
    double             sigma          = 0.05;
    bool               compact        = false;
    double             rad            = 0.25;
};

/**********************************************************************/
/**
 * @name validation of the solver using a gaussian blob
 *
 */
/**@{ */
void vtube(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int nSolve, int type, FLUPS_DiffType order, int vdir);
/**@} */
/**********************************************************************/
#endif
