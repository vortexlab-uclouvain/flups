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
// Utility for calculating the exponential integral of order n  from -\infty to x. 
// Evaluates the exponential integral E_n(x)
// This snippets of code is taken from the book Numerical Recipies, third.
//
// Here MAXIT is the maximum allowed number of iterations; 
// EULER is Euler’s constant \gamma; 
// EPS is the desired relative error, not smaller than the machine precision; 
// BIG is a number near the largest representable floating-point number.
/**********************************************************************/
static int    MAXIT = 100;
static double EULER = 0.5772156649015328606;
static double EPS   = std::numeric_limits<double>::epsilon();
static double BIG   = std::numeric_limits<double>::max()* EPS;

static double expint(const int n, const double x) {
    // ------------------------------------------------------------------------
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int    i, ii, nm1 = n - 1;
    double a, b, c, d, del, fact, h, psi, ans;
    if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1))) {
        if (rank == 0){ 
            printf(" The parameters are x = %e - n = %d \n", x, n);
        }
        throw std::runtime_error("bad arguments in expint");
    }
    if (n == 0) {
        ans = exp(-x) / x;  // Special case
    } else if (x == 0.0) {
        ans = 1.0 / nm1;
    } else if (x > 1.0) {  // Lentz's algorithm
        b = x + n;
        c = BIG;
        d = 1.0 / b;
        h = d;
        for (i = 1; i <= MAXIT; i++) {
            a   = -i * (nm1 + i);
            b   += 2.0;
            d   = 1.0 / (a * d + b);
            c   = b + a / c;
            del = c * d;
            h   *= del;
            if (abs(del - 1.0) <= EPS) {
                ans = h * exp(-x);
                return ans;
            }
        }
        throw std::runtime_error("continued fraction failed in expint");
    } else {                                              // Evaluate series.
        ans  = (nm1 != 0) ? (1.0 / nm1) : (-log(x) - EULER);  // Set first term.
        fact = 1.0;
        for (i = 1; i <= MAXIT; i++) {
            fact *= -x / i;
            if (i != nm1) {
                del = -fact / (i - nm1);
            } else {
                psi = -EULER;  // Compute \psi
                for (ii = 1; ii <= nm1; ii++) {
                    psi += 1.0 / ii;
                }
                del = fact * (-log(x) + psi);
            }
            ans += del;
            if (abs(del) < abs(ans) * EPS) {
                return ans;
            }
        }
        throw std::runtime_error("series failed in expint");
    }
    // ------------------------------------------------------------------------
    return ans;
}

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