/**
 * @file ji0.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief Bessel-integral function of order 0
 * @version
 * 
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

#ifndef _JI0_H
#define _JI0_H

#include <cmath>

/*
 * Numerical approximation of the Bessel-integral function of order zero.
 *
 * Meaningful documentation: 
 * [1] P. Humbert, "Bessel-integral function", Philosophical Magazine, 8 (1929), pp.861-898, (887)
 * [2] Y.L. Luke, "Bessel functions and their integrals", in "Mathematical Functions and their Approximations", Academic Press, Inc., 1975, pp.311-412.
 *
 * From [1], the definition of Ji0 is
 *    Ji0(x) = - \int_x^\infty t^{-1} J_0(t) dt
 * We here use an alternative representation from [1,eq.3]:
 *    Ji0(x) = \gamma + log(x/2) - \int_0^x  (1-J0(u))/u  du
 * and we define the last integral term as "Ji0c". At the bottom of 
 * the same page in [1], Ji0c is expressed as an infinite sum of 
 * polynomials. This is the formulation that we use here.
 * The sum is truncated at n=50, and the result is thus only valid
 * for (0 <= x <= ~30). 
 */

#define C_GAMMA 0.5772156649015328606
#define N_KEPT 50

static long double inv_f_sqr[N_KEPT];


//precompute and store sqr of the factorial
static void init_Ji0(){
    inv_f_sqr[0] = (long double) 1;
    for (int i = 1;i<=N_KEPT;i++){
        inv_f_sqr[i]     = inv_f_sqr[i - 1] * (long double) i;
        inv_f_sqr[i - 1] = 1 / (inv_f_sqr[i - 1] * inv_f_sqr[i - 1]);  //square the previous value
    }
    inv_f_sqr[N_KEPT] = 1 / (inv_f_sqr[N_KEPT] * inv_f_sqr[N_KEPT]);   
}

static inline double Ji0c(double x){
    long double val = 0.0;

    for(int n=N_KEPT; n > 0; n--){
        val -= pow( (long double) -(.25 * x * x), n) * inv_f_sqr[n] / n * .5;
        // out = out - (-1)^i * (.5*x).^(2*i) / f.^2 /i *.5;
    }
    return (double) val;
}

static inline double Ji0(double x){
    return -Ji0c(x) + log(x/2) + C_GAMMA;
}

#endif