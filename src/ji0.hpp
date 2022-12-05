/**
 * @file ji0.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef _JI0_H
#define _JI0_H

#include <cmath>


#define C_GAMMA 0.5772156649015328606
#define N_KEPT 50

static long double inv_f_sqr[N_KEPT+1];


/**
 * @brief static precomputation of required coefficients.
 * 
 */
static void init_Ji0(){
    inv_f_sqr[0] = (long double) 1;
    for (int i = 1;i<=N_KEPT;i++){
        inv_f_sqr[i]     = inv_f_sqr[i - 1] * (long double) i;
        inv_f_sqr[i - 1] = 1 / (inv_f_sqr[i - 1] * inv_f_sqr[i - 1]);  //square the previous value
    }
    inv_f_sqr[N_KEPT] = 1 / (inv_f_sqr[N_KEPT] * inv_f_sqr[N_KEPT]);   
}

/**
 * @brief Numerical approximation of the \f$ \int_0^x  \frac{(1-J0(u))}{u}  du \f$ (for 0 <= x <= ~30).
 * 
 * @ref init_Ji0 must be called before.
 */
static inline double Ji0c(double x){
    long double val = 0.0;

    for(int n=N_KEPT; n > 0; n--){
        val -= pow( (long double) -(.25 * x * x), n) * inv_f_sqr[n] / n * .5;
    }
    return (double) val;
}

/**
 * @brief Numerical approximation of the Bessel-integral function of order zero (for 0 <= x <= ~30).
 * 
 * @ref init_Ji0 must be called before.
 * 
 * Meaningful documentation: 
 * - [1] P. Humbert, "Bessel-integral function", Philosophical Magazine, 8 (1929), pp.861-898, (887)
 * - [2] Y.L. Luke, "Bessel functions and their integrals", in "Mathematical Functions and their Approximations", Academic Press, Inc., 1975, pp.311-412.
 *
 * From [1], the definition of Ji0 is
 *    \f[
 *      Ji0(x) = - \int_x^\infty t^{-1} J_0(t) dt
 *    \f]
 * We here use an alternative representation from [1,eq.3]:
 *    \f[
 *      Ji0(x) = \gamma + \log(x/2) - \int_0^x  \frac{(1-J0(u))}{u}  du
 *    \f]
 * and we define the last integral term as "Ji0c". At the bottom of 
 * the same page in [1], Ji0c is expressed as an infinite sum of 
 * polynomials. This is the formulation that we use here.
 * The sum is truncated at n=50, and the result is thus only valid
 * for (0 <= x <= ~30). 
 * 
 * A more general evaluation of Ji0 can be performed using Meijer-G functions.
 */
static inline double Ji0(double x){
    return -Ji0c(x) + log(x/2) + C_GAMMA;
}

#endif
