/**
 * @file ji0.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef _JI0_H
#define _JI0_H

#include <cmath>

#include "backend.hpp"

#define C_GAMMA 0.5772156649015328606
#define N_KEPT 50

// long double isn't supported in device code
#if M_HAVE_GPU
using ldouble_t = double;
#else
using ldouble_t = long double;
#endif

/**
 * @brief the coefficients 1/(n!)^2 for n = 0..N_KEPT, computed at compile time
 *
 * A constant expression is the only way to fill a constant-memory table on the GPU backends
 * (device code cannot write constant memory, and a kernel doing so is silently ignored).
 */
struct Ji0Table {
    ldouble_t inv_f_sqr[N_KEPT + 1];
    constexpr Ji0Table() : inv_f_sqr() {
        ldouble_t fact = 1;
        inv_f_sqr[0]   = 1;
        for (int i = 1; i <= N_KEPT; i++) {
            fact *= (ldouble_t)i;
            inv_f_sqr[i] = 1 / (fact * fact);
        }
    }
};
m_constant constexpr Ji0Table ji0_table_ = Ji0Table();

/**
 * @brief Numerical approximation of the \f$ \int_0^x  \frac{(1-J0(u))}{u}  du \f$ (for 0 <= x <= ~30).
 *
 */
m_device double Ji0c(double x) {
    ldouble_t val = 0.0;

    for (int n = N_KEPT; n > 0; n--) {
        val -= pow((ldouble_t)-(.25 * x * x), n) * ji0_table_.inv_f_sqr[n] / n * .5;
    }
    return (double)val;
}

/**
 * @brief Numerical approximation of the Bessel-integral function of order zero (for 0 <= x <= ~30).
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
m_device double Ji0(double x) {
    return -Ji0c(x) + log(x / 2) + C_GAMMA;
}

#endif
