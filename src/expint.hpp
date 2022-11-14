/**
 * @file exphint.hpp
 * @copyright Copyright © UCLouvain 2020
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef _H_EXPONENTIAL_INTEGRAL
#define _H_EXPONENTIAL_INTEGRAL

/**********************************************************************/
/*                                                                    */
/*                   double ExponentialIntegral()                     */
/*                                                                    */
/**********************************************************************/
/*                                                                    */
/*  DESCRIPTION:                                                      */
/*  Calculation of the exponential integral E (x) using Chebychev     */
/*                                           1                        */
/*  polynomial expansion.                                             */
/*                                                                    */
/**********************************************************************/

#include <math.h>
#include <iostream>
#include <limits>
#include <stdexcept>

/**********************************************************************/
static double c_gamma = 0.5772156649015328606;
static double expint1(double x);
static double expint2(double x);
static double expint(const int n, const double x);
/**********************************************************************/

static double expint_ei(double x) 
{
    double value;
 
    if (x >= -4. && x <= 4.)
        value = expint1(x);
    else if (x > 4.)
        value = expint2(x);
    else {
        value = 0.;
    }

    return (value);
}

/**********************************************************************/
/*                                                                    */
/*                      double expint1()                              */
/*                                                                    */
/**********************************************************************/
/*                                                                    */
/*  DESCRIPTION:                                                      */
/*  Calculation of the exponential integral  for -4<= x <= 4 using an */
/*  expansion in terms of Chebyshev polynomials.                      */
/*                                                                    */
/**********************************************************************/

static double expint1(double x)
{
    static int MAX = 23; /* The number of coefficients in a[].   */

    static double a[23] = {7.8737715392882774,
                           -8.0314874286705335,
                           3.8797325768522250,
                           -1.6042971072992259,
                           0.5630905453891458,
                           -0.1704423017433357,
                           0.0452099390015415,
                           -0.0106538986439085,
                           0.0022562638123478,
                           -0.0004335700473221,
                           0.0000762166811878,
                           -0.0000123417443064,
                           0.0000018519745698,
                           -0.0000002588698662,
                           0.0000000338604319,
                           -0.0000000041611418,
                           0.0000000004821606,
                           -0.0000000000528465,
                           0.0000000000054945,
                           -0.0000000000005433,
                           0.0000000000000512,
                           -0.0000000000000046,
                           0.0000000000000004};

    int    k;
    double arg, t, value, b0, b1, b2;

    arg = .25 * x; /* Argument in Chebyshev expansion is x/4. */
    t   = 2. * arg;

    b2 = 0.;
    b1 = 0.;
    b0 = a[MAX - 1];

    for (k = MAX - 2; k >= 0; k--) {
        b2 = b1;
        b1 = b0;
        b0 = t * b1 - b2 + a[k];
    }

    value = .5 * (b0 - b2);

    value += log(fabs(x));

    return (-value);

}

/**********************************************************************/
/*                                                                    */
/*                      double expint2()                              */
/*                                                                    */
/**********************************************************************/
/*                                                                    */
/*  DESCRIPTION:                                                      */
/*  Calculation of the exponential integral for x >= 4 using an expan-*/
/*  sionin terms of Chebyshev polynomials.                            */
/*                                                                    */
/**********************************************************************/

static double expint2(double x)
{
    static int MAX = 23; /* The number of coefficients in a[].   */

    static double a[23] = {0.2155283776715125,
                           0.1028106215227030,
                           -0.0045526707131788,
                           0.0003571613122851,
                           -0.0000379341616932,
                           0.0000049143944914,
                           -0.0000007355024922,
                           0.0000001230603606,
                           -0.0000000225236907,
                           0.0000000044412375,
                           -0.0000000009328509,
                           0.0000000002069297,
                           -0.0000000000481502,
                           0.0000000000116891,
                           -0.0000000000029474,
                           0.0000000000007691,
                           -0.0000000000002070,
                           0.0000000000000573,
                           -0.0000000000000163,
                           0.0000000000000047,
                           -0.0000000000000014,
                           0.0000000000000004,
                           -0.0000000000000001};

    int    k;
    double arg, t, value, b0, b1, b2;

    arg = 4. / x; /* Argument in the Chebyshev expansion.       */
    t   = 2. * (2. * arg - 1.);

    b2 = 0.;
    b1 = 0.;
    b0 = a[MAX - 1];

    for (k = MAX - 2; k >= 0; k--) {
        b2 = b1;
        b1 = b0;
        b0 = t * b1 - b2 + a[k];
    }

    value = .5 * (b0 - b2);

    value *= exp(-x);

    return (value);

}

/**********************************************************************/
/*                                                                    */
/*                      double expint()                              */
/*                                                                    */
/**********************************************************************/
// Utility for calculating the exponential integral of order n  from -\infty to x. 
// Evaluates the exponential integral E_n(x)
// This snippets of code is taken from the book Numerical Recipies, third.
//
// Here MAXIT is the maximum allowed number of iterations; 
// c_gamma is Euler’s constant \gamma; 
// EPS is the desired relative error, not smaller than the machine precision; 
// BIG is a number near the largest representable floating-point number.
/**********************************************************************/
static int    MAXIT = 1000;
static double EPS   = std::numeric_limits<double>::epsilon();
static double BIG   = std::numeric_limits<double>::max()* EPS;

static double expint(const int n, const double x) {
    // ------------------------------------------------------------------------
    int    i, ii, nm1 = n - 1;
    double a, b, c, d, del, fact, h, psi, ans;
    if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1))) {
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
        ans  = (nm1 != 0) ? (1.0 / nm1) : (-log(x) - c_gamma);  // Set first term.
        fact = 1.0;
        for (i = 1; i <= MAXIT; i++) {
            fact *= -x / i;
            if (i != nm1) {
                del = -fact / (i - nm1);
            } else {
                psi = -c_gamma;  // Compute \psi
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


#endif