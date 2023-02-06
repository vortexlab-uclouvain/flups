/**
 * @file green_kernels.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief defines the 3D Green functions kernels
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

#include "defines.hpp"
#include "expint.hpp"
#include "si.hpp"
#include "ji0.hpp"
#include "lgf_3unb0spe.hpp"
#include "lgf_1unb2spe.hpp"

using LGFOneUnbounded::lgf2_symbol;
using LGFOneUnbounded::lgf4_symbol;
using LGFOneUnbounded::lgf6_symbol;

/**
 * @name 3 directions unbounded - 0 direction spectral
 * 
 * @{
 */
// ----------------------------------------------------------- 3D - KERNELS ----------------------------------------------------------
//notice that these function will likely not be inlined as we have a pointer to them...
static inline double hej_2_3unb0spe_(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    return -c_1o4pi / r * (erf(r / eps * c_1osqrt2));
}
static inline double hej_4_3unb0spe_(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return -c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (rho)*exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double hej_6_3unb0spe_(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return -c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (c_7o4 * rho - c_1o4 * pow(rho, 3)) * exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double hej_8_3unb0spe_(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return -c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (c_19o8 * rho - c_2o3 * pow(rho, 3) + c_1o24 * pow(rho, 5)) * exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double hej_10_3unb0spe_(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return -c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (c_187o64 * rho - c_233o192 * pow(rho, 3) + c_29o192 * pow(rho, 5) - c_1o192 * pow(rho, 7)) * exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double hej_0_3unb0spe_(const void* params,const double* data) {
    double r       = ((double*)params)[0];
    double sigma   = ((double*)params)[1];
    double c_1osig = 1.0 / sigma;
    double rho     = r * c_1osig;
    return -c_1o2pi2 * c_1osig * Si(rho)/rho;
}
static inline double chat_2_3unb0spe_(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    return -c_1o4pi / r ;
}
static inline double lgf_2_3unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];
    double h  = ((double*)params)[8];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N + iz * N * N];
    } else {  // if not, we use the extrapolation
        const double rho  = sqrt(ix * ix + iy * iy + iz * iz);
        lgf_2_3unb0spe_expansion(&green, ix, iy, iz, rho);
    }
    return -green/(h);
}
static inline double lgf_4_3unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];
    double h  = ((double*)params)[8];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N + iz * N * N];
    } else {  // if not, we use the extrapolation
        const double rho  = sqrt(ix * ix + iy * iy + iz * iz);
        lgf_4_3unb0spe_expansion(&green, ix, iy, iz, rho);
    }
    return -green/(h);
}
static inline double lgf_6_3unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];
    double h  = ((double*)params)[8];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N + iz * N * N];
    } else {  // if not, we use the extrapolation
        const double rho  = sqrt(ix * ix + iy * iy + iz * iz);
        lgf_6_3unb0spe_expansion(&green, ix, iy, iz, rho);
    }
    return -green/(h);
}
static inline double lgf_2_2unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N];

    } else {  // if not, we use the extrapolation
        const double rho     = sqrt(ix * ix + iy * iy);
        lgf_2_2unb0spe_expansion(&green, ix, iy, rho);
    }
    return -green;
}
static inline double lgf_4_2unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N];
    } else {  // if not, we use the extrapolation
        const double rho     = sqrt(ix * ix + iy * iy);
        lgf_4_2unb0spe_expansion(&green, ix, iy, rho);
    }
    return -green;
}
static inline double lgf_6_2unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N];
    } else {  // if not, we use the extrapolation
        const double rho     = sqrt(ix * ix + iy * iy);
        lgf_6_2unb0spe_expansion(&green, ix, iy, rho);
    }
    return -green;
}
/**@} */


/**
 * @name 2 directions unbounded - 1 direction spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
static inline double hej_2_2unb1spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params)[0];
    const double sig = ((double*)params)[2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    // return -c_1o2pi * (log(r) - exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2)); //mistaken coefs in [Spietz2018]
    // return -c_1o2pi * (log(r) + .5 * expint_ei(rho2 / 2.0));
    return c_1o2pi * (log(r) + 0.5 * expint_ei(rho2 * 0.5));
    // return -c_1o2pi * (.5*log(rho*.5) + .5 * expint_ei(rho2 / 2));
}
static inline double hej_2_2unb1spe_r0_(const void* params,const double* data) {
    const double sig = ((double*)params)[2];
    return -c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig));
}

static inline double hej_4_2unb1spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho  = r/sig;
    const double rho2 = rho*rho;
    // return -c_1o2pi * (log(r) - (1 - .5 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2)); //mistaken coefs in [Spietz2018]
    // return -c_1o2pi * (log(r) - exp(-rho2 / 2.0) + .5 * expint_ei(rho2 / 2.0));
    return c_1o2pi * (log(r) - 0.5 * exp(-rho2 * 0.5) + 0.5 * expint_ei(rho2 * 0.5));
    // return -c_1o2pi * (.5*log(rho2*.5) - exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2));
}
static inline double hej_4_2unb1spe_r0_(const void* params,const double* data) {
    const double sig = ((double*)params)[2];

    return -c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig) + .5);
}

static inline double hej_6_2unb1spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    // return -c_1o2pi * (log(r) - (1 - rho2 + .125 * rho2 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2)); //mistaken coefs in [Spietz2018]
    return c_1o2pi * (log(r) - (0.75 - 0.125 * rho2) * exp(-rho2 * 0.5) + 0.5 * expint_ei(rho2 * 0.5));
    // return -c_1o2pi * (.5*log(rho2*.5) - (.75 - .125 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2));
}
static inline double hej_6_2unb1spe_r0_(const void* params,const double* data) {
    const double sig = ((double*)params)[2];

    return -c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig) + .75);
}

static inline double hej_8_2unb1spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    return c_1o2pi * (log(r) - (c_11o12 - c_7o24 * rho2 + c_1o48 * rho2*rho2) * exp(-rho2 * 0.5) + 0.5 * expint_ei(rho2 * 0.5));
}
static inline double hej_8_2unb1spe_r0_(const void* params,const double* data) {
    const double sig = ((double*)params)[2];

    return -c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig) + c_11o12);
}

static inline double hej_10_2unb1spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    return c_1o2pi * (log(r) - (c_25o24 - c_23o48 * rho2 + c_13o192 * rho2*rho2 - c_1o384 * rho2*rho2*rho2) * exp(-rho2 * 0.5) + 0.5 * expint_ei(rho2 * 0.5));
}
static inline double hej_10_2unb1spe_r0_(const void* params,const double* data) {
    const double sig = ((double*)params)[2];

    return -c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig) + c_25o24);
}
static inline double hej_0_2unb1spe_k0_(const void* params,const double* data) {
    //works as well for r0
    const double r       = ((double*)params)[0];
    const double sigma   = ((double*)params)[2]; // h/pi
    const double rho     = r / sigma;

    //Switching between the approximate G with sharp cutoff and the
    // analytical singular solution. 
    // CAUTION: the switch is done arbitrarily and abruplty, with no
    //          mathematical reason. Hence, the result is likely inaccurate.
    const double tmp = c_1o2pi * (Ji0c(rho) + log(2 * sigma) - c_gamma);
    return (rho<30.0) ? tmp : c_1o2pi * log(r); 
}

static inline double zero_(const void* params,const double* data) {   
    return 0.0;
}

static inline double chat_2_2unb1spe_(const void* params,const double* data) {
    const double r      = ((double*)params) [0];
    const double k      = ((double*)params) [1];

    return -c_1o2pi * besselk0(fabs(k) * r);
}
static inline double chat_2_2unb1spe_r0_(const void* params,const double* data) {
    const double k      = ((double*)params) [1];
    const double r_eq2D = ((double*)params) [3];

    return -(1.0 - k * r_eq2D * besselk1(k * r_eq2D)) * c_1opi / ((k * r_eq2D) * (k * r_eq2D));
}
static inline double chat_2_2unb1spe_k0_(const void* params,const double* data) {
    const double r      = ((double*)params) [0];
    // const double sig = ((double*)params)[2];
    
    return  c_1o2pi * log(r) ; //caution: mistake on the sign in [Chatelain2010]
}

/**@} */

/**
 * @name 1 direction unbounded - 2 directions spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
static inline double hej_2_1unb2spe_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;

    const double subfun = s * rho > 100. ? 0.0 : ((1.0 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1.0 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return - .25 * sig / s * subfun ;
}
static inline double hej_2_1unb2spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rosqrt2 = r*c_1osqrt2;
    // return -.5* (r * erf(rosqrt2/sig) + (exp(-r*r/(2*sig*sig)) - 1.)*sig*M_SQRT2*c_1osqrtpi) ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
    return 0.5 * r * erf(rosqrt2/sig) - (1.-exp(-rho*rho*.5)) *sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double hej_4_1unb2spe_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return - 0.25 * sig / s * subfun \
           - sig * M_SQRT2 * c_1osqrtpi * .25 * exp(-.5 * (s * s + rho * rho));
}
static inline double hej_4_1unb2spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rosqrt2 = r*c_1osqrt2;
    return 0.5 * r * erf(rosqrt2/sig) - (1.-exp(-rho*rho*.5)) *.5*sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double hej_6_1unb2spe_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return - 0.25 * sig / s * subfun \
           - sig * M_SQRT2 * c_1osqrtpi * (c_5o16 + c_1o16 * (s * s - rho * rho)) * exp(-.5 * (s * s + rho * rho));
}
static inline double hej_6_1unb2spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rosqrt2 = r*c_1osqrt2;
    return 0.5 * r * erf(rosqrt2/sig) - (3.-exp(-rho*rho*.5) * (rho*rho+3.) ) *.125*sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double hej_8_1unb2spe_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho  = r/sig;
    const double s    = k*sig;
    const double s2   = s*s;
    const double rho2 = rho*rho;
    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return - 0.25 * sig / s * subfun \
           - sig * M_SQRT2 * c_1osqrtpi * (c_11o32 + c_1o12 * s2 - 0.125 * rho2 - c_1o48 * s2*rho2 + c_1o96 * (s2*s2+rho2*rho2) ) * exp(-.5 * (s * s + rho2));
}
static inline double hej_8_1unb2spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    const double rosqrt2 = r*c_1osqrt2;
    return 0.5 * r * erf(rosqrt2/sig) - (c_5o16 - exp(-rho*rho*.5) * (-c_1o48 * rho2*rho2 + .25 * rho2 + c_5o16) ) *sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double hej_10_1unb2spe_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho  = r/sig;
    const double s    = k*sig;
    const double s2   = s*s;
    const double rho2 = rho*rho;
    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return - 0.25 * sig / s * subfun \
           - sig * M_SQRT2 * c_1osqrtpi * (c_93o256 + c_73o768 * s2 - c_47o256 * rho2 - c_17o384 * s2*rho2 + c_11o768 * s2*s2 + c_23o768 * rho2*rho2 + c_1o256 * (s2*rho2*rho2 - s2*s2*rho2) + c_1o768 * (s2*s2*s2 - rho2*rho2*rho2) ) * exp(-.5 * (s * s + rho2));
}
static inline double hej_10_1unb2spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    const double rosqrt2 = r*c_1osqrt2;
    return 0.5 * r * erf(rosqrt2/sig) - (c_35o128 - exp(-rho*rho*.5) * (c_1o384 * rho2*rho2*rho2 - c_23o384 * rho2*rho2 + c_47o128 * rho2 + c_35o128) ) *sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double chat_2_1unb2spe_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];

    return -0.5 * exp(-k * r) / k;
}
static inline double chat_2_1unb2spe_k0_(const void* params,const double* data) {
    const double r   = ((double*)params) [0];

    return 0.5 * fabs(r);
}
static inline double lgf_2_1unb2spe_(const void* params,const double* data) {
    const double r  = ((double*)params) [0];
    const double kx = ((double*)params) [3];
    const double ky = ((double*)params) [4];
    const double kz = ((double*)params) [5];
    const double h  = ((double*)params) [6];

    // one of (kx, ky, kz) is zero, so one of the symbols will be as well
    const double c = lgf2_symbol(kx * h) + lgf2_symbol(ky * h) + lgf2_symbol(kz * h);
    const int n = (int)round(abs(r) / h);
    return -h * LGFOneUnbounded::lgf2(n, c);
}
static inline double lgf_4_1unb2spe_(const void* params,const double* data) {
    const double r  = ((double*)params) [0];
    const double kx = ((double*)params) [3];
    const double ky = ((double*)params) [4];
    const double kz = ((double*)params) [5];
    const double h  = ((double*)params) [6];

    const double c = lgf4_symbol(kx * h) + lgf4_symbol(ky * h) + lgf4_symbol(kz * h);
    const int n = (int)round(abs(r) / h);
    return -h * LGFOneUnbounded::lgf4(n, c);
}
static inline double lgf_6_1unb2spe_(const void* params,const double* data) {
    const double r  = ((double*)params) [0];
    const double kx = ((double*)params) [3];
    const double ky = ((double*)params) [4];
    const double kz = ((double*)params) [5];
    const double h  = ((double*)params) [6];

    const double c = lgf6_symbol(kx * h) + lgf6_symbol(ky * h) + lgf6_symbol(kz * h);
    const int n = (int)round(abs(r) / h);
    return -h * LGFOneUnbounded::lgf6(n, c);
}
/**@} */


/**
 * @name 3 directions spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
static inline double hej_2_0unb3spe_(const void* params,const double* data) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    return - exp(-ssqr * 0.5) / (ksqr); 
}
static inline double hej_4_0unb3spe_(const void* params,const double* data) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    return - (1.0 + ssqr * 0.5) * exp(-ssqr * 0.5) / (ksqr);
}
static inline double hej_6_0unb3spe_(const void* params,const double* data) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    return - (1.0 + ssqr * 0.5 + ssqr * ssqr * 0.125 ) * exp(-ssqr * 0.5) / (ksqr);
}
static inline double hej_8_0unb3spe_(const void* params,const double* data) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    const double ssqr2 = ssqr * ssqr;
    return - (1.0 + ssqr * 0.5 + ssqr2 * 0.125 + ssqr2 * ssqr * c_1o48) * exp(-ssqr * 0.5) / (ksqr);
}
static inline double hej_10_0unb3spe_(const void* params,const double* data) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    const double ssqr2 = ssqr * ssqr;
    return - (1.0 + ssqr * 0.5 + ssqr2 * 0.125 + ssqr2 * ssqr * c_1o48 + ssqr2 * ssqr2 * c_1o384) * exp(-ssqr * 0.5) / (ksqr);
}

static inline double chat_2_0unb3spe_(const void* params,const double* data) {
    const double ksqr   = ((double*)params) [0];

    return - 1.0 / ksqr;
}
static inline double lgf_2_0unb3spe_(const void* params, const double* data) {
    const double k[3] = {((double*)params)[2], ((double*)params)[3], ((double*)params)[4]};
    const double h  = ((double*)params)[5];
    return - h * h / (lgf2_symbol(k[0] * h) + lgf2_symbol(k[1] * h) + lgf2_symbol(k[2] * h));
}
static inline double lgf_4_0unb3spe_(const void* params, const double* data) {
    const double k[3] = {((double*)params)[2], ((double*)params)[3], ((double*)params)[4]};
    const double h  = ((double*)params)[5];
    return - h * h / (lgf4_symbol(k[0] * h) + lgf4_symbol(k[1] * h) + lgf4_symbol(k[2] * h));
}
static inline double lgf_6_0unb3spe_(const void* params, const double* data) {
    const double k[3] = {((double*)params)[2], ((double*)params)[3], ((double*)params)[4]};
    const double h  = ((double*)params)[5];
    return - h * h / (lgf6_symbol(k[0] * h) + lgf6_symbol(k[1] * h) + lgf6_symbol(k[2] * h));
}
/**@} */