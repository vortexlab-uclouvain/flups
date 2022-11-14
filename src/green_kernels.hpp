/**
 * @file green_kernels.hpp
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

#include "defines.hpp"
#include "expint.hpp"
#include "si.hpp"
#include "ji0.hpp"

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

/**
 * @brief LGF 3D
 * 
 * @param params 
 * @param data 
 * @return double 
 */
static inline double lgf_2_3unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];
    double h  = ((double*)params)[8];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = - data[ix + iy * N + iz * N * N];

    } else {  // if not, we use the extrapolation
        const double rho     = sqrt(ix * ix + iy * iy + iz * iz);
        const double rho_2   = rho * rho;
        const double oorho_6 = 1.0 / std::pow(rho, 6.0);
        const double oorho_7 = 1.0 / std::pow(rho, 7.0);
        // ix
        const double ix_2  = std::pow(ix, 2.0);
        const double ix_4  = std::pow(ix, 4.0);
        const double ix_6  = std::pow(ix, 6.0);
        const double ix_8  = std::pow(ix, 8.0);
        const double ix_10 = std::pow(ix, 10.0);
        const double ix_12 = std::pow(ix, 12.0);
        // iy
        const double iy_2  = std::pow(iy, 2.0);
        const double iy_4  = std::pow(iy, 4.0);
        const double iy_6  = std::pow(iy, 6.0);
        const double iy_8  = std::pow(iy, 8.0);
        const double iy_10 = std::pow(iy, 10.0);
        const double iy_12 = std::pow(iy, 12.0);
        //iz
        const double iz_2  = std::pow(iz, 2.0);
        const double iz_4  = std::pow(iz, 4.0);
        const double iz_6  = std::pow(iz, 6.0);
        const double iz_8  = std::pow(iz, 8.0);
        const double iz_10 = std::pow(iz, 10.0);
        const double iz_12 = std::pow(iz, 12.0);

        green = - c_1o4pi / rho \
            - 1.0/(  16.0 * M_PI) * (ix_4 + iy_4 + iz_4 - 3.0 * (ix_2 * iy_2 + iy_2 * iz_2 + ix_2 * iz_2)) * oorho_7 \
            - 1.0/( 128.0 * M_PI) * (23.0 * (ix_8 + iy_8 + iz_8) - 244.0 * (ix_6 * (iy_2 + iz_2) + iy_6 * (ix_2 + iz_2) + iz_6 * (ix_2 + iy_2)) - 228.0 * ix_2 * iy_2 * iz_2 * rho_2 + 621.0 * (ix_4 * iy_4 + ix_4 * iz_4 + iy_4 * iz_4)) * oorho_7 * oorho_6 \
            - 1.0/(2048.0 * M_PI) * (2588.0 * (ix_12 + iy_12 + iz_12) - 65676.0 * (ix_10 * iy_2 + ix_10 * iz_2 + ix_2 * iy_10 + iy_10 * iz_2 + ix_2 * iz_10 + iy_2 * iz_10) + 426144.0 * (ix_8 * iy_4 + ix_4 * iy_8 + ix_8 * iz_4 + iy_8 * iz_4 + ix_4 * iz_8 + iy_4 * iz_8) - 712884.0 * (ix_6 * iy_6 + iy_6 * iz_6 + ix_6 * iz_6) - 62892.0 * (ix_8 * iy_2 * iz_2 + ix_2 * iy_8 * iz_2 + ix_2 * iy_2 * iz_8) - 297876.0 * (ix_6 * iy_4 * iz_2 + ix_4 * iy_6 * iz_2 + ix_4 * iy_2 * iz_6 + ix_2 * iy_4 * iz_6 + ix_6 * iy_2 * iz_4 + ix_2 * iy_6 * iz_4) + 2507340.0 * ix_4 * iy_4 * iz_4) * oorho_7 * oorho_6 * oorho_6;
    }
    
    return green/(h);
}
/**
 * @brief LGF 2D
 * 
 * @param params 
 * @param data 
 * @return double 
 */
static inline double lgf_2_2unb0spe_(const void* params,const double* data) {
    int    ix = (int)((double*)params)[4];
    int    iy = (int)((double*)params)[5];
    int    iz = (int)((double*)params)[6];
    int    N  = (int)((double*)params)[7];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = - data[ix + iy * N];

    } else {  // if not, we use the extrapolation
        const double rho     = sqrt(ix * ix + iy * iy);
        const double oorho_6 = 1.0 / std::pow(rho, 6.0);
        // const double ix_1     = ix;
        const double ix_2  = std::pow(ix, 2.0);
        const double ix_4  = std::pow(ix, 4.0);
        const double ix_6  = std::pow(ix, 6.0);
        const double ix_8  = std::pow(ix, 8.0);
        const double ix_10 = std::pow(ix, 10.0);
        const double ix_12 = std::pow(ix, 12.0);
        const double ix_14 = std::pow(ix, 14.0);
        const double ix_16 = std::pow(ix, 16.0);
        // const double iy_1     = iy;
        const double iy_2  = std::pow(iy, 2.0);
        const double iy_4  = std::pow(iy, 4.0);
        const double iy_6  = std::pow(iy, 6.0);
        const double iy_8  = std::pow(iy, 8.0);
        const double iy_10 = std::pow(iy, 10.0);
        const double iy_12 = std::pow(iy, 12.0);
        const double iy_14 = std::pow(iy, 14.0);
        const double iy_16 = std::pow(iy, 16.0);

        green =   1.0 / (   2.0 * M_PI) * (log(rho) + c_gamma + log(8.0) * c_1o2)\
                - 1.0 / (  24.0 * M_PI) * (ix_4 - 6.0 * ix_2 * iy_2 + iy_4) * oorho_6\
                - 1.0 / ( 480.0 * M_PI) * (43.0 * (ix_8 + iy_8) - 772.0 * (ix_6 * iy_2 + ix_2 * iy_6) + 1570.0 * ix_4 * iy_4) * oorho_6 * oorho_6\
                - 1.0 / (2016.0 * M_PI) * (609.0 * (ix_12 + iy_12) - 24234.0 * (ix_10 * iy_2 + ix_2 * iy_10) + 109935.0 * (ix_8 * iy_4 + ix_4 * iy_8) - 160524.0 * ix_6 * iy_6) * oorho_6 * oorho_6 * oorho_6\
                - 1.0 / (2880.0 * M_PI) * (63139.0 * (ix_16 + iy_16) - 4467336.0 * (ix_14 * iy_2 + ix_2 * iy_14) + 38334996.0 * (ix_12 * iy_4 + ix_4 * iy_12) - 98512568.0 * (ix_10 * iy_6 + ix_6 * iy_10) + 122747922.0 * ix_8 * iy_8) * oorho_6 * oorho_6 * oorho_6 * oorho_6;
    }
    return green;
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
    const double kx = ((double*)params)[2];
    const double ky = ((double*)params)[3];
    const double kz = ((double*)params)[4];
    const double h  = ((double*)params)[5];

    return - h * h / (4.0 * pow(sin(kx * h / 2.0), 2.0) + 4.0 * pow(sin(ky * h / 2.0), 2.0) + 4.0 * pow(sin(kz * h / 2.0), 2.0));
}
/**@} */