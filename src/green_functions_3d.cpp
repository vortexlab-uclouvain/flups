/**
 * @file green_functions_3d.cpp
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

#include "green_functions_3d.hpp"

// **Symmetry computation:**
// 
// We have to take the symmetry around symstart. e.g. in X direction: `symstart[0] - (ix - symstart[0]) = 2 symstart[0] - ix`
// 
// In some cases when we have an R2C transform, it ask for 2 additional doubles.
// The value is meaningless but we would like to avoid segfault and nan's.
// To do so, we use 2 tricks:
// - The `abs` is used to stay on the positivie side and hence avoid negative memory access
// - The `max` is used to prevent the computation of the value in 0, which is never used in the symmetry.
// 
// As an example, the final formula is then ( in the X direction):
// `max( abs(2 symstart[0] - ix) , 1)`

/**
 * @brief generic type for Green kernel, takes a table of parameters that can be used depending on the kernel
 * 
 */
typedef double (*GreenKernel)(const void*);
/**
 * @brief generic type for unbounded Green kernel, alike the @ref GreenKernel but with an additional data array for the LGF case
 * 
 */
typedef double (*GreenKernel_unb)(const void*,const double*);


/**
 * @name 3 directions unbounded - 0 direction spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
//notice that these function will likely not be inlined as we have a pointer to them...
static inline double _hej_2_3unb0spe(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    return c_1o4pi / r * (erf(r / eps * c_1osqrt2));
}
static inline double _hej_4_3unb0spe(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (rho)*exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double _hej_6_3unb0spe(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (c_7o4 * rho - c_1o4 * pow(rho, 3)) * exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double _chat_2_3unb0spe(const void* params,const double* data) {
    double r   = ((double*)params) [0];
    return c_1o4pi / r ;
}
static inline double _lgf_2_3unb0spe(const void* params,const double* data) {
    int    ix = (int)((double*)params)[2];
    int    iy = (int)((double*)params)[3];
    int    iz = (int)((double*)params)[4];
    int    N  = (int)((double*)params)[5];
    double h  = ((double*)params)[6];

    // if the point is close enough, it will be already precomputed
    double green;
    if (ix < N && iy < N && iz < N) {
        green = data[ix + iy * N + iz * N * N];

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

        green = c_1o4pi / rho \
            + 1.0/(  16.0 * M_PI) * (ix_4 + iy_4 + iz_4 - 3.0 * (ix_2 * iy_2 + iy_2 * iz_2 + ix_2 * iz_2)) * oorho_7 \
            + 1.0/( 128.0 * M_PI) * (23.0 * (ix_8 + iy_8 + iz_8) - 244.0 * (ix_6 * (iy_2 + iz_2) + iy_6 * (ix_2 + iz_2) + iz_6 * (ix_2 + iy_2)) - 228.0 * ix_2 * iy_2 * iz_2 * rho_2 + 621.0 * (ix_4 * iy_4 + ix_4 * iz_4 + iy_4 * iz_4)) * oorho_7 * oorho_6 \
            + 1.0/(2048.0 * M_PI) * (2588.0 * (ix_12 + iy_12 + iz_12) - 65676.0 * (ix_10 * iy_2 + ix_10 * iz_2 + ix_2 * iy_10 + iy_10 * iz_2 + ix_2 * iz_10 + iy_2 * iz_10) + 426144.0 * (ix_8 * iy_4 + ix_4 * iy_8 + ix_8 * iz_4 + iy_8 * iz_4 + ix_4 * iz_8 + iy_4 * iz_8) - 712884.0 * (ix_6 * iy_6 + iy_6 * iz_6 + ix_6 * iz_6) - 62892.0 * (ix_8 * iy_2 * iz_2 + ix_2 * iy_8 * iz_2 + ix_2 * iy_2 * iz_8) - 297876.0 * (ix_6 * iy_4 * iz_2 + ix_4 * iy_6 * iz_2 + ix_4 * iy_2 * iz_6 + ix_2 * iy_4 * iz_6 + ix_6 * iy_2 * iz_4 + ix_2 * iy_6 * iz_4) + 2507340.0 * ix_4 * iy_4 * iz_4) * oorho_7 * oorho_6 * oorho_6;
    }
    
    return green/(h);
}

/**
 * @brief Compute the Green function for 3dirunbounded
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param eps the smoothing length (only used for HEJ kernels)
 * 
 */
void cmpt_Green_3D_3dirunbounded_0dirspectral(const Topology *topo, const double hfact[3], const double symstart[3], double *green, GreenType typeGreen, const double eps){
    BEGIN_FUNC;

    FLUPS_CHECK(!(topo->isComplex()),"Green topology cannot been complex with 0 dir spectral", LOCATION);

    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(hfact[0] != 0.0, "grid spacing cannot be 0", LOCATION);
    FLUPS_CHECK(hfact[1] != 0.0, "grid spacing cannot be 0", LOCATION);
    FLUPS_CHECK(hfact[2] != 0.0, "grid spacing cannot be 0", LOCATION);

    double      G0;  //value of G in 0
    GreenKernel_unb G;
    int GN = 0;
    double* Gdata = NULL;

    switch (typeGreen) {
        case HEJ_2:
            G  = &_hej_2_3unb0spe;
            G0 =       M_SQRT2 / (4.0 * eps * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_4:
            G  = &_hej_4_3unb0spe;
            G0 = 3.0 * M_SQRT2 / (8.0 * eps * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_6:
            G  = &_hej_6_3unb0spe;
            G0 = 15.0 * M_SQRT2 / (32.0 * eps * sqrt(M_PI * M_PI * M_PI));
            break;
        case CHAT_2:
            G  = &_chat_2_3unb0spe;
            G0 = .5 * pow(1.5 * c_1o2pi * hfact[0] * hfact[1] * hfact[2], 2. / 3.);
            break;
        case LGF_2:
            FLUPS_CHECK(hfact[0] == hfact[1],"the grid has to be isotropic to use the LGFs",LOCATION);
            FLUPS_CHECK(hfact[1] == hfact[2],"the grid has to be isotropic to use the LGFs",LOCATION);
            // read the LGF data and store it
            _lgf_readfile(&GN,&Gdata);
            // associate the Green's function
            G  = &_lgf_2_3unb0spe;
            G0 = Gdata[0];
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.", LOCATION);
    }

    int istart[3];
    topo->get_istart_glob(istart);

    const int nf = topo->nf();
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            //local indexes start
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf);

            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                int is[3];
                cmpt_symID(ax0, i0, i1, i2, istart, symstart, 0, is);

                // symmetrized position
                const double x0 = (is[ax0]) * hfact[ax0];
                const double x1 = (is[ax1]) * hfact[ax1];
                const double x2 = (is[ax2]) * hfact[ax2];

                // green function value
                const double r2 = x0 * x0 + x1 * x1 + x2 * x2;
                const double r  = sqrt(r2);

                // the first two arguments are used in standard kernels and the others 5 ones are aimed for LGFs only
                const double tmp[7] = {r, eps, is[ax0], is[ax1], is[ax2],GN,hfact[ax0]};
                green[id + i0 * nf] = -G(tmp,Gdata);
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0) {
        green[0] = -G0;
    }
    // free Gdata if needed
    if (Gdata != NULL) {
        flups_free(Gdata);
    }

    END_FUNC;
}
/**@} */


/**
 * @name 2 directions unbounded - 1 direction spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
static inline double _hej_2_2unb1spe_k0(const void* params) {
    const double r   = ((double*)params)[0];
    const double sig = ((double*)params)[2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    // return -c_1o2pi * (log(r) - exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2)); //mistaken coefs in [Spietz2018]
    return -c_1o2pi * (log(r) + .5 * expint_ei(rho2 / 2));
    // return -c_1o2pi * (.5*log(rho*.5) + .5 * expint_ei(rho2 / 2));
}
static inline double _hej_2_2unb1spe_r0(const void* params) {
    const double sig = ((double*)params)[2];

    return c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig));
}

static inline double _hej_4_2unb1spe_k0(const void* params) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho  = r/sig;
    const double rho2 = rho*rho;
    // return -c_1o2pi * (log(r) - (1 - .5 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2)); //mistaken coefs in [Spietz2018]
    return -c_1o2pi * (log(r) - exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2));
    // return -c_1o2pi * (.5*log(rho2*.5) - exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2));
}
static inline double _hej_4_2unb1spe_r0(const void* params) {
    const double sig = ((double*)params)[2];

    return c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig) + .5);
}

static inline double _hej_6_2unb1spe_k0(const void* params) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rho2 = rho*rho;
    // return -c_1o2pi * (log(r) - (1 - rho2 + .125 * rho2 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2)); //mistaken coefs in [Spietz2018]
    return -c_1o2pi * (log(r) - (.75 - .125 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2));
    // return -c_1o2pi * (.5*log(rho2*.5) - (.75 - .125 * rho2) * exp(-rho2 / 2) + .5 * expint_ei(rho2 / 2));
}
static inline double _hej_6_2unb1spe_r0(const void* params) {
    const double sig = ((double*)params)[2];

    return c_1o2pi * (c_gamma * .5 - log(M_SQRT2 * sig) + .75);
}
static inline double _zero(const void* params) {   
    return - 0.0;
}

static inline double _chat_2_2unb1spe(const void* params) {
    const double r      = ((double*)params) [0];
    const double k      = ((double*)params) [1];

    return c_1o2pi * besselk0(fabs(k) * r);
}
static inline double _chat_2_2unb1spe_r0(const void* params) {
    const double k      = ((double*)params) [1];
    const double r_eq2D = ((double*)params) [3];

    return (1.0 - k * r_eq2D * besselk1(k * r_eq2D)) * c_1opi / ((k * r_eq2D) * (k * r_eq2D));
}
static inline double _chat_2_2unb1spe_k0(const void* params) {
    const double r      = ((double*)params) [0];
    // const double sig = ((double*)params)[2];
    
    return  - c_1o2pi * log(r) ; //caution: mistake on the sign in [Chatelain2010]
}

/**
 * @brief Compute the Green function for 2dirunbounded and 1dirspectral
 * 
 * The wave number in each direction is obtained as k_i = (i_s + koffset_i) * kfact_i, where is the global (potentially symmetric) index
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors (must be 0 in the spedtral dir)
 * @param kfact the k multiplicative factor (must be 0 in the unbounded dir)
 * @param koffset the k additive factor
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param eps the smoothing length (only used for HEJ kernels)
 */
void cmpt_Green_3D_2dirunbounded_1dirspectral(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps) {
    BEGIN_FUNC;
    
    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[0] != hfact[0], "grid spacing[0] cannot be = to dk[0]", LOCATION);
    FLUPS_CHECK(kfact[1] != hfact[1], "grid spacing[1] cannot be = to dk[1]", LOCATION);
    FLUPS_CHECK(kfact[2] != hfact[2], "grid spacing[2] cannot be = to dk[2]", LOCATION);

    // @Todo For Helmolz, we need Green to be complex 
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.", LOCATION);
    // opt_double_ptr mygreen = green; //casting of the Green function to be able to access real and complex part
    //Implementation note: if you want to do Helmolz, you need Hankel functions (3rd order Bessel) which are not implemented in stdC. Consider the use of boost lib.
    //notice that bessel_k has been introduced in c++17
    
    GreenKernel G;    // the Green kernel (general expression in the whole domain)
    GreenKernel Gk0;  // the Green kernel (particular expression in k=0)
    GreenKernel Gr0;  // the Green kernel (particular expression in r=0)

    switch (typeGreen) {
        case HEJ_2:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation.", LOCATION);
            
            // Note: 
            // According to [Spietz2018], we can obtain the **approximate** Green kernel by using the 2D unbounded kernel 
            // for mode 0 in the spectral direction, and the rest of the Green kernel is the same as in full spectral.
            // We here fill with zero the greatest part of Green: we are actually interested only in doing the FFT
            // of _hej_*_2unb1spe_k0 in the 2 remaining spatial directions. We will complete the Green function with the
            // full spectral part afterwards. 
            G   = &_zero;
            Gk0 = &_hej_2_2unb1spe_k0;
            Gr0 = &_hej_2_2unb1spe_r0;
            break;
        case HEJ_4:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation.",LOCATION);
            G   = &_zero;
            Gk0 = &_hej_4_2unb1spe_k0;
            Gr0 = &_hej_4_2unb1spe_r0;
            break;
        case HEJ_6:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation.",LOCATION);
            G   = &_zero;
            Gk0 = &_hej_6_2unb1spe_k0;
            Gr0 = &_hej_6_2unb1spe_r0;
            break;
        case CHAT_2:
            G   = &_chat_2_2unb1spe;
            Gk0 = &_chat_2_2unb1spe_k0;
            Gr0 = &_chat_2_2unb1spe_r0;
            // caution: the value of G in k=r=0 is specified at the end of this routine
            break;
        case LGF_2:
            FLUPS_ERROR("Lattice Green Function not implemented yet.", LOCATION);
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.", LOCATION);
    }

    int istart[3];
    topo->get_istart_glob(istart);

    const int    nf      = topo->nf();
    const int    ax0     = topo->axis();
    const int    ax1     = (ax0 + 1) % 3;
    const int    ax2     = (ax0 + 2) % 3;
    const int    nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
    const double r_eq2D  = c_1osqrtpi * sqrt(hfact[ax0] * hfact[ax1] + hfact[ax1] * hfact[ax2] + hfact[ax2] * hfact[ax0]);

    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            //local indexes start
            const size_t id = localIndex(ax0,0, i1, i2, ax0, nmem,nf);
        
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                
                // global indexes
                int is[3];
                cmpt_symID(ax0,i0,i1,i2,istart,symstart,0,is);

                // (symmetrized) wave number : only one kfact is non-zero
                const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];
                const double k = k0 + k1 + k2;

                //(symmetrized) position : only one hfact is zero
                const double x0 = (is[ax0]) * hfact[ax0];
                const double x1 = (is[ax1]) * hfact[ax1];
                const double x2 = (is[ax2]) * hfact[ax2];
                const double r  = sqrt(x0 * x0 + x1 * x1 + x2 * x2);

                const double tmp[4] = {r, k, eps, r_eq2D};

                // green function value
                // Implementation note: having a 'if' in a loop is highly discouraged... however, this is the init so we prefer having a
                // this routine with a high readability and lower efficency than the opposite.
                if (r <= (hfact[ax0] + hfact[ax1] + hfact[ax2]) * .2) {
                    green[id + i0 * topo->nf()] = -Gr0(tmp);
                } else if (k <= (kfact[ax0] + kfact[ax1] + kfact[ax2]) * 0.2) {
                    green[id + i0 * topo->nf()] = -Gk0(tmp);
                } else {
                    green[id + i0 * topo->nf()] = -G(tmp);
                }
            }
        }
    }
    
    // reset the value in x=y=0.0 and k=0
    if (typeGreen == CHAT_2 && istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0) {
        // green[0] = -2.0 * log(1 + sqrt(2)) * c_1opiE3o2 / r_eq2D;
        green[0] = .25 * c_1o2pi * (M_PI - 6.0 + 2. * log(.5 * M_PI * r_eq2D));  //caution: mistake in [Chatelain2010]
    }
    END_FUNC;
}
/**@} */


/**
 * @name 1 direction unbounded - 2 directions spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
static inline double _hej_2_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;

    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return .25 * sig / s * subfun ;
}
static inline double _hej_2_1unb2spe_k0(const void* params) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rosqrt2 = r*c_1osqrt2;
    // return -.5* (r * erf(rosqrt2/sig) + (exp(-r*r/(2*sig*sig)) - 1.)*sig*M_SQRT2*c_1osqrtpi) ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
    return -.5* r * erf(rosqrt2/sig) + (1.-exp(-rho*rho*.5)) *sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double _hej_4_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return .25 * sig / s * subfun + \
           sig * M_SQRT2 * c_1osqrtpi * .25 * exp(-.5 * (s * s + rho * rho));
}
static inline double _hej_4_1unb2spe_k0(const void* params) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rosqrt2 = r*c_1osqrt2;
    return -.5* r * erf(rosqrt2/sig) + (1.-exp(-rho*rho*.5)) *.5*sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double _hej_6_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    const double subfun = s * rho > 100. ? 0 : ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
    return .25 * sig / s * subfun + \
           sig * M_SQRT2 * c_1osqrtpi * (c_5o16 + c_1o16 * (s * s - rho * rho)) * exp(-.5 * (s * s + rho * rho));
}
static inline double _hej_6_1unb2spe_k0(const void* params) {
    const double r   = ((double*)params) [0];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double rosqrt2 = r*c_1osqrt2;
    return -.5* r * erf(rosqrt2/sig) + (3.-exp(-rho*rho*.5) * (rho*rho+3.) ) *.125*sig*c_1osqrt2*c_1osqrtpi ; //mistakenly 0.0 in [Hejlesen:2013] and [Spietz:2018]
}

static inline double _chat_2_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];

    return .5 * exp(-k * r) / k;
}
static inline double _chat_2_1unb2spe_k0(const void* params) {
    const double r   = ((double*)params) [0];

    return -.5 * fabs(r);
}


/**
 * @brief Compute the Green function for 1dirunbounded and 2dirspectral
 * 
 * The wave number in each direction is obtained as k_i = (i_s + koffset_i) * kfact_i, where is the global (potentially symmetric) index
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors (must be 0 in the spedtral dir)
 * @param kfact the k multiplicative factor (must be 0 in the unbounded dir)
 * @param koffset the k additive factor
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param eps the smoothing length (only used for HEJ kernels)
 */
void cmpt_Green_3D_1dirunbounded_2dirspectral(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps) {
    BEGIN_FUNC;

    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[0] != hfact[0], "grid spacing[0] cannot be = to dk[0]", LOCATION);
    FLUPS_CHECK(kfact[1] != hfact[1], "grid spacing[1] cannot be = to dk[1]", LOCATION);
    FLUPS_CHECK(kfact[2] != hfact[2], "grid spacing[2] cannot be = to dk[2]", LOCATION);

    // @Todo For Helmolz, we need Green to be complex 
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.", LOCATION);
    // double* mygreen = green; //casting of the Green function to be able to access real and complex part

    GreenKernel G;   // the Green kernel (general expression in the whole domain)
    GreenKernel G0;  // the Green kernel (particular expression in k=0)

    switch (typeGreen) {
        case HEJ_2:
            G  = &_hej_2_1unb2spe;
            G0 = &_hej_2_1unb2spe_k0;
            break;
        case HEJ_4:
            G  = &_hej_4_1unb2spe;
            G0 = &_hej_4_1unb2spe_k0;
            break;
        case HEJ_6:
            G  = &_hej_6_1unb2spe;
            G0 = &_hej_6_1unb2spe_k0;
            break;
        case CHAT_2:
            G  = &_chat_2_1unb2spe;
            G0 = &_chat_2_1unb2spe_k0;
            break;
        case LGF_2:
            FLUPS_ERROR("Lattice Green Function not implemented yet.", LOCATION);
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.", LOCATION);
    }

    int istart[3];
    topo->get_istart_glob(istart);

    const int nf      = topo->nf();
    const int ax0     = topo->axis();
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
    // FLUPS_INFO("KFAC= %lf %lf %lf", kfact[0],kfact[1],kfact[2]);
    // FLUPS_INFO("Koff= %lf %lf %lf", koffset[0],koffset[1],koffset[2]);
    // FLUPS_INFO("HFAC= %lf %lf %lf", hfact[0],hfact[1],hfact[2]);
    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            //local indexes start
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf);

            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                int is[3];
                cmpt_symID(ax0,i0,i1,i2,istart,symstart,0,is);

                // (symmetrized) wave number : only 1 kfact is zero
                const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];
                const double k  = sqrt(k0 * k0 + k1 * k1 + k2 * k2);

                //(symmetrized) position : only 1 hfact is non-zero
                const double x = is[ax0] * hfact[ax0] + is[ax1] * hfact[ax1] + is[ax2] * hfact[ax2];

                const double tmp[3] = {x, k, eps};

                // green function value
                // Implementation note: having a 'if' in a loop is highly discouraged... however, this is the init so we prefer having a
                // this routine with a high readability and lower efficency than the opposite.
                if (k <= (kfact[ax0] + kfact[ax1] + kfact[ax2]) * 0.2) {
                    green[id + i0 * nf] = -G0(tmp);
                }
                else{
                    green[id + i0 * nf] = -G(tmp);
                }
            }
        }
    }
    END_FUNC;
}

/**@} */


/**
 * @name 3 directions spectral
 * 
 * @{
 */
// ----------------------------------------------------------- KERNELS ----------------------------------------------------------
static inline double _hej_2_0unb3spe(const void* params) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    return exp(-ssqr / 2) / (ksqr); 
}
static inline double _hej_4_0unb3spe(const void* params) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    return (1 + ssqr / 2) * exp(-ssqr / 2) / (ksqr);
}
static inline double _hej_6_0unb3spe(const void* params) {
    const double ksqr = ((double*)params)[0];
    const double sig  = ((double*)params)[1];

    const double ssqr = ksqr * (sig * sig);
    return (1 + ssqr / 2 + ssqr * ssqr / 8) * exp(-ssqr / 2) / (ksqr);
}

static inline double _chat_2_0unb3spe(const void* params) {
    const double ksqr   = ((double*)params) [0];

    return 1 / ksqr;
}

/**
 * @brief Compute the Green function for 3dirspectral (in the whole spectral domain)
 * 
 * __Note on performance__: obviously, the Green function in full spectral
 * is \f$\-frac{1}{k^2}\f$ (at least for CHAT_2). We could perform that operation directly in the
 * loop of `dothemagic`. We here choose to still precompute and store it.
 * We burn more memory, but we should fasten `dothemagic` as we replace a
 * (expensive) evaluation of \f$\frac{1}{k^2}\f$ by a memory access.
 *
 * The wave number in each direction is obtained as k_i = (i_s + koffset_i) * kfact_i, where is the global (potentially symmetric) index.
 * 
 * @param topo the topology associated to the Green's function
 * @param kfact the k multiplicative factor
 * @param koffset the k additive factor
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param eps the smoothing length (only used for HEJ kernels)
 */
void cmpt_Green_3D_0dirunbounded_3dirspectral(const Topology *topo, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps){
    cmpt_Green_3D_0dirunbounded_3dirspectral(topo, kfact, koffset, symstart, green, typeGreen, eps, NULL, NULL);
}


/**
 * @brief Compute the Green function for 3dirspectral (in a portion of the spectral domain)
 * 
 * The wave number in each direction is obtained as k_i = (i_s + koffset_i) * kfact_i, where is the global (potentially symmetric) index.
 * 
 * It is possible to fill only partially the spectral space, by specifying a global istart_custom and iend_custom which dictate the
 * span of global index in each direction that will be filled.
 * 
 * @param topo the topology associated to the Green's function
 * @param kfact the k multiplicative factor
 * @param koffset the k additive factor
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param eps the smoothing length (only used for HEJ kernels)
 * @param istart_custom global index where we start to fill data, in each dir. If NULL, we start at the beginning of the spectral space.
 * @param iend_custom global index where we end to fill data, in each dir. If NULL, we end at the end of the spectral space.
 */
void cmpt_Green_3D_0dirunbounded_3dirspectral(const Topology *topo, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps, const int istart_custom[3], const int iend_custom[3]){
    BEGIN_FUNC;

    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(kfact[0] != 0.0, "dk cannot be 0", LOCATION);
    FLUPS_CHECK(kfact[1] != 0.0, "dk cannot be 0", LOCATION);
    FLUPS_CHECK(kfact[2] != 0.0, "dk cannot be 0", LOCATION);

    GreenKernel G;   // the Green kernel (general expression in the whole domain)

    switch (typeGreen) {
        case HEJ_2:
            G = &_hej_2_0unb3spe;
            break;
        case HEJ_4:
            G = &_hej_4_0unb3spe;
            break;
        case HEJ_6:
            G = &_hej_6_0unb3spe;
            break;
        case CHAT_2:
            G = &_chat_2_0unb3spe;
            break;
        case LGF_2:
            FLUPS_ERROR("Lattice Green Function not implemented yet.", LOCATION);
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.", LOCATION);
    }

    int istart[3];
    topo->get_istart_glob(istart);

    int is_[3] = {0, 0, 0};
    int ie_[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)};

    if (istart_custom != NULL) {
        //switch to local index
        for (int ip = 0; ip < 3; ip++) {
            is_[ip] = fmin( fmax(istart_custom[ip]- istart[ip], 0), topo->nloc(ip));
        }
    } 
    if (iend_custom != NULL) {
        //switch to local index
        for (int ip = 0; ip < 3; ip++) {
            ie_[ip] = fmin( fmax(iend_custom[ip]- istart[ip], 0), topo->nloc(ip));
        }
    }
    const int is[3] = {is_[0], is_[1], is_[2]};
    const int ie[3] = {ie_[0], ie_[1], ie_[2]};

    // printf("IS_0 : %d,%d,%d \n",is[0],is[1],is[2]);
    // printf("IS_E : %d,%d,%d \n",ie[0],ie[1],ie[2]);
    // printf("ISTART : %d,%d,%d \n",istart[0],istart[1],istart[2]);
    // printf("K_OFFSET : %lf,%lf,%lf \n",koffset[0],koffset[1],koffset[2]);
    // FLUPS_INFO("KFAC= %lf %lf %lf", kfact[0],kfact[1],kfact[2]);
    // FLUPS_INFO("HFAC= %lf %lf %lf", hfact[0],hfact[1],hfact[2]);
    // printf("IEND : %d,%d,%d \n",topo->nloc(0),topo->nloc(1),topo->nloc(2));

    const int nf      = topo->nf();
    const int ax0     = topo->axis();
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
    for (int i2 = is[ax2]; i2 < ie[ax2]; i2++) {
        for (int i1 = is[ax1]; i1 < ie[ax1]; i1++) {
            //local indexes start
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem,nf);
            for (int i0 = is[ax0]; i0 < ie[ax0]; i0++) {
                int il[3];
                cmpt_symID(ax0, i0, i1, i2, istart, symstart, 0, il);
                //the previous call works with koffset below because there is never a shiftgreen AND a symstart together

                // (symmetrized) wave number
                const double k0 = (il[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (il[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (il[ax2] + koffset[ax2]) * kfact[ax2];

                // green function value
                const double ksqr = k0 * k0 + k1 * k1 + k2 * k2;

                const double tmp[2] = {ksqr, eps};

                green[id + i0 * nf] = -G(tmp);
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0 \
        && koffset[0]+koffset[1]+koffset[2]<0.2 ) {
        green[0] = -0.0;
    }
    END_FUNC;
}
/**@} */
