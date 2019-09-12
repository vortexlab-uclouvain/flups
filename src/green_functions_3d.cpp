/**
 * @file Green_functions_3d.cpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * -------------------------------
 * **Symmetry computation:**
 * 
 * We have to take the symmetry around symstart. e.g. in X direction: `symstart[0] - (ix - symstart[0]) = 2 symstart[0] - ix`
 * 
 * In some cases when we have an R2C transform, it ask for 2 additional doubles.
 * The value is meaningless but we would like to avoid segfault and nan's.
 * To do so, we use 2 tricks:
 * - The `abs` is used to stay on the positivie side and hence avoid negative memory access
 * - The `max` is used to prevent the computation of the value in 0, which is never used in the symmetry.
 * 
 * As an example, the final formula is then ( in the X direction):
 * `max( abs(2 symstart[0] - ix) , 1)`
 * 
 */

#include "green_functions_3d.hpp"


using namespace FLUPS;

/**
 * @brief generic type for Green kernel, takes a table of parameters that can be used depending on the kernel
 * 
 */
typedef double (*GreenKernel)(const void* );

//notice that these function will likely not be inlined as we have a pointer to them...
static inline double _hej_2_3unb0spe(const void* params) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    return c_1o4pi / r * (erf(r / eps * c_1osqrt2));
}
static inline double _hej_4_3unb0spe(const void* params) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (rho)*exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double _hej_6_3unb0spe(const void* params) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (c_7o4 * rho - c_1o4 * pow(rho, 3)) * exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double _chat_2_3unb0spe(const void* params) {
    double r   = ((double*)params) [0];
    return c_1o4pi / r ;
}

/**
 * @brief Compute the Green function for 3dirunbounded
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors
 * @param symstart the symmetry plan asked
 * @param green the Green function array
 * @param typeGreen the type of Green function to be 
 * @param eps the smoothing length (only used for HEJ kernels)
 * 
 */
void cmpt_Green_3D_3dirunbounded_0dirspectral(const Topology *topo, const double hfact[3], const double symstart[3], double *green, GreenType typeGreen, const double eps){
    BEGIN_FUNC;

    FLUPS_CHECK(!(topo->isComplex()),"Green topology cannot been complex with 0 dir spectral");

    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(hfact[0] != 0.0, "grid spacing cannot be 0");
    FLUPS_CHECK(hfact[1] != 0.0, "grid spacing cannot be 0");
    FLUPS_CHECK(hfact[2] != 0.0, "grid spacing cannot be 0");

    int ax0 = topo->axis();
    int ax1 = (ax0 + 1) % 3;
    int ax2 = (ax0 + 2) % 3;
    
    double      G0;  //value of G in 0
    GreenKernel G;

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
            FLUPS_ERROR("Lattice Green Function not implemented yet.");
            //please add the parameters you need to params
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.");
    }

    int istart[3];
    get_istart_glob(istart,topo);

    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                int is[3];
                cmpt_symID(ax0,i0,i1,i2,istart,symstart,0,is);

                // symmetrized position
                const double x0 = (is[ax0])*hfact[ax0];
                const double x1 = (is[ax1])*hfact[ax1];
                const double x2 = (is[ax2])*hfact[ax2];

                // green function value
                const double r2 = x0 * x0 + x1 * x1 + x2 * x2;
                const double r  = sqrt(r2);
                const size_t id = i0 + topo->nloc(ax0) * (i1 + i2 * topo->nloc(ax1));

                const double tmp[2] = {r, eps};
                green[id]           = -G(tmp);
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0) {
        green[0] = -G0;
    }
}


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
    const double sig = ((double*)params)[2];
    
    return  - c_1o2pi * log(r) ; //caution: mistake on the sign in [Chatelain2010]
}

/**
 * @brief 
 * 
 * @param topo must be the topo in which ax0 is the spectral dir
 * @param hfact 
 * @param kfact 
 * @param koffset 
 * @param symstart 
 * @param green 
 * @param typeGreen 
 * @param eps 
 */
void cmpt_Green_3D_2dirunbounded_1dirspectral(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps) {
    const int ax0 = topo->axis();  
    const int ax1 = (ax0 + 1) % 3; 
    const int ax2 = (ax0 + 2) % 3; 

    // printf("kfact - hfact : %lf,%lf,%lf - %lf,%lf,%lf\n",kfact[ax0],kfact[ax1],kfact[ax2],hfact[ax0],hfact[ax1],hfact[ax2]);

    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[ax0] != hfact[ax0], "grid spacing[0] cannot be = to dk[0]");
    FLUPS_CHECK(kfact[ax1] != hfact[ax1], "grid spacing[1] cannot be = to dk[1]");
    FLUPS_CHECK(kfact[ax2] != hfact[ax2], "grid spacing[2] cannot be = to dk[2]");

    // @Todo For Helmolz, we need Green to be complex 
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.");
    // opt_double_ptr mygreen = green; //casting of the Green function to be able to access real and complex part
    //Implementation note: if you want to do Helmolz, you need Hankel functions (3rd order Bessel) which are not implemented in stdC. Consider the use of boost lib.
    //notice that bessel_k has been introduced in c++17
    
    GreenKernel G;    // the Green kernel (general expression in the whole domain)
    GreenKernel Gk0;  // the Green kernel (particular expression in k=0)
    GreenKernel Gr0;  // the Green kernel (particular expression in r=0)

    switch (typeGreen) {
        case HEJ_2:
            // FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail a approximation.");
            
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
            // FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail a approximation.");
            G   = &_zero;
            Gk0 = &_hej_4_2unb1spe_k0;
            Gr0 = &_hej_4_2unb1spe_r0;
            break;
        case HEJ_6:
            // FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail a approximation.");
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
            FLUPS_ERROR("Lattice Green Function not implemented yet.");
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.");
    }

    int istart[3];
    get_istart_glob(istart,topo);

    const double r_eq2D = c_1osqrtpi * sqrt( hfact[ax0]*hfact[ax1]+hfact[ax1]*hfact[ax2]+hfact[ax2]*hfact[ax0] );

    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            
            //local indexes start
            size_t id = localindex_ao(0, i1, i2, topo);
        
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
}



static inline double _hej_2_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    return .25 * sig / s * ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho));
}
static inline double _hej_4_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    return .25 * sig / s * ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho)) + \
           sig * M_SQRT2 * c_1osqrtpi * .25 * exp(-.5 * (s * s + rho * rho));
}
static inline double _hej_6_1unb2spe(const void* params) {
    const double r   = ((double*)params) [0];
    const double k   = ((double*)params) [1];
    const double sig = ((double*)params) [2];

    const double rho = r/sig;
    const double s   = k*sig;
    return .25 * sig / s * ((1 - erf(c_1osqrt2 * (s - rho))) * exp(-s * rho) + (1 - erf(c_1osqrt2 * (s + rho))) * exp(s * rho)) + \
           sig * M_SQRT2 * c_1osqrtpi * (c_5o16 + c_1o16 * (s * s - rho * rho)) * exp(-.5 * (s * s + rho * rho));
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
 * @brief 
 * 
 * @param topo 
 * @param hfact 
 * @param kfact 
 * @param koffset 
 * @param symstart 
 * @param green 
 * @param typeGreen 
 * @param eps 
 */
void cmpt_Green_3D_1dirunbounded_2dirspectral(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps) {
    const int ax0 = topo->axis();  
    const int ax1 = (ax0 + 1) % 3; 
    const int ax2 = (ax0 + 2) % 3; 

    // printf("kfact - hfact : %lf,%lf,%lf - %lf,%lf,%lf\n",kfact[ax0],kfact[ax1],kfact[ax2],hfact[ax0],hfact[ax1],hfact[ax2]);

    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[ax0] != hfact[ax0], "grid spacing[0] cannot be = to dk[0]");
    FLUPS_CHECK(kfact[ax1] != hfact[ax1], "grid spacing[1] cannot be = to dk[1]");
    FLUPS_CHECK(kfact[ax2] != hfact[ax2], "grid spacing[2] cannot be = to dk[2]");

    // @Todo For Helmolz, we need Green to be complex 
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.");
    // opt_double_ptr mygreen = green; //casting of the Green function to be able to access real and complex part

    GreenKernel G;   // the Green kernel (general expression in the whole domain)
    GreenKernel G0;  // the Green kernel (particular expression in k=0)

    switch (typeGreen) {
        case HEJ_2:
            G  = &_hej_2_1unb2spe;
            G0 = &_hej_2_1unb2spe;
            break;
        case HEJ_4:
            G  = &_hej_4_1unb2spe;
            G0 = &_hej_4_1unb2spe;
            break;
        case HEJ_6:
            G  = &_hej_6_1unb2spe;
            G0 = &_hej_6_1unb2spe;
            break;
        case CHAT_2:
            G  = &_chat_2_1unb2spe;
            G0 = &_chat_2_1unb2spe_k0;
            break;
        case LGF_2:
            FLUPS_ERROR("Lattice Green Function not implemented yet.");
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.");
    }

    int istart[3];
    get_istart_glob(istart,topo);


    //Note: i0 (ax0) is the only spatial (i.e. non spectral) axis
    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            //local indexes start
            size_t id = localindex_ao(0, i1, i2, topo);

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
                if (k <= (kfact[ax0]+kfact[ax1]+kfact[ax2])*0.2 ){
                    green[id + i0 * topo->nf()] = -G0(tmp);
                }
                else{
                    green[id + i0 * topo->nf()] = -G(tmp);
                }
            }
        }
    }
}

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
 * @brief Compute the Green function for 3dirspectral
 * 
 * __Note on performance__: obviously, the Green function in full spectral
 * is \f$\-frac{1}{k^2}\f$ (at least for CHAT_2). We could perform that operation directly in the
 * loop of `dothemagic`. We here choose to still precompute and store it.
 * We burn more memory, but we should fasten `dothemagic` as we replace a
 * (expensive) evaluation of \f$\frac{1}{k^2}\f$ by a memory access.
 * 
 * @param topo 
 * @param kfact 
 * @param koffset
 * @param symstart 
 * @param green 
 * @param typeGreen 
 * @param alpha 
 */
void cmpt_Green_3D_0dirunbounded_3dirspectral(const Topology *topo, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps){
    const int istart0[3] = {0, 0, 0};
    const int ishift[3]  = {0, 0, 0};
    cmpt_Green_3D_0dirunbounded_3dirspectral(topo, kfact, koffset, symstart, green, typeGreen, eps, istart0, ishift);
}


void cmpt_Green_3D_0dirunbounded_3dirspectral(const Topology *topo, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double eps, const int istart0[3], const int ishift[3]){
    BEGIN_FUNC;

    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(kfact[0] != 0.0, "dk cannot be 0");
    FLUPS_CHECK(kfact[1] != 0.0, "dk cannot be 0");
    FLUPS_CHECK(kfact[2] != 0.0, "dk cannot be 0");

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    printf("kfact : %lf,%lf,%lf \n",kfact[ax0],kfact[ax1],kfact[ax2]);
    // printf("koff  : %lf,%lf,%lf \n",koffset[ax0],koffset[ax1],koffset[ax2]);   

    printf("eps : %lf \n",eps);
 
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
            FLUPS_ERROR("Lattice Green Function not implemented yet.");
            break;
        default:
            FLUPS_ERROR("Green Function type unknow.");
    }

    int istart[3];
    get_istart_glob(istart,topo);

    //forgetting a part of the domain
    const int is0 = istart[ax0] == 0 ? istart0[ax0] : 0;
    const int is1 = istart[ax1] == 0 ? istart0[ax1] : 0;
    const int is2 = istart[ax2] == 0 ? istart0[ax2] : 0;
    printf("IS0 : %d,%d,%d \n",is0,is1,is2);

    for (int ip; ip<3;ip++){
        istart[ip] += ishift[ip];
    }
    printf("ISTART : %d,%d,%d \n",istart[0],istart[1],istart[2]);


    for (int i2 = is2; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = is1; i1 < topo->nloc(ax1); i1++) {
            //local indexes start
            size_t id = localindex_ao(0, i1, i2, topo);

            for (int i0 = is0; i0 < topo->nloc(ax0); i0++) {
                int is[3];
                cmpt_symID(ax0,i0,i1,i2,istart,symstart,0,is);

                // (symmetrized) wave number
                const double k0 = (is[ax0]+koffset[ax0])*kfact[ax0];
                const double k1 = (is[ax1]+koffset[ax1])*kfact[ax1];
                const double k2 = (is[ax2]+koffset[ax2])*kfact[ax2];

                // green function value
                const double ksqr = k0 * k0 + k1 * k1 + k2 * k2;

                const double tmp[2] = {ksqr, eps};
                
                green[id + i0*topo->nf()] = -G(tmp);
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0 && koffset[0]+koffset[1]+koffset[2]<0.2 && ishift[0]==0 && ishift[1]==0 && ishift[2]==0) {
        green[0] = -0.0;
    }
}