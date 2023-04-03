/**
 * @file green_functions.cpp
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

#include "green_functions.hpp"
#include "green_kernels.hpp"
#include "ji0.hpp"

/**
 * @brief generic type for Green kernel, takes a table of parameters that can be used depending on the kernel
 * 
 */
typedef double (*GreenKernel)(const void*,const double*);

/**
 * @brief Compute the Green function for 0 dir spectral (i.e. 3 dir unbounded or 2 dirunbounded)
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param length the characteristic length (only used for HEJ kernels = epsilon)
 */
void cmpt_Green_3dirunbounded(const Topology *topo, const double hfact[3], const double symstart[3], double *green, GreenType typeGreen, const double length){
    BEGIN_FUNC;
    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(hfact[0] != 0.0, "hfact[0] cannot be 0");
    FLUPS_CHECK(hfact[1] != 0.0, "hfact[1] cannot be 0");
    FLUPS_CHECK(hfact[2] != 0.0, "hfact[2] cannot be 0");

    // FLUPS_INFO("K_OFFSET : %lf,%lf,%lf \n",koffset[0],koffset[1],koffset[2]);
    // FLUPS_INFO("KFAC= %lf %lf %lf", kfact[0],kfact[1],kfact[2]);
    // FLUPS_INFO("HFAC= %lf %lf %lf", hfact[0],hfact[1],hfact[2]);
    
    GreenKernel G;

    double  G0;  //value of G in 0
    int     GN    = 0;
    double *Gdata = NULL;

    //==========================    3D  =================================
    switch (typeGreen) {
        case HEJ_2:
            G  = &hej_2_3unb0spe_;
            G0 = - M_SQRT2 / (4.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_4:
            G  = &hej_4_3unb0spe_;
            G0 = - 3.0 * M_SQRT2 / (8.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_6:
            G  = &hej_6_3unb0spe_;
            G0 = - 15.0 * M_SQRT2 / (32.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_8:
            G  = &hej_8_3unb0spe_;
            G0 = - 35.0 * M_SQRT2 / (64.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_10:
            G  = &hej_10_3unb0spe_;
            G0 = - 315.0 * M_SQRT2 / (512.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_0:
            G  = &hej_0_3unb0spe_;
            G0 = - 1.0/(2.0*M_PI*M_PI*length);
            break;
        case CHAT_2:
            G  = &chat_2_3unb0spe_;
            G0 = - 0.5 * pow(1.5 * c_1o2pi * hfact[0] * hfact[1] * hfact[2], 2. / 3.);
            break;
        case LGF_2:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_2, 3, &GN, &Gdata);
            // associate the Green's function
            G = &lgf_2_3unb0spe_;
            break;
        case LGF_4:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_4, 3, &GN, &Gdata);
            G = &lgf_4_3unb0spe_;
            break;
        case LGF_6:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_6, 3, &GN, &Gdata);
            G = &lgf_6_3unb0spe_;
            break;
        case LGF_8:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_8, 3, &GN, &Gdata);
            G = &lgf_8_3unb0spe_;
            break;
        case MEHR_4L:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_4L, 3, &GN, &Gdata);
            G = &mehr_4l_3unb0spe_;
            break;
        case MEHR_6L:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_6L, 3, &GN, &Gdata);
            G = &mehr_6l_3unb0spe_;
            break;
        case MEHR_4F:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_4F, 3, &GN, &Gdata);
            G = &mehr_4f_3unb0spe_;
            break;
        case MEHR_6F:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_6F, 3, &GN, &Gdata);
            G = &mehr_6f_3unb0spe_;
            break;
        default:
            FLUPS_CHECK(false, "Green Function type unknown.");
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
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

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

                // the first two arguments are used in standard kernels, the two zeros are for compatibility with the 2dirunbounded function,
                // and the others 5 ones are aimed for LGFs only
                // the symmetrized indexes will be negative!!
                const double tmp[9] = {r, length, 0, 0, static_cast<double>(std::abs(is[ax0])), static_cast<double>(std::abs(is[ax1])), static_cast<double>(std::abs(is[ax2])), static_cast<double>(GN), hfact[ax0]};
                green[id + i0 * nf] = G(tmp,Gdata);
            }
        }
    }
    // reset the value in 0.0 but not for LGF's since we have already pre-computed its value
    const bool is_lgf = (typeGreen == LGF_2 || typeGreen == LGF_4 || typeGreen == LGF_6 || typeGreen == LGF_8);
    const bool is_mehr = (typeGreen == MEHR_4F || typeGreen == MEHR_6F || typeGreen == MEHR_4L || typeGreen == MEHR_6L);
    const bool correct_istart = (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0);
    if (!is_lgf && !is_mehr && correct_istart) {
        green[0] = G0;
    }
    // free Gdata if needed
    if (Gdata != NULL) {
        m_free(Gdata);
    }

    END_FUNC;
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
 * @param typeGreen  the type of Green function 
 * @param length the characteristic length (only used for HEJ kernels = epsilon)
 * 
 * @warning For 3D kernels: According to [Spietz2018], we can obtain the **approximate** Green kernel by using the 2D unbounded kernel 
            for mode 0 in the spectral direction, and the rest of the Green kernel is the same as in full spectral.
            We here fill with zero most part of Green data. Indeed, we are interested only in doing the FFT
            of hej__*2unb1spe_k0_ in the 2 remaining spatial directions. We will complete the Green function with the
            full spectral part afterwards, while going through Solver::cmptGreenFunction_.
 * 
 */
void cmpt_Green_2dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length) {
    BEGIN_FUNC;
    
    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[0] != hfact[0], "grid spacing[0] cannot be = to dk[0]");
    FLUPS_CHECK(kfact[1] != hfact[1], "grid spacing[1] cannot be = to dk[1]");
    // check that if hfact or kfact != 0, they are not the same
    FLUPS_CHECK(!(kfact[2] == hfact[2] && (kfact[2]!= 0.0 || hfact[2] != 0.0)), "grid spacing[2] cannot be = to dk[2]");

    // @Todo For Helmolz, we need Green to be complex 
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.");
    // opt_double_ptr mygreen = green; //casting of the Green function to be able to access real and complex part
    //Implementation note: if you want to do Helmolz, you need Hankel functions (3rd order Bessel) which are not implemented in stdC. Consider the use of boost lib.
    //notice that bessel_k has been introduced in c++17

    GreenKernel G;    // the Green kernel (general expression in the whole domain)
    GreenKernel Gk0;  // the Green kernel (particular expression in k=0)
    GreenKernel Gr0;  // the Green kernel (particular expression in r=0)

    int     GN    = 0;
    double *Gdata = NULL;

    switch (typeGreen) {
        case HEJ_2:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation in 3D.");
            // see warning in the function description
            G   = &zero_;
            Gk0 = &hej_2_2unb1spe_k0_;
            Gr0 = &hej_2_2unb1spe_r0_;
            break;
        case HEJ_4:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation in 3D.");
            G   = &zero_;
            Gk0 = &hej_4_2unb1spe_k0_;
            Gr0 = &hej_4_2unb1spe_r0_;
            break;
        case HEJ_6:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation in 3D.");
            G   = &zero_;
            Gk0 = &hej_6_2unb1spe_k0_;
            Gr0 = &hej_6_2unb1spe_r0_;
            break;
        case HEJ_8:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation in 3D.");
            G   = &zero_;
            Gk0 = &hej_8_2unb1spe_k0_;
            Gr0 = &hej_8_2unb1spe_r0_;
            break;
        case HEJ_10:
            FLUPS_WARNING("HEJ kernels in 2dirunbounded 1dirspectral entail an approximation in 3D.");
            G   = &zero_;
            Gk0 = &hej_10_2unb1spe_k0_;
            Gr0 = &hej_10_2unb1spe_r0_;
            break;        
        case HEJ_0:
            FLUPS_WARNING("HEJ0 (theoretically spectral) kernel for 2D unbounded entails an approximation greatly affecting accuracy.");
            init_Ji0();
            G   = &zero_;
            Gk0 = &hej_0_2unb1spe_k0_;
            Gr0 = &hej_0_2unb1spe_k0_;
            break;        
        case CHAT_2:
            G   = &chat_2_2unb1spe_;
            Gk0 = &chat_2_2unb1spe_k0_;
            Gr0 = &chat_2_2unb1spe_r0_;
            // caution: the value of G in k=r=0 is specified at the end of this routine
            break;
        case LGF_2:
            FLUPS_CHECK(hfact[2] < 1.0e-14, "This LGF cannot be called in a 3D problem -> h[3] = %e",hfact[3]);
            FLUPS_CHECK(hfact[0] == hfact[1], "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_2, 2, &GN, &Gdata);
            // associate the Green's function
            G   = &zero_;
            Gk0 = &lgf_2_2unb0spe_;
            Gr0 = &lgf_2_2unb0spe_;
            break;
        case LGF_4:
            FLUPS_CHECK(hfact[2] < 1.0e-14, "This LGF cannot be called in a 3D problem -> h[3] = %e",hfact[3]);
            FLUPS_CHECK(hfact[0] == hfact[1], "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_4, 2, &GN, &Gdata);
            G   = &zero_;
            Gk0 = &lgf_4_2unb0spe_;
            Gr0 = &lgf_4_2unb0spe_;
            break;
        case LGF_6:
            FLUPS_CHECK(hfact[2] < 1.0e-14, "This LGF cannot be called in a 3D problem -> h[3] = %e",hfact[3]);
            FLUPS_CHECK(hfact[0] == hfact[1], "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_6, 2, &GN, &Gdata);
            G   = &zero_;
            Gk0 = &lgf_6_2unb0spe_;
            Gr0 = &lgf_6_2unb0spe_;
            break;
        case LGF_8:
            FLUPS_CHECK(hfact[2] < 1.0e-14, "This LGF cannot be called in a 3D problem -> h[3] = %e",hfact[3]);
            FLUPS_CHECK(hfact[0] == hfact[1], "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_8, 2, &GN, &Gdata);
            G   = &zero_;
            Gk0 = &lgf_8_2unb0spe_;
            Gr0 = &lgf_8_2unb0spe_;
            break;
        case MEHR_4L:
            FLUPS_CHECK(false, "MEHR4 kernel not available for 2D unbounded problems.");
            break;
        case MEHR_6L:
            FLUPS_CHECK(false, "MEHR6 kernel not available for 2D unbounded problems.");
            break;
        case MEHR_4F:
            FLUPS_CHECK(false, "MEHR4 kernel not available for 2D unbounded problems.");
            break;
        case MEHR_6F:
            FLUPS_CHECK(false, "MEHR6 kernel not available for 2D unbounded problems.");
            break;
        default:
            FLUPS_CHECK(false, "Green Function type unknown.");
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
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                // global indexes
                int is[3];
                cmpt_symID(ax0, i0, i1, i2, istart, symstart, 0, is);

                // (symmetrized) wave number : only one kfact is non-zero
                const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];
                const double k  = k0 + k1 + k2;

                //(symmetrized) position : only one hfact is zero
                const double x0 = (is[ax0]) * hfact[ax0];
                const double x1 = (is[ax1]) * hfact[ax1];
                const double x2 = (is[ax2]) * hfact[ax2];
                const double r  = sqrt(x0 * x0 + x1 * x1 + x2 * x2);

                // the symmetrized indexes will be negative!!
                const double tmp[9] = {r, k, length, r_eq2D, static_cast<double>(std::abs(is[ax0])), static_cast<double>(std::abs(is[ax1])), static_cast<double>(std::abs(is[ax2])), static_cast<double>(GN), hfact[ax0]};

                // green function value
                // Implementation note: having a 'if' in a loop is highly discouraged... however, this is the init so we prefer having a
                // this routine with a high readability and lower efficency than the opposite.
                if (r <= (hfact[ax0] + hfact[ax1] + hfact[ax2]) * .2) {
                    // we should enter this case for 2d and 3d cases
                    green[id + i0 * topo->nf()] = Gr0(tmp, Gdata);
                } else if (k <= (kfact[ax0] + kfact[ax1] + kfact[ax2]) * 0.2) {
                    // we should always enter this routine for 2d case and sometimes for 3d cases
                    green[id + i0 * topo->nf()] = Gk0(tmp, Gdata);
                } else {
                    green[id + i0 * topo->nf()] = G(tmp, Gdata);
                }
            }
        }
    }
    
    // reset the value in x=y=0.0 and k=0 for singular expressions
    if ((typeGreen == CHAT_2) && istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0) {
        // green[0] = -2.0 * log(1 + sqrt(2)) * c_1opiE3o2 / r_eq2D;
        green[0] = - 0.25 * c_1o2pi * (M_PI - 6.0 + 2.0 * log(0.5 * M_PI * r_eq2D));  //caution: mistake in [Chatelain2010]
    }
    END_FUNC;
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
 * @param length the characteristic length (only used for HEJ kernels = epsilon)
 */
void cmpt_Green_1dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length) {
    BEGIN_FUNC;

    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[0] != hfact[0], "grid spacing[0] cannot be = to dk[0]");
    FLUPS_CHECK(kfact[1] != hfact[1], "grid spacing[1] cannot be = to dk[1]");
    // check that if hfact or kfact != 0, they are not the same
    FLUPS_CHECK(!(kfact[2] == hfact[2] && (kfact[2]!= 0.0 || hfact[2] != 0.0)), "grid spacing[2] cannot be = to dk[2]");

    // @Todo For Helmolz, we need Green to be complex 
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.");
    // double* mygreen = green; //casting of the Green function to be able to access real and complex part

    GreenKernel G;   // the Green kernel (general expression in the whole domain)
    GreenKernel G0;  // the Green kernel (particular expression in k=0)

    switch (typeGreen) {
        case HEJ_2:
            G  = &hej_2_1unb2spe_;
            G0 = &hej_2_1unb2spe_k0_;
            break;
        case HEJ_4:
            G  = &hej_4_1unb2spe_;
            G0 = &hej_4_1unb2spe_k0_;
            break;
        case HEJ_6:
            G  = &hej_6_1unb2spe_;
            G0 = &hej_6_1unb2spe_k0_;
            break;
        case HEJ_8:
            G  = &hej_8_1unb2spe_;
            G0 = &hej_8_1unb2spe_k0_;
            break;
        case HEJ_10:
            G  = &hej_10_1unb2spe_;
            G0 = &hej_10_1unb2spe_k0_;
            break;        
        case HEJ_0:
            FLUPS_CHECK(false, "HEJ0 kernel not available for 1D unbounded problems.");
            break;        
        case CHAT_2:
            G  = &chat_2_1unb2spe_;
            G0 = &chat_2_1unb2spe_k0_;
            break;
        case LGF_2:
            // FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            G  = &lgf_2_1unb2spe_;
            G0 = &lgf_2_1unb2spe_;
            break;
        case LGF_4:
            G  = &lgf_4_1unb2spe_;
            G0 = &lgf_4_1unb2spe_;
            break;
        case LGF_6:
            G  = &lgf_6_1unb2spe_;
            G0 = &lgf_6_1unb2spe_;
            break;
        case LGF_8:
            G  = &lgf_6_1unb2spe_;
            G0 = &lgf_6_1unb2spe_;
            break;
        case MEHR_4L:
            G  = &mehr_4l_1unb2spe_;
            G0 = &mehr_4l_1unb2spe_;
            break;
        case MEHR_6L:
            G  = &mehr_6l_1unb2spe_;
            G0 = &mehr_6l_1unb2spe_;
            break;
            break;
        case MEHR_4F:
            FLUPS_CHECK(false, "MEHR_4F kernel not available for 1 dir unbounded problems.");
            // G  = &mehr_4f_1unb2spe_;
            // G0 = &mehr_4f_1unb2spe_;
            break;
        case MEHR_6F:
            FLUPS_CHECK(false, "MEHR_6F kernel not available for 1 dir unbounded problems.");
            // G  = &mehr_6f_1unb2spe_;
            // G0 = &mehr_6f_1unb2spe_;
            break;
        default:
            FLUPS_CHECK(false, "Green Function type unknown.");
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
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                int is[3];
                cmpt_symID(ax0,i0,i1,i2,istart,symstart,0,is);

                // (symmetrized) wave number : only 1 kfact is zero
                const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];
                const double k  = sqrt(k0 * k0 + k1 * k1 + k2 * k2);

                //(symmetrized) position : only 1 hfact is non-zero
                const double x0 = (is[ax0]) * hfact[ax0];
                const double x1 = (is[ax1]) * hfact[ax1];
                const double x2 = (is[ax2]) * hfact[ax2];
                const double r  = sqrt(x0 * x0 + x1 * x1 + x2 * x2);
                const double h = hfact[ax0] + hfact[ax1] + hfact[ax2];

                // GF parameters: k0, k1, k2, and h are for LGFs only
                const double tmp[7] = {r, k, length, k0, k1, k2, h};

                // green function value
                // Implementation note: having a 'if' in a loop is highly discouraged... however, this is the init so we prefer having a
                // this routine with a high readability and lower efficency than the opposite.
                if (k <= (kfact[ax0] + kfact[ax1] + kfact[ax2]) * 0.2) {
                    green[id + i0 * nf] = G0(tmp,NULL);
                }
                else{
                    green[id + i0 * nf] = G(tmp,NULL);
                }
            }
        }
    }
    END_FUNC;
}

/**
 * @brief Compute the Green function for 3dirspectral (in the whole spectral domain)
 * 
 * _Note_ on performance__: obviously, the Green function in full spectral
 * is \f$\-frac{1}{k^2}\f$ (at least for CHAT_2). We could perform that operation directly in the
 * loop of `dothemagic`. We here choose to still precompute and store it.
 * We burn more memory, but we should fasten `dothemagic` as we replace a
 * (expensive) evaluation of \f$\frac{1}{k^2}\f$ by a memory access.
 *
 * The wave number in each direction is obtained as k_i = (i_s + koffset_i) * kfact_i, where is the global (potentially symmetric) index.
 * 
 * @param topo the topology associated to the Green's function
 * @param hgrid the grid spacing h = hx = hy = hz, used only for the LGF
 * @param kfact the k multiplicative factor
 * @param koffset the k additive factor
 * @param symstart index of the symmetry in each direction
 * @param green the Green function array
 * @param typeGreen the type of Green function 
 * @param length the characteristic length (only used for HEJ kernels = epsilon)
 */
void cmpt_Green_0dirunbounded(const Topology *topo, const double hgrid, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length) {
    cmpt_Green_0dirunbounded(topo, hgrid, kfact, koffset, symstart, green, typeGreen, length, NULL, NULL);
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
 * @param length the characteristic length (only used for HEJ kernels = epsilon)
 * @param istart_custom global index where we start to fill data, in each dir. If NULL, we start at the beginning of the spectral space.
 * @param iend_custom global index where we end to fill data, in each dir. If NULL, we end at the end of the spectral space.
 */
void cmpt_Green_0dirunbounded(const Topology *topo, const double hgrid, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length, const int istart_custom[3], const int iend_custom[3]) {
    BEGIN_FUNC;

    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(kfact[0] != 0.0, "dk cannot be 0");
    FLUPS_CHECK(kfact[1] != 0.0, "dk cannot be 0");
    // FLUPS_CHECK(kfact[2] != 0.0, "dk cannot be 0");

    GreenKernel G;   // the Green kernel (general expression in the whole domain)

    switch (typeGreen) {
        case HEJ_2:
            G = &hej_2_0unb3spe_;
            break;
        case HEJ_4:
            G = &hej_4_0unb3spe_;
            break;
        case HEJ_6:
            G = &hej_6_0unb3spe_;
            break;
        case HEJ_8:
            G = &hej_8_0unb3spe_;
            break;
        case HEJ_10:
            G = &hej_10_0unb3spe_;
            break; 
        case HEJ_0:
            //spectral solution is here given by 1/k^2, i.e. same 
            //as CHAT_2 kernel
            G = &chat_2_0unb3spe_;
            break;                           
        case CHAT_2:
            G = &chat_2_0unb3spe_;
            break;
        case LGF_2:
            G = &lgf_2_0unb3spe_;
            break;
        case LGF_4:
            G = &lgf_4_0unb3spe_;
            break;
        case LGF_6:
            G = &lgf_6_0unb3spe_;
            break;
        case MEHR_4L:
            G = &mehr_4l_0unb3spe_;
            break;
        case MEHR_6L:
            G = &mehr_6l_0unb3spe_;
            break;
        case MEHR_4F:
            FLUPS_CHECK(false, "MEHR_4F not implemented in spectral form.");
            // G = &mehr_4f_0unb3spe_;
            break;
        case MEHR_6F:
            FLUPS_CHECK(false, "MEHR_6F not implemented in spectral form.");
            // G = &mehr_6f_0unb3spe_;
            break;
        default:
            FLUPS_CHECK(false, "Green Function type unknown.");
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
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);
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

                // const double tmp[2] = {ksqr, eps};
                const double tmp[6] = {ksqr, length, k0, k1, k2, hgrid};

                green[id + i0 * nf] = G(tmp, NULL);
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0 \
        && koffset[0]+koffset[1]+koffset[2]<0.2 ) {
        green[0] = 0.0;
    }
    END_FUNC;
}

