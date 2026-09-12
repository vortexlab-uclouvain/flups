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
#include "lgf_2unb1spe.hpp"
#include "hej_2unb1spe.hpp"

// Helper function to get the function pointer, depending on if it's on the GPU or not
inline void get_kernel_ptr(GreenKernel *host_ptr, GreenKernel *device_ptr, stream_t stream) {
#if defined(BACKEND_CUDA)
    m_cudart_check_errors(cudaMemcpyFromSymbolAsync(host_ptr, *device_ptr, sizeof(GreenKernel), 0, cudaMemcpyDeviceToHost, stream));
#elif defined(BACKEND_HIP)
    m_hiprt_check_errors(hipMemcpyFromSymbolAsync(host_ptr, *device_ptr, sizeof(GreenKernel), 0, hipMemcpyDeviceToHost, stream));
#else
    *host_ptr = *device_ptr;
#endif
}

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
void cmpt_Green_3dirunbounded(const Topology *topo, const double hfact[3], const double symstart[3], double *green, GreenType typeGreen, const double length, stream_t stream) {
    BEGIN_FUNC;
    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(hfact[0] != 0.0, "hfact[0] cannot be 0");
    FLUPS_CHECK(hfact[1] != 0.0, "hfact[1] cannot be 0");
    FLUPS_CHECK(hfact[2] != 0.0, "hfact[2] cannot be 0");

    GreenKernel G;

    double  G0;  // value of G in 0
    int     GN    = 0;
    double *Gdata = NULL;

    //==========================    3D  =================================
    switch (typeGreen) {
        case HEJ_2:
            get_kernel_ptr(&G, &hej_2_3unb0spe, stream);
            G0 = -M_SQRT2 / (4.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_4:
            get_kernel_ptr(&G, &hej_4_3unb0spe, stream);
            G0 = -3.0 * M_SQRT2 / (8.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_6:
            get_kernel_ptr(&G, &hej_6_3unb0spe, stream);
            G0 = -15.0 * M_SQRT2 / (32.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_8:
            get_kernel_ptr(&G, &hej_8_3unb0spe, stream);
            G0 = -35.0 * M_SQRT2 / (64.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_10:
            get_kernel_ptr(&G, &hej_10_3unb0spe, stream);
            G0 = -315.0 * M_SQRT2 / (512.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_0:
            get_kernel_ptr(&G, &hej_0_3unb0spe, stream);
            G0 = -1.0 / (2.0 * M_PI * M_PI * length);
            break;
        case CHAT_2:
            get_kernel_ptr(&G, &chat_2_3unb0spe, stream);
            G0 = -0.5 * pow(1.5 * c_1o2pi * hfact[0] * hfact[1] * hfact[2], 2. / 3.) / (hfact[0] * hfact[1] * hfact[2]);
            break;
        case LGF_2:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_2, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &lgf_2_3unb0spe, stream);
            break;
        case LGF_4:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_4, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &lgf_4_3unb0spe, stream);
            break;
        case LGF_6:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_6, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &lgf_6_3unb0spe, stream);
            break;
        case LGF_8:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_8, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &lgf_8_3unb0spe, stream);
            break;
        case MEHR_4L:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_4L, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &mehr_4l_3unb0spe, stream);
            break;
        case MEHR_6L:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_6L, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &mehr_6l_3unb0spe, stream);
            break;
        case MEHR_4F:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_4F, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &mehr_4f_3unb0spe, stream);
            break;
        case MEHR_6F:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_6F, 3, &GN, &Gdata, stream);
            get_kernel_ptr(&G, &mehr_6f_3unb0spe, stream);
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
    const int nloc[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)};
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};

    // redefine kfact, koffset & symstart for the capture by value
    const double hfact_[3]    = {hfact[0], hfact[1], hfact[2]};
    const double symstart_[3] = {symstart[0], symstart[1], symstart[2]};

    auto op = m_backend_lambda(const int idx) {
        // unroll the index
        const int i0 = idx % nloc[ax0];
        const int i1 = (idx / nloc[ax0]) % nloc[ax1];
        const int i2 = idx / (nloc[ax0] * nloc[ax1]);

        // local index start
        const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

        int is[3];
        cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, is);

        // symmetrized position
        const double x0 = (is[ax0]) * hfact_[ax0];
        const double x1 = (is[ax1]) * hfact_[ax1];
        const double x2 = (is[ax2]) * hfact_[ax2];

        // green function value
        const double r2 = x0 * x0 + x1 * x1 + x2 * x2;
        const double r  = sqrt(r2);

        // the first two arguments are used in standard kernels, the two zeros are for compatibility with the 2dirunbounded function,
        // and the others 5 ones are aimed for LGFs only
        // the symmetrized indexes will be negative!!
        const double tmp[9] = {r, length, 0, 0, static_cast<double>(std::abs(is[ax0])), static_cast<double>(std::abs(is[ax1])), static_cast<double>(std::abs(is[ax2])), static_cast<double>(GN), hfact_[ax0]};
        green[id + i0 * nf] = G(tmp, Gdata);
    };
    parallel_for(op, nloc[ax0] * nloc[ax1] * nloc[ax2], "3unb green", stream);

    // reset the value in 0.0 but not for LGF's since we have already pre-computed its value
    const bool is_lgf         = (typeGreen == LGF_2 || typeGreen == LGF_4 || typeGreen == LGF_6 || typeGreen == LGF_8);
    const bool is_mehr        = (typeGreen == MEHR_4F || typeGreen == MEHR_6F || typeGreen == MEHR_4L || typeGreen == MEHR_6L);
    const bool correct_istart = (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0);
    if (!is_lgf && !is_mehr && correct_istart) {
        auto op = m_backend_lambda(const int idx) {
            green[idx] = G0;
        };
        parallel_for(op, 1, "3unb green reset", stream);
    }
    // free Gdata if needed
    if (Gdata != NULL) {
        backend_free(Gdata);
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
 * @param hgrid the grid spacing in each direction (unlike hfact, it is also set in the spectral dir)
 * @param green the Green function array
 * @param typeGreen  the type of Green function
 * @param length the characteristic length (only used for HEJ kernels = epsilon)
 *
 * @warning For the HEJ kernels in 3D: According to [Spietz2018], we can obtain the **approximate** Green kernel
            by using the 2D unbounded kernel for mode 0 in the spectral direction, and the rest of the Green kernel
            is the same as in full spectral. We here fill with zero most part of Green data. Indeed, we are
            interested only in doing the FFT of hej__*2unb1spe_k0_ in the 2 remaining spatial directions. We will
            complete the Green function with the full spectral part afterwards, while going through
            Solver::cmptGreenFunction_.
 *
 * @note For the LGF and MEHR kernels in 3D the treatment is **exact**: the k = 0 plane still gets the shipped 2D
         lattice Green's function, but every other plane is filled here, in physical space, with the true screened
         slice LGF of that plane (see @ref lgf_2unb1spe.hpp), so that the usual Hockney-Eastwood mirror + FFT that
         follows produces the exact kernel. Solver::cmptGreenFunction_ must therefore *not* overwrite those planes
         with the sampled symbol (which would periodise the slice kernel over the doubled box).
 *
 */
void cmpt_Green_2dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], const double hgrid[3], double *green, GreenType typeGreen, const double length, stream_t stream) {
    BEGIN_FUNC;

    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[0] != hfact[0], "grid spacing[0] cannot be = to dk[0]");
    FLUPS_CHECK(kfact[1] != hfact[1], "grid spacing[1] cannot be = to dk[1]");
    // check that if hfact or kfact != 0, they are not the same
    FLUPS_CHECK(!(kfact[2] == hfact[2] && (kfact[2] != 0.0 || hfact[2] != 0.0)), "grid spacing[2] cannot be = to dk[2]");

    // @Todo For Helmolz, we need Green to be complex
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.");
    // opt_double_ptr mygreen = green; //casting of the Green function to be able to access real and complex part
    // Implementation note: if you want to do Helmolz, you need Hankel functions (3rd order Bessel) which are not implemented in stdC. Consider the use of boost lib.
    // notice that bessel_k has been introduced in c++17

    // check for isotropy
    const bool unb_dirs_are_isotropic = ((hfact[0] == 0.0) && (hfact[1] == hfact[2])) ||
                                        ((hfact[1] == 0.0) && (hfact[0] == hfact[2])) ||
                                        ((hfact[2] == 0.0) && (hfact[0] == hfact[1]));

    GreenKernel G;    // the Green kernel (general expression in the whole domain)
    GreenKernel Gk0;  // the Green kernel (particular expression in k=0)
    GreenKernel Gr0;  // the Green kernel (particular expression in r=0)

    int     GN    = 0;
    double *Gdata = NULL;

    switch (typeGreen) {
        case HEJ_2:
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &hej_2_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &hej_2_2unb1spe_r0, stream);
            break;
        case HEJ_4:
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &hej_4_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &hej_4_2unb1spe_r0, stream);
            break;
        case HEJ_6:
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &hej_6_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &hej_6_2unb1spe_r0, stream);
            break;
        case HEJ_8:
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &hej_8_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &hej_8_2unb1spe_r0, stream);
            break;
        case HEJ_10:
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &hej_10_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &hej_10_2unb1spe_r0, stream);
            break;
        case HEJ_0:
            FLUPS_WARNING("HEJ0 (theoretically spectral) kernel for 2D unbounded entails an approximation greatly affecting accuracy.");
            parallel_for(m_backend_lambda(int i) { init_Ji0(); }, 1, "init Ji0", stream);
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &hej_0_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &hej_0_2unb1spe_k0, stream);
            break;
        case CHAT_2:
            get_kernel_ptr(&G, &chat_2_2unb1spe, stream);
            get_kernel_ptr(&Gk0, &chat_2_2unb1spe_k0, stream);
            get_kernel_ptr(&Gr0, &chat_2_2unb1spe_r0, stream);
            // caution: the value of G in k=r=0 is specified at the end of this routine
            break;
        case LGF_2:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_2, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &lgf_2_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &lgf_2_2unb0spe, stream);
            break;
        case LGF_4:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_4, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &lgf_4_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &lgf_4_2unb0spe, stream);
            break;
        case LGF_6:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_6, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &lgf_6_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &lgf_6_2unb0spe, stream);
            break;
        case LGF_8:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_8, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &lgf_8_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &lgf_8_2unb0spe, stream);
            break;
        case MEHR_4L:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_4L, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &mehr_4l6l_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &mehr_4l6l_2unb0spe, stream);
            break;
        case MEHR_6L:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_6L, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &mehr_4l6l_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &mehr_4l6l_2unb0spe, stream);
            break;
        case MEHR_4F:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_4F, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &mehr_4f_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &mehr_4f_2unb0spe, stream);
            break;
        case MEHR_6F:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_6F, 2, &GN, &Gdata, stream);
            // associate the Green's function
            get_kernel_ptr(&G, &zero, stream);
            get_kernel_ptr(&Gk0, &mehr_6f_2unb0spe, stream);
            get_kernel_ptr(&Gr0, &mehr_6f_2unb0spe, stream);
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
    const int    nloc[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)};
    const int    nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
    const double r_eq2D  = c_1osqrtpi * sqrt(hfact[ax0] * hfact[ax1] + hfact[ax1] * hfact[ax2] + hfact[ax2] * hfact[ax0]);

    // redefine kfact, koffset & symstart for the capture by value
    const double hfact_[3]    = {hfact[0], hfact[1], hfact[2]};
    const double kfact_[3]    = {kfact[0], kfact[1], kfact[2]};
    const double koffset_[3]  = {koffset[0], koffset[1], koffset[2]};
    const double symstart_[3] = {symstart[0], symstart[1], symstart[2]};

    //-------------------------------------------------------------------------
    /** - find the spectral direction, if any
     *
     * It is the only one with a non-zero kfact. A pure 2D case (the third direction being EMPTY)
     * has none: every point is then a k = 0 point and the standard loop below is used as is.
     */
    //-------------------------------------------------------------------------
    int ax_spec = -1;
    int n_spec  = 0;
    for (int id = 0; id < 3; ++id) {
        if (kfact[id] != 0.0) {
            ax_spec = id;
            n_spec += 1;
        }
    }
    FLUPS_CHECK(n_spec <= 1, "cmpt_Green_2dirunbounded got %d spectral directions", n_spec);

    // the LGF/MEHR kernels get the exact screened slice LGF on every k != 0 plane, the HEJ ones the
    // exact screened slice of their regularised kernel; both are filled in physical space below
    const bool is_lgf_slice   = LGFTwoUnbounded::handles(typeGreen);
    const bool is_hej_slice   = HEJTwoUnbounded::handles(typeGreen);
    const bool do_exact_slice = (is_lgf_slice || is_hej_slice) && (n_spec == 1);

    if (!do_exact_slice) {
        auto op = m_backend_lambda(const int idx) {
            // unroll the index
            const int i0 = idx % nloc[ax0];
            const int i1 = (idx / nloc[ax0]) % nloc[ax1];
            const int i2 = idx / (nloc[ax0] * nloc[ax1]);

            // local index start
            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

            int is[3];
            cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, is);

            // (symmetrized) wave number : only one kfact is non-zero
            const double k0 = (is[ax0] + koffset_[ax0]) * kfact_[ax0];
            const double k1 = (is[ax1] + koffset_[ax1]) * kfact_[ax1];
            const double k2 = (is[ax2] + koffset_[ax2]) * kfact_[ax2];
            const double k  = k0 + k1 + k2;

            //(symmetrized) position : only one hfact is zero
            const double x0 = (is[ax0]) * hfact_[ax0];
            const double x1 = (is[ax1]) * hfact_[ax1];
            const double x2 = (is[ax2]) * hfact_[ax2];
            const double r  = sqrt(x0 * x0 + x1 * x1 + x2 * x2);

            // the symmetrized indexes will be negative!!
            const double tmp[9] = {r, k, length, r_eq2D, static_cast<double>(std::abs(is[ax0])), static_cast<double>(abs(is[ax1])), static_cast<double>(abs(is[ax2])), static_cast<double>(GN), hfact_[ax0]};

            // green function value
            // Implementation note: having a 'if' in a loop is highly discouraged... however, this is the init so we prefer having a
            // this routine with a high readability and lower efficency than the opposite.
            if (r <= (hfact_[ax0] + hfact_[ax1] + hfact_[ax2]) * .2) {
                // we should enter this case for 2d and 3d cases
                green[id + i0 * nf] = Gr0(tmp, Gdata);
            } else if (k <= (kfact_[ax0] + kfact_[ax1] + kfact_[ax2]) * 0.2) {
                // we should always enter this routine for 2d case and sometimes for 3d cases
                green[id + i0 * nf] = Gk0(tmp, Gdata);
            } else {
                green[id + i0 * nf] = G(tmp, Gdata);
            }
        };
        parallel_for(op, nloc[ax0] * nloc[ax1] * nloc[ax2], "2unb green", stream);
    } else {
        //---------------------------------------------------------------------
        /** - LGF/MEHR: the k = 0 plane keeps the 2D table, every other plane is filled below with
         *    the exact screened slice LGF, see @ref lgf_2unb1spe.hpp
         */
        //---------------------------------------------------------------------
        const int    ax_s     = ax_spec;
        const double kfact_s  = kfact[ax_s];
        const double koff_s   = koffset[ax_s];
        const double symsta_s = symstart[ax_s];

        // ------------------------------------------------------------------ the k ~ 0 plane
        const bool hej_k0 = is_hej_slice;
        auto op_k0 = m_backend_lambda(const int idx) {
            const int i0 = idx % nloc[ax0];
            const int i1 = (idx / nloc[ax0]) % nloc[ax1];
            const int i2 = idx / (nloc[ax0] * nloc[ax1]);

            const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

            int is[3];
            cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, is);

            // the folded spectral index may be negative, the kernel is even in k
            const double k = (is[ax_s] + koff_s) * kfact_s;

            double value = 0.0;
            if (fabs(k) <= kfact_s * 0.2) {
                // the two-dimensional Green's function of that plane.
                // Note: for the LGF/MEHR kernels we do NOT test r ~ 0 first, because the r = 0
                // point of a k != 0 slice is a regular point of the slice table and the r = 0
                // point of the k = 0 plane is handled by Gk0 itself. The regularised k = 0
                // kernels carry a log(r) and do need their r ~ 0 companion.
                const double x0     = (is[ax0]) * hfact_[ax0];
                const double x1     = (is[ax1]) * hfact_[ax1];
                const double x2     = (is[ax2]) * hfact_[ax2];
                const double r      = sqrt(x0 * x0 + x1 * x1 + x2 * x2);
                const double tmp[9] = {r, k, length, r_eq2D, static_cast<double>(std::abs(is[ax0])), static_cast<double>(abs(is[ax1])), static_cast<double>(abs(is[ax2])), static_cast<double>(GN), hfact_[ax0]};
                const bool   at_r0  = (r <= (hfact_[ax0] + hfact_[ax1] + hfact_[ax2]) * 0.2);
                value               = (hej_k0 && at_r0) ? Gr0(tmp, Gdata) : Gk0(tmp, Gdata);
            }
            green[id + i0 * nf] = value;
        };
        parallel_for(op_k0, nloc[ax0] * nloc[ax1] * nloc[ax2], "2unb green k0 plane", stream);

        // ------------------------------------------------------------------ the k != 0 planes
        // the loop levels: level 0 is ax0 (the fastest one), level 2 is ax2
        const int axl[3] = {ax0, ax1, ax2};
        int       lvl_s  = 0;
        for (int l = 0; l < 3; ++l) {
            if (axl[l] == ax_s) {
                lvl_s = l;
            }
        }
        // keep the fastest remaining level first, so that the writes (and the reads in R) stay
        // as contiguous as the slice allows
        const int lvl_a = (lvl_s == 0) ? 1 : 0;
        const int lvl_b = (lvl_s == 2) ? 1 : 2;
        const int nlvl[3] = {nloc[ax0], nloc[ax1], nloc[ax2]};
        const int na      = nlvl[lvl_a];
        const int nb      = nlvl[lvl_b];

        // the unbounded direction handled by the closed form (x) and by the quadrature (y)
        const int ax_x = axl[lvl_a];
        const int ax_y = axl[lvl_b];

        // the largest |n| reachable in each unbounded direction: symstart when the direction is
        // folded (UNB-UNB), the full extent otherwise (UNB-EVEN/ODD, which is not folded)
        const int nglob[3] = {topo->nglob(0), topo->nglob(1), topo->nglob(2)};
        const int nx_max   = (symstart[ax_x] > 0.0) ? (int)symstart[ax_x] : (nglob[ax_x] - 1);
        const int ny_max   = (symstart[ax_y] > 0.0) ? (int)symstart[ax_y] : (nglob[ax_y] - 1);

        // the lattice closed forms assume an isotropic lattice in the three directions; the HEJ ones
        // are continuum kernels and take the physical position and wavenumber, so they do not
        FLUPS_CHECK(!is_lgf_slice || (hgrid[0] == hgrid[1] && hgrid[1] == hgrid[2]),
                    "the grid has to be isotropic to use the LGF/MEHR kernels in 2 unbounded + 1 spectral directions (h = %e %e %e)", hgrid[0], hgrid[1], hgrid[2]);
        const double h_spec = hgrid[ax_s];
        const double h_x    = hgrid[ax_x];
        const double h_y    = hgrid[ax_y];

        //---------------------------------------------------------------------
        /** - stage 1: the Gauss-Legendre rule and the Cw table, built once */
        //---------------------------------------------------------------------
        // the smallest non-zero |k| of the grid: |is + koffset| is at least 1 (koffset = 0) or
        // koffset itself (koffset = 0.5), and the quadrature must resolve that slice too
        const double kzmin   = kfact_s * ((koff_s == 0.0) ? 1.0 : koff_s);
        const double kzmin_h = kzmin * h_spec;
        // the lattice rule integrates over the Brillouin zone [0, pi] in lattice units, the
        // regularised one over [0, kmax] in physical units
        const double hej_kmax = is_hej_slice ? HEJTwoUnbounded::quad_kmax(length) : 0.0;
        const int    nq       = is_hej_slice ? HEJTwoUnbounded::quad_size(ny_max, h_y, hej_kmax, length, kzmin)
                                             : LGFTwoUnbounded::quad_size(ny_max, kzmin_h);
        const int nnx = nx_max + 1;
        const int nny = ny_max + 1;

        double *kq_host = (double *)m_calloc(sizeof(double) * nq);
        double *wq_host = (double *)m_calloc(sizeof(double) * nq);
        if (is_hej_slice) {
            HEJTwoUnbounded::gauss_legendre_0kmax(nq, hej_kmax, kq_host, wq_host);
        } else {
            LGFTwoUnbounded::gauss_legendre_0pi(nq, kq_host, wq_host);
        }

        double *kq = (double *)backend_malloc(sizeof(double) * nq);
        double *wq = (double *)backend_malloc(sizeof(double) * nq);
        backend_memcpy(kq, kq_host, sizeof(double) * nq, memcpyH2D, stream);
        backend_memcpy(wq, wq_host, sizeof(double) * nq, memcpyH2D, stream);
        backend_sync_stream(stream);
        m_free(kq_host);
        m_free(wq_host);

        double *cw = (double *)backend_malloc(sizeof(double) * nny * nq);
        double *rr = (double *)backend_malloc(sizeof(double) * nnx * nq);
        if (is_hej_slice) {
            HEJTwoUnbounded::build_cw(cw, kq, wq, ny_max, nq, h_y, stream);
        } else {
            LGFTwoUnbounded::build_cw(cw, kq, wq, ny_max, nq, stream);
        }

        //---------------------------------------------------------------------
        /** - stage 2 & 3: for each locally owned k != 0 plane, rebuild R and fill the plane */
        //---------------------------------------------------------------------
        for (int isl = 0; isl < nlvl[lvl_s]; ++isl) {
            // the global, folded, index of that plane in the spectral direction
            const int ie = istart[ax_s] + isl;
            const int is_s = (symsta_s == 0.0 || ie <= symsta_s) ? ie : -std::max((int)fabs(2.0 * symsta_s - ie), 1);
            const double k = (is_s + koff_s) * kfact_s;
            if (fabs(k) <= kfact_s * 0.2) {
                continue;  // the k = 0 plane has already been filled with the 2D table
            }
            // the closed form is even in every wavenumber, so only |k| matters
            if (is_hej_slice) {
                HEJTwoUnbounded::build_r(rr, kq, (int)typeGreen, nx_max, nq, fabs(k), h_x, length, stream);
            } else {
                LGFTwoUnbounded::build_r(rr, kq, (int)typeGreen, nx_max, nq, fabs(k) * h_spec, stream);
            }

            const int  isl_      = isl;
            const bool hej_slice = is_hej_slice;
            auto op_slice  = m_backend_lambda(const int idx) {
                int ii[3];
                ii[lvl_s] = isl_;
                ii[lvl_a] = idx % na;
                ii[lvl_b] = idx / na;
                const int i0 = ii[0];
                const int i1 = ii[1];
                const int i2 = ii[2];

                const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

                int is[3];
                cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, is);
                const int nx = abs(is[ax_x]);
                const int ny = abs(is[ax_y]);
                assert(nx <= nx_max && ny <= ny_max);

                // the lattice tables are in the positive-operator convention and FLUPS stores -G;
                // the regularised ones already carry the stored sign
                green[id + i0 * nf] = hej_slice ? HEJTwoUnbounded::slice_point(rr, cw, nx, ny, nnx, nny, nq)
                                                : -LGFTwoUnbounded::slice_point(rr, cw, nx, ny, nnx, nny, nq);
            };
            parallel_for(op_slice, na * nb, "2unb1spe green slice", stream);
        }

        backend_sync_stream(stream);
        backend_free(kq);
        backend_free(wq);
        backend_free(cw);
        backend_free(rr);
    }

    // reset the value in x=y=0.0 and k=0 for singular expressions
    if ((typeGreen == CHAT_2) && istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0) {
        auto op = m_backend_lambda(const int idx) {
            // caution: mistake in [Chatelain2010]
            green[idx] = 0.25 * c_1o2pi * (M_PI - 6.0 + 2.0 * log(0.5 * M_PI * r_eq2D * r_eq2D));
        };
        parallel_for(op, 1, "2unb green reset", stream);
    }

    // free Gdata if needed
    if (Gdata != NULL) {
        backend_free(Gdata);
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
void cmpt_Green_1dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length, stream_t stream) {
    BEGIN_FUNC;

    // assert that the green spacing and dk is not 0.0 - this is also a way to check that ax0 will be spectral, and the others are still to be transformed
    FLUPS_CHECK(kfact[0] != hfact[0], "grid spacing[0] cannot be = to dk[0]");
    FLUPS_CHECK(kfact[1] != hfact[1], "grid spacing[1] cannot be = to dk[1]");
    // check that if hfact or kfact != 0, they are not the same
    FLUPS_CHECK(!(kfact[2] == hfact[2] && (kfact[2] != 0.0 || hfact[2] != 0.0)), "grid spacing[2] cannot be = to dk[2]");

    // @Todo For Helmolz, we need Green to be complex
    // FLUPS_CHECK(topo->isComplex(), "I can't fill a non complex topo with a complex green function.");
    // double* mygreen = green; //casting of the Green function to be able to access real and complex part

    GreenKernel G;   // the Green kernel (general expression in the whole domain)
    GreenKernel G0;  // the Green kernel (particular expression in k=0)

    switch (typeGreen) {
        case HEJ_2:
            get_kernel_ptr(&G, &hej_2_1unb2spe, stream);
            get_kernel_ptr(&G0, &hej_2_1unb2spe_k0, stream);
            break;
        case HEJ_4:
            get_kernel_ptr(&G, &hej_4_1unb2spe, stream);
            get_kernel_ptr(&G0, &hej_4_1unb2spe_k0, stream);
            break;
        case HEJ_6:
            get_kernel_ptr(&G, &hej_6_1unb2spe, stream);
            get_kernel_ptr(&G0, &hej_6_1unb2spe_k0, stream);
            break;
        case HEJ_8:
            get_kernel_ptr(&G, &hej_8_1unb2spe, stream);
            get_kernel_ptr(&G0, &hej_8_1unb2spe_k0, stream);
            break;
        case HEJ_10:
            get_kernel_ptr(&G, &hej_10_1unb2spe, stream);
            get_kernel_ptr(&G0, &hej_10_1unb2spe_k0, stream);
            break;
        case HEJ_0:
            FLUPS_CHECK(false, "HEJ0 kernel not available for 1D unbounded problems.");
            break;
        case CHAT_2:
            get_kernel_ptr(&G, &chat_2_1unb2spe, stream);
            get_kernel_ptr(&G0, &chat_2_1unb2spe_k0, stream);
            break;
        case LGF_2:
            get_kernel_ptr(&G, &lgf_2_1unb2spe, stream);
            get_kernel_ptr(&G0, &lgf_2_1unb2spe, stream);
            break;
        case LGF_4:
            get_kernel_ptr(&G, &lgf_4_1unb2spe, stream);
            get_kernel_ptr(&G0, &lgf_4_1unb2spe, stream);
            break;
        case LGF_6:
            get_kernel_ptr(&G, &lgf_6_1unb2spe, stream);
            get_kernel_ptr(&G0, &lgf_6_1unb2spe, stream);
            break;
        case LGF_8:
            get_kernel_ptr(&G, &lgf_8_1unb2spe, stream);
            get_kernel_ptr(&G0, &lgf_8_1unb2spe, stream);
            break;
        case MEHR_4L:
            get_kernel_ptr(&G, &mehr_4l_1unb2spe, stream);
            get_kernel_ptr(&G0, &mehr_4l_1unb2spe, stream);
            break;
        case MEHR_6L:
            get_kernel_ptr(&G, &mehr_6l_1unb2spe, stream);
            get_kernel_ptr(&G0, &mehr_6l_1unb2spe, stream);
            break;
        case MEHR_4F:
            get_kernel_ptr(&G, &mehr_4f_1unb2spe, stream);
            get_kernel_ptr(&G0, &mehr_4f_1unb2spe, stream);
            break;
        case MEHR_6F:
            get_kernel_ptr(&G, &mehr_6f_1unb2spe, stream);
            get_kernel_ptr(&G0, &mehr_6f_1unb2spe, stream);
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
    const int nloc[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)};
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};

    // redefine kfact, koffset & symstart for the capture by value
    const double hfact_[3]    = {hfact[0], hfact[1], hfact[2]};
    const double kfact_[3]    = {kfact[0], kfact[1], kfact[2]};
    const double koffset_[3]  = {koffset[0], koffset[1], koffset[2]};
    const double symstart_[3] = {symstart[0], symstart[1], symstart[2]};

    auto op = m_backend_lambda(const int idx) {
        // unroll the index
        const int i0 = idx % nloc[ax0];
        const int i1 = (idx / nloc[ax0]) % nloc[ax1];
        const int i2 = idx / (nloc[ax0] * nloc[ax1]);

        // local index start
        const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

        int is[3];
        cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, is);

        // (symmetrized) wave number : only 1 kfact is zero
        const double k0 = (is[ax0] + koffset_[ax0]) * kfact_[ax0];
        const double k1 = (is[ax1] + koffset_[ax1]) * kfact_[ax1];
        const double k2 = (is[ax2] + koffset_[ax2]) * kfact_[ax2];
        const double k  = sqrt(k0 * k0 + k1 * k1 + k2 * k2);

        //(symmetrized) position : only 1 hfact is non-zero
        const double x0 = (is[ax0]) * hfact_[ax0];
        const double x1 = (is[ax1]) * hfact_[ax1];
        const double x2 = (is[ax2]) * hfact_[ax2];
        const double r  = sqrt(x0 * x0 + x1 * x1 + x2 * x2);
        const double h  = hfact_[ax0] + hfact_[ax1] + hfact_[ax2];

        // GF parameters: k0, k1, k2, and h are for LGFs only
        const double tmp[7] = {r, k, length, k0, k1, k2, h};

        // green function value
        // Implementation note: having a 'if' in a loop is highly discouraged... however, this is the init so we prefer having a
        // this routine with a high readability and lower efficency than the opposite.
        if (k <= (kfact_[ax0] + kfact_[ax1] + kfact_[ax2]) * 0.2) {
            green[id + i0 * nf] = G0(tmp, NULL);
        } else {
            green[id + i0 * nf] = G(tmp, NULL);
        }
    };
    parallel_for(op, nloc[ax0] * nloc[ax1] * nloc[ax2], "1unb green", stream);
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
void cmpt_Green_0dirunbounded(const Topology *topo, const double hgrid, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length, stream_t stream) {
    cmpt_Green_0dirunbounded(topo, hgrid, kfact, koffset, symstart, green, typeGreen, length, NULL, NULL, stream);
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
void cmpt_Green_0dirunbounded(const Topology *topo, const double hgrid, const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length, const int istart_custom[3], const int iend_custom[3], stream_t stream) {
    BEGIN_FUNC;

    // assert that the green spacing is not 0.0 everywhere
    FLUPS_CHECK(kfact[0] != 0.0, "dk cannot be 0");
    FLUPS_CHECK(kfact[1] != 0.0, "dk cannot be 0");
    // FLUPS_CHECK(kfact[2] != 0.0, "dk cannot be 0");

    // Select the Green function to use
    GreenKernel G;
    switch (typeGreen) {
        case HEJ_2:
            get_kernel_ptr(&G, &hej_2_0unb3spe, stream);
            break;
        case HEJ_4:
            get_kernel_ptr(&G, &hej_4_0unb3spe, stream);
            break;
        case HEJ_6:
            get_kernel_ptr(&G, &hej_6_0unb3spe, stream);
            break;
        case HEJ_8:
            get_kernel_ptr(&G, &hej_8_0unb3spe, stream);
            break;
        case HEJ_10:
            get_kernel_ptr(&G, &hej_10_0unb3spe, stream);
            break;
        case HEJ_0:
            // spectral solution is here given by 1/k^2, i.e. same as the CHAT_2 kernel
            get_kernel_ptr(&G, &chat_2_0unb3spe, stream);
            break;
        case CHAT_2:
            get_kernel_ptr(&G, &chat_2_0unb3spe, stream);
            break;
        case LGF_2:
            get_kernel_ptr(&G, &lgf_2_0unb3spe, stream);
            break;
        case LGF_4:
            get_kernel_ptr(&G, &lgf_4_0unb3spe, stream);
            break;
        case LGF_6:
            get_kernel_ptr(&G, &lgf_6_0unb3spe, stream);
            break;
        case LGF_8:
            get_kernel_ptr(&G, &lgf_8_0unb3spe, stream);
            break;
        case MEHR_4L:
            get_kernel_ptr(&G, &mehr_4l_0unb3spe, stream);
            break;
        case MEHR_6L:
            get_kernel_ptr(&G, &mehr_6l_0unb3spe, stream);
            break;
        case MEHR_4F:
            get_kernel_ptr(&G, &mehr_4f_0unb3spe, stream);
            break;
        case MEHR_6F:
            get_kernel_ptr(&G, &mehr_6f_0unb3spe, stream);
            break;
        default:
            FLUPS_CHECK(false, "Green Function type unknown.");
    }

    int istart[3];
    topo->get_istart_glob(istart);

    int is_[3] = {0, 0, 0};
    int ie_[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)};

    if (istart_custom != NULL) {
        // switch to local index
        for (int ip = 0; ip < 3; ip++) {
            is_[ip] = fmin(fmax(istart_custom[ip] - istart[ip], 0), topo->nloc(ip));
        }
    }
    if (iend_custom != NULL) {
        // switch to local index
        for (int ip = 0; ip < 3; ip++) {
            ie_[ip] = fmin(fmax(iend_custom[ip] - istart[ip], 0), topo->nloc(ip));
        }
    }
    const int is[3] = {is_[0], is_[1], is_[2]};
    const int ie[3] = {ie_[0], ie_[1], ie_[2]};
    const int ni[3] = {ie[0] - is[0], ie[1] - is[1], ie[2] - is[2]};

    const int nf      = topo->nf();
    const int ax0     = topo->axis();
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};

    // redefine kfact, koffset & symstart for the capture by value
    const double kfact_[3]    = {kfact[0], kfact[1], kfact[2]};
    const double koffset_[3]  = {koffset[0], koffset[1], koffset[2]};
    const double symstart_[3] = {symstart[0], symstart[1], symstart[2]};

    auto op = m_backend_lambda(const int idx) {
        // unroll the index
        const int i0 = idx % ni[ax0] + is[ax0];
        const int i1 = (idx / ni[ax0]) % ni[ax1] + is[ax1];
        const int i2 = idx / (ni[ax0] * ni[ax1]) + is[ax2];

        // local index start
        const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

        int il[3];
        cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, il);

        // (symmetrized) wave number
        const double k0 = (il[ax0] + koffset_[ax0]) * kfact_[ax0];
        const double k1 = (il[ax1] + koffset_[ax1]) * kfact_[ax1];
        const double k2 = (il[ax2] + koffset_[ax2]) * kfact_[ax2];

        // green function value
        const double ksqr   = k0 * k0 + k1 * k1 + k2 * k2;
        const double tmp[6] = {ksqr, length, k0, k1, k2, hgrid};

        green[id + i0 * nf] = G(tmp, NULL);
    };
    parallel_for(op, ni[0] * ni[1] * ni[2], "0unb green", stream);

    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0 && koffset[0] + koffset[1] + koffset[2] < 0.2) {
        auto op = m_backend_lambda(const int idx) {
            green[idx] = 0.0;
        };
        parallel_for(op, 1, "0unb green reset", stream);
    }
    END_FUNC;
}
