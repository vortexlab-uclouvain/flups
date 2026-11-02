/**
 * @file green_functions.cpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#include "green_functions.hpp"

#include <algorithm>
#include <climits>

#include "green_kernels.hpp"
#include "ji0.hpp"
#include "hej_2unb1spe.hpp"
#include "lgf_2unb1spe.hpp"

/**
 * @brief reads the value of a kernel symbol into host memory
 *
 * On the GPU backends the kernels are `__device__` function pointers, whose value is known to the
 * device only: it has to be fetched through the symbol.
 *
 * @param symbol the address of the kernel symbol, as declared in @ref green_kernels.hpp
 * @return the kernel, usable from the host to build a device lambda
 */
static inline GreenKernel kernel_from_symbol(GreenKernel *symbol) {
#if defined(BACKEND_CUDA)
    GreenKernel kernel;
    m_cudart_check_errors(cudaMemcpyFromSymbol(&kernel, *symbol, sizeof(GreenKernel)));
    return kernel;
#elif defined(BACKEND_HIP)
    GreenKernel kernel;
    m_hiprt_check_errors(hipMemcpyFromSymbol(&kernel, *symbol, sizeof(GreenKernel)));
    return kernel;
#else
    return *symbol;
#endif
}

/**
 * @brief the interval of @f$ |i_s| @f$ reached by the local indexes of one direction
 *
 * @ref cmpt_symID folds a global index @p ie onto a signed one; the values of @f$ |i_s| @f$ taken
 * over a local range `[istart, istart + nloc)` form a **contiguous** interval, also when the range
 * straddles @p symstart (the folded branch decreases while the direct one increases). A table
 * indexed by `|i_s| - n_lo` therefore covers the rank exactly, with at most `nloc + 1` entries.
 *
 * @param istart the first global index owned by this rank
 * @param nloc the number of indexes owned by this rank
 * @param symstart the symmetry index of the direction, 0 if the direction is not folded
 * @param[out] n_lo the smallest @f$ |i_s| @f$ reached
 * @param[out] n_cnt the number of values of @f$ |i_s| @f$ reached, i.e. the table size
 */
static inline void cmpt_local_absID_(const int istart, const int nloc, const double symstart, int *n_lo, int *n_cnt) {
    int lo = INT_MAX;
    int hi = -1;
    for (int i = 0; i < nloc; ++i) {
        const int ie = istart + i;
        const int is = (symstart == 0.0 || ie <= symstart) ? ie : -std::max((int)fabs(2.0 * symstart - ie), 1);
        const int n  = std::abs(is);
        lo           = std::min(lo, n);
        hi           = std::max(hi, n);
    }
    if (hi < lo) {
        // no local point at all: keep a one-entry table so that the allocation stays valid
        lo = 0;
        hi = 0;
    }
    FLUPS_CHECK(hi - lo + 1 <= nloc + 1, "the local |i_s| interval [%d,%d] is not contiguous over %d indexes", lo, hi, nloc);
    *n_lo  = lo;
    *n_cnt = hi - lo + 1;
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
            G  = kernel_from_symbol(&hej_2_3unb0spe);
            G0 = -M_SQRT2 / (4.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_4:
            G  = kernel_from_symbol(&hej_4_3unb0spe);
            G0 = -3.0 * M_SQRT2 / (8.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_6:
            G  = kernel_from_symbol(&hej_6_3unb0spe);
            G0 = -15.0 * M_SQRT2 / (32.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_8:
            G  = kernel_from_symbol(&hej_8_3unb0spe);
            G0 = -35.0 * M_SQRT2 / (64.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_10:
            G  = kernel_from_symbol(&hej_10_3unb0spe);
            G0 = -315.0 * M_SQRT2 / (512.0 * length * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_0:
            G  = kernel_from_symbol(&hej_0_3unb0spe);
            G0 = -1.0 / (2.0 * M_PI * M_PI * length);
            break;
        case CHAT_2:
            G  = kernel_from_symbol(&chat_2_3unb0spe);
            G0 = -0.5 * pow(1.5 * c_1o2pi * hfact[0] * hfact[1] * hfact[2], 2. / 3.) / (hfact[0] * hfact[1] * hfact[2]);
            break;
        case LGF_2:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_2, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&lgf_2_3unb0spe);
            break;
        case LGF_4:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_4, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&lgf_4_3unb0spe);
            break;
        case LGF_6:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_6, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&lgf_6_3unb0spe);
            break;
        case LGF_8:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(LGF_8, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&lgf_8_3unb0spe);
            break;
        case MEHR_4L:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_4L, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&mehr_4l_3unb0spe);
            break;
        case MEHR_6L:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_6L, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&mehr_6l_3unb0spe);
            break;
        case MEHR_4F:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_4F, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&mehr_4f_3unb0spe);
            break;
        case MEHR_6F:
            FLUPS_CHECK((hfact[0] == hfact[1]) && (hfact[1] == hfact[2]), "the grid has to be isotropic to use the LGFs");
            lgf_readfile_(MEHR_6F, 3, &GN, &Gdata, stream);
            G = kernel_from_symbol(&mehr_6f_3unb0spe);
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
 * @param prof the profiler (may be NULL): the construction of the exact slices is timed as `green_kernel`
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
void cmpt_Green_2dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], const double hgrid[3], double *green, GreenType typeGreen, const double length, stream_t stream, H3LPR::Profiler *prof) {
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
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&hej_2_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&hej_2_2unb1spe_r0);
            break;
        case HEJ_4:
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&hej_4_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&hej_4_2unb1spe_r0);
            break;
        case HEJ_6:
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&hej_6_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&hej_6_2unb1spe_r0);
            break;
        case HEJ_8:
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&hej_8_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&hej_8_2unb1spe_r0);
            break;
        case HEJ_10:
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&hej_10_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&hej_10_2unb1spe_r0);
            break;
        case HEJ_0:
            FLUPS_WARNING("HEJ0 (theoretically spectral) kernel for 2D unbounded entails an approximation greatly affecting accuracy.");
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&hej_0_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&hej_0_2unb1spe_k0);
            break;
        case CHAT_2:
            G   = kernel_from_symbol(&chat_2_2unb1spe);
            Gk0 = kernel_from_symbol(&chat_2_2unb1spe_k0);
            Gr0 = kernel_from_symbol(&chat_2_2unb1spe_r0);
            // caution: the value of G in k=r=0 is specified at the end of this routine
            break;
        case LGF_2:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_2, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&lgf_2_2unb0spe);
            Gr0 = kernel_from_symbol(&lgf_2_2unb0spe);
            break;
        case LGF_4:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_4, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&lgf_4_2unb0spe);
            Gr0 = kernel_from_symbol(&lgf_4_2unb0spe);
            break;
        case LGF_6:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_6, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&lgf_6_2unb0spe);
            Gr0 = kernel_from_symbol(&lgf_6_2unb0spe);
            break;
        case LGF_8:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the LGFs");
            // read the LGF data and store it
            lgf_readfile_(LGF_8, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&lgf_8_2unb0spe);
            Gr0 = kernel_from_symbol(&lgf_8_2unb0spe);
            break;
        case MEHR_4L:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_4L, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&mehr_4l6l_2unb0spe);
            Gr0 = kernel_from_symbol(&mehr_4l6l_2unb0spe);
            break;
        case MEHR_6L:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_6L, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&mehr_4l6l_2unb0spe);
            Gr0 = kernel_from_symbol(&mehr_4l6l_2unb0spe);
            break;
        case MEHR_4F:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_4F, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&mehr_4f_2unb0spe);
            Gr0 = kernel_from_symbol(&mehr_4f_2unb0spe);
            break;
        case MEHR_6F:
            FLUPS_CHECK(unb_dirs_are_isotropic, "the grid has to be isotropic to use the MEHR kernels");
            // read the LGF data and store it
            lgf_readfile_(MEHR_6F, 2, &GN, &Gdata, stream);
            // associate the Green's function
            G   = kernel_from_symbol(&zero);
            Gk0 = kernel_from_symbol(&mehr_6f_2unb0spe);
            Gr0 = kernel_from_symbol(&mehr_6f_2unb0spe);
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
        // the construction of the slices is timed on its own: it is the cost of the kernel itself,
        // without the plans and the transforms that surround it
        backend_sync_stream(stream);
        m_profStarti(prof, "green_kernel");
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
        auto       op_k0  = m_backend_lambda(const int idx) {
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
        const int lvl_a   = (lvl_s == 0) ? 1 : 0;
        const int lvl_b   = (lvl_s == 2) ? 1 : 2;
        const int nlvl[3] = {nloc[ax0], nloc[ax1], nloc[ax2]};
        const int na      = nlvl[lvl_a];
        const int nb      = nlvl[lvl_b];

        // the largest |n| reachable in each unbounded direction: symstart when the direction is
        // folded (UNB-UNB), the full extent otherwise (UNB-EVEN/ODD, which is not folded)
        const int nglob[3] = {topo->nglob(0), topo->nglob(1), topo->nglob(2)};
        const int nmax_a   = (symstart[axl[lvl_a]] > 0.0) ? (int)symstart[axl[lvl_a]] : (nglob[axl[lvl_a]] - 1);
        const int nmax_b   = (symstart[axl[lvl_b]] > 0.0) ? (int)symstart[axl[lvl_b]] : (nglob[axl[lvl_b]] - 1);

        // the unbounded direction handled by the closed form (x) and by the quadrature (y). The
        // number of closed-form evaluations is proportional to the extent of the former, so the
        // narrower direction goes to the closed form. The choice is made on the *global* extents
        // (and so is the size of the quadrature rule below): the Green's function then stays a
        // function of the grid only, not of the decomposition.
        const bool x_is_a = (nmax_a <= nmax_b);
        const int  ax_x   = x_is_a ? axl[lvl_a] : axl[lvl_b];
        const int  ax_y   = x_is_a ? axl[lvl_b] : axl[lvl_a];
        const int  nx_max = x_is_a ? nmax_a : nmax_b;
        const int  ny_max = x_is_a ? nmax_b : nmax_a;

        // the interval of |n| this rank actually reads in each of them: the two tables are built
        // on that interval only, which divides their size -- and the number of closed-form
        // evaluations -- by the number of ranks in the corresponding direction
        int nx_lo = 0, nnx = 0, ny_lo = 0, nny = 0;
        cmpt_local_absID_(istart[ax_x], nloc[ax_x], symstart[ax_x], &nx_lo, &nnx);
        cmpt_local_absID_(istart[ax_y], nloc[ax_y], symstart[ax_y], &ny_lo, &nny);
        FLUPS_CHECK(nx_lo + nnx - 1 <= nx_max && ny_lo + nny - 1 <= ny_max,
                    "the local tables (%d:%d, %d:%d) exceed the global half-widths (%d, %d)",
                    nx_lo, nx_lo + nnx - 1, ny_lo, ny_lo + nny - 1, nx_max, ny_max);

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
        // the rule is sized per slice below; `kzmin` gives the largest of them, i.e. the size of
        // the buffers. `quad_size` is non-increasing in its wavenumber argument, so no slice can
        // ask for more than that.
        const int nq_max = is_hej_slice ? HEJTwoUnbounded::quad_size(ny_max, h_y, hej_kmax, length, kzmin)
                                        : LGFTwoUnbounded::quad_size(ny_max, kzmin_h);

        // the slices are not assembled one by one but by blocks of at most `nblk` contiguous
        // planes sharing the same rule: a per-slice product is far too skinny for the vendor BLAS
        // (see the stage 2 & 3 comment below). The block is capped both by a fixed number of
        // slices and by a memory budget on the two tall tables R and T.
        const int    nb_max       = 64;
        const size_t bytes_budget = (size_t)256 << 20;
        int          nblk         = (int)std::min<size_t>((size_t)nb_max, std::max<size_t>((size_t)1, bytes_budget / ((size_t)nnx * ((size_t)nq_max + (size_t)nny) * sizeof(double))));
        nblk                      = std::max(1, std::min(nblk, nlvl[lvl_s]));

        double *kq      = (double *)backend_malloc(sizeof(double) * nq_max);
        double *wq      = (double *)backend_malloc(sizeof(double) * nq_max);
        double *cw      = (double *)backend_malloc(sizeof(double) * nny * nq_max);
        // the tall tables of one block: R is column-major ((nblk * nnx) x Q), slice `b` occupying
        // the rows [b * nnx, (b+1) * nnx), and T = alpha * R Cw^T is column-major
        // ((nblk * nnx) x nny). Both keep the leading dimension of a *full* block, whatever the
        // number of slices the block ends up holding.
        double *rr = (double *)backend_malloc(sizeof(double) * (size_t)nblk * nnx * nq_max);
        double *tt = (double *)backend_malloc(sizeof(double) * (size_t)nblk * nnx * nny);
        // the lattice closed forms are split in a per-node part (the roots and the prefactors,
        // which only depend on the two wavenumbers) and a per-index part; the former lives in
        // this table, one entry per quadrature node, see @ref LGFTwoUnbounded::build_params
        void *qp = is_hej_slice ? NULL : backend_malloc(LGFTwoUnbounded::params_bytes() * nq_max);

        //---------------------------------------------------------------------
        /** - stage 2 & 3: rebuild R for every locally owned k != 0 plane, assemble the planes by
         *    blocks as one tall dgemm each and scatter them
         *
         * With R stored column-major as (nnx x Q) and Cw column-major as (nny x Q), one slice is
         * the single rank-Q product @f$ T = \alpha\, R\, Cw^T @f$ (`alpha = 1/pi`, with the sign
         * convention of the kernel family folded in), i.e. one `backend_dgemm_nt`, and the kernel
         * that follows is a pure scatter of T through the folded indexes.
         *
         * Such a per-slice product is however far too skinny to be of any use to the vendor BLAS:
         * with a few tens of rows and a few hundreds of columns it selects poor kernels and runs
         * at a couple of percent of peak. The slices are therefore accumulated into *blocks* of
         * contiguous planes sharing the same quadrature rule (hence the same Cw): slice `b` of
         * the block writes its R into the rows `[b * nnx, (b+1) * nnx)` of a tall table of leading
         * dimension `nblk * nnx`, and the whole block is assembled by a *single* product of
         * `blk_cnt * nnx` rows and scattered by a *single* kernel. The block size is the smallest
         * of a 64-slice cap, a memory budget on the two tall tables, and the number of local
         * planes. A block is flushed as soon as a slice cannot join it -- a rule change, a gap in
         * the plane indexes (the k = 0 plane), or a full block.
         *
         * The local index `isl` already visits the slices by increasing @f$ |k| @f$: the spectral
         * direction of the Green's function is always *unfolded* here (see the check below), so
         * @f$ k = (i_e + k_{off})\, k_{fact} @f$ grows with `isl`. That ordering is what makes the
         * per-slice sizing of the quadrature rule cheap: the endpoint near-singularity that sets
         * the floor of @ref LGFTwoUnbounded::quad_size has the width of the slice's own
         * wavenumber, so only the lowest modes need the largest rule; the rule size is then
         * non-increasing along the loop and the rule -- and the Cw table with it -- is rebuilt
         * only when it actually drops.
         */
        //---------------------------------------------------------------------
        // The plans are sorted by type (@ref FFTW_plan_dim::operator< ) and the spectral ones come
        // first (PERPER = 6 and SYMSYM = 2, both < UNBUNB = 8): the spectral direction of the
        // Green's function is therefore executed on real data, i.e. by an r2c/r2r plan, whose
        // output is never folded. Every local slice then carries its own, distinct, |k|.
        FLUPS_CHECK(symsta_s == 0.0,
                    "the spectral direction of the Green's function is expected unfolded (r2c/r2r plan), got symstart = %e", symsta_s);

        int nq       = 0;   // the size of the rule currently held by kq, wq and cw
        int blk_isl0 = -1;  // the first plane index of the block currently open
        int blk_cnt  = 0;   // the number of slices it holds
        // the lattice tables are in the positive-operator convention and FLUPS stores -G;
        // the regularised ones already carry the stored sign
        const double alpha = is_hej_slice ? c_1opi : -c_1opi;
        // the leading dimension of the two tall tables, that of a *full* block
        const int ldt = nblk * nnx;
        // one extra iteration, with no slice, to flush the last block
        for (int isl = 0; isl <= nlvl[lvl_s]; ++isl) {
            // the global index of that plane in the spectral direction, and its wavenumber
            bool   active = false;
            double kabs   = 0.0;
            int    nq_j   = 0;
            if (isl < nlvl[lvl_s]) {
                const int    ie = istart[ax_s] + isl;
                const double k  = (ie + koff_s) * kfact_s;
                kabs            = fabs(k);
                // the k = 0 plane has already been filled with the 2D table
                active = (kabs > kfact_s * 0.2);
                if (active) {
                    // the rule this slice needs, and the tables that go with it
                    nq_j = is_hej_slice ? HEJTwoUnbounded::quad_size(ny_max, h_y, hej_kmax, length, kabs)
                                        : LGFTwoUnbounded::quad_size(ny_max, kabs * h_spec);
                    FLUPS_CHECK(nq_j <= nq_max, "the rule of the slice |k| = %e (%d nodes) exceeds the allocated one (%d nodes)", kabs, nq_j, nq_max);
                }
            }

            // a slice joins the open block only if it shares its rule, follows it without a gap
            // and the block is not full yet
            const bool joins = active && (blk_cnt > 0) && (nq_j == nq) && (isl == blk_isl0 + blk_cnt) && (blk_cnt < nblk);

            //----------------------------------------------------------------- flush the block
            if (blk_cnt > 0 && !joins) {
                const int m = blk_cnt * nnx;
                backend_dgemm_nt(m, nny, nq, alpha, rr, ldt, cw, nny, tt, ldt, stream);

                const int isl0_    = blk_isl0;
                const int cnt_     = blk_cnt;
                auto      op_slice = m_backend_lambda(const int idx) {
                    const int per = na * nb;
                    const int bb  = idx / per;
                    const int rem = idx % per;

                    int ii[3];
                    ii[lvl_s]    = isl0_ + bb;
                    ii[lvl_a]    = rem % na;
                    ii[lvl_b]    = rem / na;
                    const int i0 = ii[0];
                    const int i1 = ii[1];
                    const int i2 = ii[2];

                    const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

                    int is[3];
                    cmpt_symID(ax0, i0, i1, i2, istart, symstart_, 0, is);
                    // the tables only span the local interval of |n|
                    const int nx = abs(is[ax_x]) - nx_lo;
                    const int ny = abs(is[ax_y]) - ny_lo;
                    assert(0 <= nx && nx < nnx && 0 <= ny && ny < nny);

                    // the block has been assembled by the dgemm above, this is a pure scatter
                    green[id + i0 * nf] = tt[(bb * nnx + nx) + ldt * ny];
                };
                parallel_for(op_slice, cnt_ * na * nb, "2unb1spe green slice", stream);

                blk_isl0 = -1;
                blk_cnt  = 0;
            }

            // the trailing iteration and the k = 0 plane stop here, the block is flushed
            if (!active) {
                continue;
            }

            //----------------------------------------------------------------- the rule
            // the block is empty here, so rebuilding Cw cannot invalidate a pending product
            if (nq_j != nq) {
                if (is_hej_slice) {
                    HEJTwoUnbounded::gauss_legendre_0kmax(nq_j, hej_kmax, kq, wq, stream);
                } else {
                    LGFTwoUnbounded::gauss_legendre_0pi(nq_j, 1.0, kq, wq, stream);
                }
                if (is_hej_slice) {
                    HEJTwoUnbounded::build_cw(cw, kq, wq, ny_lo, nny, nq_j, h_y, stream);
                } else {
                    LGFTwoUnbounded::build_cw(cw, kq, wq, ny_lo, nny, nq_j, stream);
                }
                nq = nq_j;
            }

            //----------------------------------------------------------------- R of that slice
            if (blk_cnt == 0) {
                blk_isl0 = isl;
            }
            const int b = blk_cnt;
            // the closed form is even in every wavenumber, so only |k| matters
            if (is_hej_slice) {
                HEJTwoUnbounded::build_r(rr + (size_t)b * nnx, kq, (int)typeGreen, nx_lo, nnx, ldt, nq, kabs, h_x, length, stream);
            } else {
                LGFTwoUnbounded::build_params(qp, kq, (int)typeGreen, nq, kabs * h_spec, stream);
                LGFTwoUnbounded::build_r(rr + (size_t)b * nnx, qp, (int)typeGreen, nx_lo, nnx, ldt, nq, stream);
            }
            blk_cnt += 1;
        }

        backend_sync_stream(stream);
        m_profStopi(prof, "green_kernel");

        backend_free(kq);
        backend_free(wq);
        backend_free(cw);
        if (qp != NULL) {
            backend_free(qp);
        }
        backend_free(rr);
        backend_free(tt);
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
            G  = kernel_from_symbol(&hej_2_1unb2spe);
            G0 = kernel_from_symbol(&hej_2_1unb2spe_k0);
            break;
        case HEJ_4:
            G  = kernel_from_symbol(&hej_4_1unb2spe);
            G0 = kernel_from_symbol(&hej_4_1unb2spe_k0);
            break;
        case HEJ_6:
            G  = kernel_from_symbol(&hej_6_1unb2spe);
            G0 = kernel_from_symbol(&hej_6_1unb2spe_k0);
            break;
        case HEJ_8:
            G  = kernel_from_symbol(&hej_8_1unb2spe);
            G0 = kernel_from_symbol(&hej_8_1unb2spe_k0);
            break;
        case HEJ_10:
            G  = kernel_from_symbol(&hej_10_1unb2spe);
            G0 = kernel_from_symbol(&hej_10_1unb2spe_k0);
            break;
        case HEJ_0:
            FLUPS_CHECK(false, "HEJ0 kernel not available for 1D unbounded problems.");
            break;
        case CHAT_2:
            G  = kernel_from_symbol(&chat_2_1unb2spe);
            G0 = kernel_from_symbol(&chat_2_1unb2spe_k0);
            break;
        case LGF_2:
            G  = kernel_from_symbol(&lgf_2_1unb2spe);
            G0 = kernel_from_symbol(&lgf_2_1unb2spe);
            break;
        case LGF_4:
            G  = kernel_from_symbol(&lgf_4_1unb2spe);
            G0 = kernel_from_symbol(&lgf_4_1unb2spe);
            break;
        case LGF_6:
            G  = kernel_from_symbol(&lgf_6_1unb2spe);
            G0 = kernel_from_symbol(&lgf_6_1unb2spe);
            break;
        case LGF_8:
            G  = kernel_from_symbol(&lgf_8_1unb2spe);
            G0 = kernel_from_symbol(&lgf_8_1unb2spe);
            break;
        case MEHR_4L:
            G  = kernel_from_symbol(&mehr_4l_1unb2spe);
            G0 = kernel_from_symbol(&mehr_4l_1unb2spe);
            break;
        case MEHR_6L:
            G  = kernel_from_symbol(&mehr_6l_1unb2spe);
            G0 = kernel_from_symbol(&mehr_6l_1unb2spe);
            break;
        case MEHR_4F:
            G  = kernel_from_symbol(&mehr_4f_1unb2spe);
            G0 = kernel_from_symbol(&mehr_4f_1unb2spe);
            break;
        case MEHR_6F:
            G  = kernel_from_symbol(&mehr_6f_1unb2spe);
            G0 = kernel_from_symbol(&mehr_6f_1unb2spe);
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
            G = kernel_from_symbol(&hej_2_0unb3spe);
            break;
        case HEJ_4:
            G = kernel_from_symbol(&hej_4_0unb3spe);
            break;
        case HEJ_6:
            G = kernel_from_symbol(&hej_6_0unb3spe);
            break;
        case HEJ_8:
            G = kernel_from_symbol(&hej_8_0unb3spe);
            break;
        case HEJ_10:
            G = kernel_from_symbol(&hej_10_0unb3spe);
            break;
        case HEJ_0:
            // spectral solution is here given by 1/k^2, i.e. same as the CHAT_2 kernel
            G = kernel_from_symbol(&chat_2_0unb3spe);
            break;
        case CHAT_2:
            G = kernel_from_symbol(&chat_2_0unb3spe);
            break;
        case LGF_2:
            G = kernel_from_symbol(&lgf_2_0unb3spe);
            break;
        case LGF_4:
            G = kernel_from_symbol(&lgf_4_0unb3spe);
            break;
        case LGF_6:
            G = kernel_from_symbol(&lgf_6_0unb3spe);
            break;
        case LGF_8:
            G = kernel_from_symbol(&lgf_8_0unb3spe);
            break;
        case MEHR_4L:
            G = kernel_from_symbol(&mehr_4l_0unb3spe);
            break;
        case MEHR_6L:
            G = kernel_from_symbol(&mehr_6l_0unb3spe);
            break;
        case MEHR_4F:
            G = kernel_from_symbol(&mehr_4f_0unb3spe);
            break;
        case MEHR_6F:
            G = kernel_from_symbol(&mehr_6f_0unb3spe);
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

/**
 * @brief returns the stencil that a lattice Green's function inverts exactly
 *
 * The weights are those of the symbols used by the kernels (see green_kernels.hpp and lgf_1unb2spe.hpp):
 * with \f$ \delta^2_{d,m} u = u(x + m h e_d) - 2 u(x) + u(x - m h e_d) \f$ (and \f$ \delta^2_d = \delta^2_{d,1} \f$),
 * - LGF_p: \f$ h^2 L_h = \sum_d \sum_m c_m \delta^2_{d,m} \f$ and \f$ M_h = 1 \f$;
 * - MEHR_4L/4F: \f$ h^2 L_h = \sum_d \delta^2_d + 1/6 \sum_{i<j} \delta^2_i \delta^2_j \f$ and
 *   \f$ M_h = 1 + 1/12 \sum_d \delta^2_d \f$;
 * - MEHR_6L/6F: \f$ h^2 L_h \f$ of MEHR_4L \f$ + 1/30\, \delta^2_0 \delta^2_1 \delta^2_2 \f$ and
 *   \f$ M_h = 1 + 1/12 \sum_d \delta^2_d + 1/90 \sum_{i<j} \delta^2_i \delta^2_j - 1/240 \sum_d \delta^4_d \f$.
 *
 * @param typeGreen the Green's function
 * @param lhs if true, the left-hand side stencil \f$ L_h \f$, else the right-hand side one \f$ M_h \f$
 * @param offsets if not NULL, the offsets of the n points, offsets[3*i + d] (size 3*n)
 * @param weights if not NULL, the weights of the n points, the left-hand side ones are in units of 1/h^2 (size n)
 * @return int the number of points n, 0 if the Green's function is not the inverse of a stencil
 */
int cmpt_Green_stencil(const GreenType typeGreen, const bool lhs, int *offsets, double *weights) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    // dense cube of weights, the widest stencil (LGF_8) reaches 4 points
    constexpr int r               = 4;
    constexpr int w               = 2 * r + 1;
    double        cube[w * w * w] = {0.0};

    const double c_d2[3] = {1.0, -2.0, 1.0};
    // add c times the product of unit second differences along the directions dirs[0 .. ndir-1]
    // (a repeated direction gives the fourth difference)
    auto add_d2_product = [&](const int ndir, const int dirs[], const double c) {
        int ncomb = 1;
        for (int id = 0; id < ndir; ++id) ncomb *= 3;
        for (int ic = 0; ic < ncomb; ++ic) {
            int    o[3] = {0, 0, 0};
            double wgt  = c;
            for (int id = 0, rest = ic; id < ndir; ++id, rest /= 3) {
                o[dirs[id]] += rest % 3 - 1;
                wgt *= c_d2[rest % 3];
            }
            cube[(o[0] + r) + w * ((o[1] + r) + w * (o[2] + r))] += wgt;
        }
    };
    // add c times the second difference of step m along d
    auto add_d2_step = [&](const int d, const int m, const double c) {
        for (int s = -1; s <= 1; ++s) {
            int o[3] = {0, 0, 0};
            o[d]     = s * m;
            cube[(o[0] + r) + w * ((o[1] + r) + w * (o[2] + r))] += c * c_d2[s + 1];
        }
    };
    auto add_pairs = [&](const double c) {
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                const int dirs[2] = {i, j};
                add_d2_product(2, dirs, c);
            }
        }
    };

    const double  c_lgf2[1] = {1.0};
    const double  c_lgf4[2] = {4.0 / 3.0, -1.0 / 12.0};
    const double  c_lgf6[3] = {3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0};
    const double  c_lgf8[4] = {8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0};
    const double *c_lgf     = NULL;
    int           n_lgf     = 0;
    switch (typeGreen) {
        case LGF_2:
            c_lgf = c_lgf2;
            n_lgf = 1;
            break;
        case LGF_4:
            c_lgf = c_lgf4;
            n_lgf = 2;
            break;
        case LGF_6:
            c_lgf = c_lgf6;
            n_lgf = 3;
            break;
        case LGF_8:
            c_lgf = c_lgf8;
            n_lgf = 4;
            break;
        default:
            break;
    }
    const bool is_mehr4 = (typeGreen == MEHR_4L || typeGreen == MEHR_4F);
    const bool is_mehr6 = (typeGreen == MEHR_6L || typeGreen == MEHR_6F);

    if (c_lgf != NULL) {
        if (lhs) {
            for (int d = 0; d < 3; ++d) {
                for (int m = 0; m < n_lgf; ++m) add_d2_step(d, m + 1, c_lgf[m]);
            }
        } else {
            cube[r + w * (r + w * r)] = 1.0;
        }
    } else if (is_mehr4 || is_mehr6) {
        if (lhs) {
            for (int d = 0; d < 3; ++d) add_d2_step(d, 1, 1.0);
            add_pairs(c_1o6);
            if (is_mehr6) {
                const int dirs[3] = {0, 1, 2};
                add_d2_product(3, dirs, c_1o30);
            }
        } else {
            cube[r + w * (r + w * r)] = 1.0;
            for (int d = 0; d < 3; ++d) add_d2_step(d, 1, c_1o12);
            if (is_mehr6) {
                add_pairs(1.0 / 90.0);
                for (int d = 0; d < 3; ++d) {
                    const int dirs[2] = {d, d};
                    add_d2_product(2, dirs, -1.0 / 240.0);
                }
            }
        }
    }

    int n = 0;
    for (int i2 = -r; i2 <= r; ++i2) {
        for (int i1 = -r; i1 <= r; ++i1) {
            for (int i0 = -r; i0 <= r; ++i0) {
                const double wgt = cube[(i0 + r) + w * ((i1 + r) + w * (i2 + r))];
                if (wgt == 0.0) continue;
                if (offsets != NULL) {
                    offsets[3 * n + 0] = i0;
                    offsets[3 * n + 1] = i1;
                    offsets[3 * n + 2] = i2;
                }
                if (weights != NULL) weights[n] = wgt;
                ++n;
            }
        }
    }
    //-------------------------------------------------------------------------
    END_FUNC;
    return n;
}
