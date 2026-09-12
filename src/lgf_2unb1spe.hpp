/**
 * @file lgf_2unb1spe.hpp
 * @brief exact screened slice Lattice Green's Functions for 2 unbounded + 1 spectral directions
 *
 * For an operator @f$ L @f$ that is periodic (or symmetric) in @f$ z @f$ and whose symbol is a
 * polynomial in @f$ s(k_z) @f$, the discrete Fourier transform along @f$ z @f$ alone
 * block-diagonalises the problem: the mode @f$ k_z @f$ obeys a two-dimensional lattice equation
 * with symbol @f$ \sigma_{k_z}(k_x, k_y) = \sigma(k_x, k_y, k_z) @f$. For @f$ k_z \ne 0 @f$ that
 * two-dimensional operator is a *screened* (modified Helmholtz) one, its Green's function decays
 * exponentially, and it can be evaluated to machine precision as
 *
 * @f[ G_{k_z}[n_x, n_y] = \frac{1}{\pi} \int_0^{\pi} \cos(n_y k)\; I(n_x; k, k_z)\; dk @f]
 *
 * where @f$ I(n; k_y, k_z) = \frac{1}{2\pi}\int_{-\pi}^{\pi} \frac{e^{i n k_x}}{\sigma} dk_x @f$
 * is the one-unbounded/two-spectral LGF already implemented in @ref lgf_1unb2spe.hpp. For
 * @f$ k_z \ne 0 @f$ the outer integrand is analytic on the *closed* interval @f$[0,\pi]@f$
 * (nothing vanishes, since @f$ \sigma_{k_z} @f$ is bounded away from zero) so a Gauss-Legendre
 * rule converges exponentially; @f$ Q \simeq 2 N @f$ nodes give machine precision.
 *
 * The evaluation is cached in two stages (both tables are tiny):
 *   -# @f$ Cw[q, n_y] = w_q \cos(n_y k_q) @f$, built once for the whole Green's function,
 *   -# @f$ R[q, n_x] = I(n_x; k_q, |k_z| h) @f$, rebuilt for each spectral slice,
 * and each point of a slice is then a length-@f$Q@f$ dot product, see @ref slice_point.
 *
 * Both tables only span the interval of @f$ |n| @f$ the calling rank actually reads: they are
 * built from a first index @f$ n_{lo} @f$ and a count, and the entry @f$ n @f$ lives at
 * @f$ n - n_{lo} @f$. The quadrature rule itself is a global property of the grid (it is sized
 * from the *global* half-width of the quadrature direction), so every rank gets the same values.
 *
 * The inner closed form is itself split in two (see @ref lgf_1unb2spe.hpp): its roots, prefactors
 * and branch selection depend on @f$ (k_q, k_z) @f$ only and are evaluated once per quadrature
 * node into a third, tiny, table (@ref build_params), leaving @f$ r^n @f$ as the only per-entry
 * work of @ref build_r. That table also carries the index beyond which the entry is below
 * @f$ 10^{-17} @f$ and is simply written as zero.
 *
 * All wavenumbers handled here are **dimensionless** (i.e. @f$ k h @f$) and all values are given
 * in the POSITIVE-operator convention (@f$ G > 0 @f$, decaying). FLUPS stores @f$ -G @f$, exactly
 * like the @f$ k_z = 0 @f$ two-dimensional tables of `lgf_*_2unb0spe_`.
 *
 * @copyright Copyright © UCLouvain 2020
 *
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 */
#ifndef LGF_TWO_UNBOUNDED_HPP_
#define LGF_TWO_UNBOUNDED_HPP_

#include <algorithm>
#include <cmath>

#include "backend.hpp"
#include "defines.hpp"
#include "lgf_1unb2spe.hpp"

namespace LGFTwoUnbounded {

//==============================================================================
// which kernels are supported
//==============================================================================
/**
 * @brief true if the exact screened-slice construction is available for this kernel
 *
 * The construction only needs the one-dimensional closed form @ref slice_inner, which exists for
 * every LGF and Mehrstellen kernel of FLUPS.
 */
inline bool handles(const GreenType typeGreen) {
    return (typeGreen == LGF_2 || typeGreen == LGF_4 || typeGreen == LGF_6 || typeGreen == LGF_8 ||
            typeGreen == MEHR_4L || typeGreen == MEHR_6L || typeGreen == MEHR_4F || typeGreen == MEHR_6F);
}

//==============================================================================
// the inner (closed-form) integral
//==============================================================================
/**
 * @brief the one-dimensional screened LGF @f$ I(n; k_y, k_z) @f$, POSITIVE operator convention
 *
 * This is exactly the kernel used by the `*_1unb2spe_` Green's functions, with the two spectral
 * wavenumbers being here the Gauss-Legendre node @p ky and the (dimensionless, positive) spectral
 * wavenumber @p kz. It is the reference form of what @ref build_params and @ref build_r do in two
 * steps; the tables are not built through it.
 *
 * @param typeGreen the kernel, must satisfy @ref handles
 * @param n the (non-negative) index along the closed-form direction
 * @param ky the dimensionless wavenumber of the quadrature direction, in [0, pi]
 * @param kz the dimensionless wavenumber of the spectral direction, >= 0
 */
m_device double slice_inner(const int typeGreen, const int n, const double ky, const double kz) {
    switch (typeGreen) {
        case LGF_2:
            return LGFOneUnbounded::lgf2(n, LGFOneUnbounded::lgf2_symbol(ky) + LGFOneUnbounded::lgf2_symbol(kz));
        case LGF_4:
            return LGFOneUnbounded::lgf4(n, LGFOneUnbounded::lgf4_symbol(ky) + LGFOneUnbounded::lgf4_symbol(kz));
        case LGF_6:
            return LGFOneUnbounded::lgf6(n, LGFOneUnbounded::lgf6_symbol(ky) + LGFOneUnbounded::lgf6_symbol(kz));
        case LGF_8:
            return LGFOneUnbounded::lgf8(n, LGFOneUnbounded::lgf8_symbol(ky) + LGFOneUnbounded::lgf8_symbol(kz));
        // for the Mehrstellen kernels the two transverse wavenumbers are ordered |k1| >= |k2|,
        // exactly as the mehr_*_1unb2spe_ kernels do (the expressions are symmetric in the two)
        case MEHR_4L:
            return LGFOneUnbounded::meh4_left(n, fmax(ky, kz), fmin(ky, kz));
        case MEHR_6L:
            return LGFOneUnbounded::meh6_left(n, fmax(ky, kz), fmin(ky, kz));
        case MEHR_4F:
            return LGFOneUnbounded::meh4_full(n, fmax(ky, kz), fmin(ky, kz));
        case MEHR_6F:
            return LGFOneUnbounded::meh6_full(n, fmax(ky, kz), fmin(ky, kz));
        default:
            assert(false && "[LGFTwoUnbounded] unsupported Green kernel");
            return 0.0;
    }
}

/**
 * @brief assemble one point of a slice from the two cached tables, POSITIVE operator convention
 *
 * @param R the (Q x nnx) table, `R[q * nnx + nx]`
 * @param Cw the (Q x nny) table, `Cw[q * nny + ny]`
 * @param nx the table index along the closed-form direction, i.e. `|n_x| - nx_lo`, in [0, nnx)
 * @param ny the table index along the quadrature direction, i.e. `|n_y| - ny_lo`, in [0, nny)
 * @param nnx the leading dimension of R, i.e. the local count of the closed-form direction
 * @param nny the leading dimension of Cw, i.e. the local count of the quadrature direction
 * @param nq the number of quadrature nodes
 */
m_device double slice_point(const double* R, const double* Cw, const int nx, const int ny,
                            const int nnx, const int nny, const int nq) {
    double sum = 0.0;
    for (int q = 0; q < nq; ++q) {
        sum += R[q * nnx + nx] * Cw[q * nny + ny];
    }
    return c_1opi * sum;
}

//==============================================================================
// the Gauss-Legendre rule on [0, pi] (host)
//==============================================================================
/**
 * @brief the number of equal panels used for a composite rule of @p nq nodes
 *
 * The integrand is analytic on [0, pi] so panels are only there to keep the per-panel order (and
 * hence the Newton solve for the nodes) moderate.
 */
inline int quad_npanel(const int nq) { return (nq + 511) / 512; }

/**
 * @brief the number of Gauss-Legendre nodes used for a table of half-width @p n
 *
 * The baseline is `Q = 2 n + 32`, which resolves the `cos(n_y k)` oscillation of the outer
 * integrand. On top of that, the inner integral `I(n_x; k, k_z)` has a near-singularity of width
 * `~ k_z h` sitting at the endpoint `k = 0` (the operator becomes singular at `k = k_z = 0`), so
 * the smallest spectral wavenumber of the grid also sets a floor: measured on LGF_2 with n = 16,
 * `Q ~ 2 pi / (k_z h)` (i.e. `Q ~ N_z`) is needed to stay at 1e-14, and `Q = 64` already gives
 * 2.8e-11 at `k_z h = 2 pi / 128`. Pass @p kzmin_h to take that into account (it only bites when
 * the spectral direction is much longer than the unbounded ones, which is exactly where the old
 * sampled-symbol treatment was worst).
 *
 * The result is rounded up so that it splits evenly over the panels.
 *
 * @param n the half-width of the table in the quadrature direction
 * @param kzmin_h the smallest non-zero dimensionless spectral wavenumber of the grid, 0 to ignore
 */
inline int quad_size(const int n, const double kzmin_h = 0.0) {
    int q0 = 2 * n + 32;
    if (kzmin_h > 0.0) {
        q0 = std::max(q0, (int)ceil(2.0 * M_PI / kzmin_h) + 32);
    }
    const int np = quad_npanel(q0);
    return ((q0 + np - 1) / np) * np;
}

/**
 * @brief the Gauss-Legendre nodes and weights of order @p n on [-1, 1] (Newton on P_n)
 */
inline void gauss_legendre_ref_(const int n, double* x, double* w) {
    const int nhalf = (n + 1) / 2;
    for (int i = 0; i < nhalf; ++i) {
        // Chebyshev-like initial guess, accurate to O(1/n^2)
        double z  = cos(M_PI * (i + 0.75) / (n + 0.5));
        double pp = 0.0;
        for (int it = 0; it < 200; ++it) {
            // Legendre recurrence: p1 = P_n(z), p2 = P_{n-1}(z)
            double p1 = 1.0, p2 = 0.0;
            for (int j = 0; j < n; ++j) {
                const double p3 = p2;
                p2              = p1;
                p1              = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            pp              = n * (z * p1 - p2) / (z * z - 1.0);
            const double dz = p1 / pp;
            z -= dz;
            if (fabs(dz) <= 2.0e-16 * (fabs(z) + 1.0)) {
                break;
            }
        }
        x[i]         = -z;
        x[n - 1 - i] = z;
        w[i]         = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

/**
 * @brief the composite Gauss-Legendre rule of @p nq nodes on [0, pi]
 *
 * @param nq the total number of nodes, must be a multiple of @ref quad_npanel(nq)
 * @param[out] k the nodes, of size nq
 * @param[out] w the weights, of size nq (they sum to pi)
 */
inline void gauss_legendre_0pi(const int nq, double* k, double* w) {
    const int np = quad_npanel(nq);
    FLUPS_CHECK(nq % np == 0, "the number of quadrature nodes (%d) must split evenly over %d panels", nq, np);
    const int nper = nq / np;

    double* xr = (double*)m_calloc(sizeof(double) * nper);
    double* wr = (double*)m_calloc(sizeof(double) * nper);
    gauss_legendre_ref_(nper, xr, wr);

    for (int p = 0; p < np; ++p) {
        const double a = M_PI * ((double)p) / ((double)np);
        const double b = M_PI * ((double)(p + 1)) / ((double)np);
        for (int i = 0; i < nper; ++i) {
            k[p * nper + i] = 0.5 * (xr[i] + 1.0) * (b - a) + a;
            w[p * nper + i] = 0.5 * wr[i] * (b - a);
        }
    }
    m_free(xr);
    m_free(wr);
}

//==============================================================================
// the cached tables (device buffers, filled by parallel_for)
//==============================================================================
/**
 * @brief fill `Cw[q * nny + j] = w_q cos((ny_lo + j) k_q)` on the device
 *
 * @param ny_lo the first (non-negative) index of the quadrature direction held by this rank
 * @param nny the number of indexes held by this rank, i.e. the leading dimension of Cw
 */
inline void build_cw(double* cw, const double* kq, const double* wq, const int ny_lo, const int nny,
                     const int nq, stream_t stream) {
    auto op = m_backend_lambda(const int idx) {
        const int j = idx % nny;
        const int q = idx / nny;
        cw[idx]     = wq[q] * cos(((double)(ny_lo + j)) * kq[q]);
    };
    parallel_for(op, nny * nq, "2unb1spe green Cw", stream);
}

//==============================================================================
// the closed form, split in a per-wavenumber part and a per-index part
//==============================================================================
/**
 * @brief the params/eval pair of one kernel, see the split in @ref lgf_1unb2spe.hpp
 *
 * The roots, the prefactors and the branch selection of a closed form depend on the two
 * wavenumbers only: they are evaluated once per quadrature node (@ref build_params) and stored in
 * a small device table, so that the `(Nx + 1) x Q` entries of R cost nothing but @f$ r^n @f$.
 *
 * @tparam TYPE the kernel, must satisfy @ref handles
 */
template <int TYPE>
struct kernel_traits_;

template <>
struct kernel_traits_<LGF_2> {
    using params_t = LGFOneUnbounded::Lgf2Params;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::lgf2_params(LGFOneUnbounded::lgf2_symbol(ky) + LGFOneUnbounded::lgf2_symbol(kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::lgf2_eval(p, n); }
};

template <>
struct kernel_traits_<LGF_4> {
    using params_t = LGFOneUnbounded::Lgf4Params;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::lgf4_params(LGFOneUnbounded::lgf4_symbol(ky) + LGFOneUnbounded::lgf4_symbol(kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::lgf4_eval(p, n); }
};

template <>
struct kernel_traits_<LGF_6> {
    using params_t = LGFOneUnbounded::Lgf6Params;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::lgf6_params(LGFOneUnbounded::lgf6_symbol(ky) + LGFOneUnbounded::lgf6_symbol(kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::lgf6_eval(p, n); }
};

template <>
struct kernel_traits_<LGF_8> {
    using params_t = LGFOneUnbounded::Lgf8Params;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::lgf8_params(LGFOneUnbounded::lgf8_symbol(ky) + LGFOneUnbounded::lgf8_symbol(kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::lgf8_eval(p, n); }
};

// the Mehrstellen closed forms take the two transverse wavenumbers ordered |k1| >= |k2|, exactly
// as the mehr_*_1unb2spe_ kernels do (the expressions are symmetric in the two)
template <>
struct kernel_traits_<MEHR_4L> {
    using params_t = LGFOneUnbounded::MehParams;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::meh4_left_params(fmax(ky, kz), fmin(ky, kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::meh4_left_eval(p, n); }
};

template <>
struct kernel_traits_<MEHR_6L> {
    using params_t = LGFOneUnbounded::MehParams;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::meh6_left_params(fmax(ky, kz), fmin(ky, kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::meh6_left_eval(p, n); }
};

template <>
struct kernel_traits_<MEHR_4F> {
    using params_t = LGFOneUnbounded::MehParams;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::meh4_full_params(fmax(ky, kz), fmin(ky, kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::meh4_full_eval(p, n); }
};

template <>
struct kernel_traits_<MEHR_6F> {
    using params_t = LGFOneUnbounded::MehParams;
    static m_device params_t params(const double ky, const double kz) {
        return LGFOneUnbounded::meh6_full_params(fmax(ky, kz), fmin(ky, kz));
    }
    static m_device double eval(const params_t& p, const int n) { return LGFOneUnbounded::meh6_full_eval(p, n); }
};

/**
 * @brief the size of one entry of the params table, i.e. the largest of the kernels
 *
 * The table is allocated once, for the largest kernel, and cast to the right type by
 * @ref build_params_ and @ref build_r_.
 */
inline constexpr size_t params_bytes() {
    return std::max({sizeof(LGFOneUnbounded::Lgf2Params), sizeof(LGFOneUnbounded::Lgf4Params),
                     sizeof(LGFOneUnbounded::Lgf6Params), sizeof(LGFOneUnbounded::Lgf8Params),
                     sizeof(LGFOneUnbounded::MehParams)});
}

/**
 * @brief fill the params table `P[q]` of one slice on the device, for one kernel
 *
 * @tparam TYPE the kernel, must satisfy @ref handles
 * @param pbuf the params table, of `nq * params_bytes()` bytes
 * @param kzh the dimensionless (and non-negative) wavenumber of the spectral direction
 */
template <int TYPE>
inline void build_params_(void* pbuf, const double* kq, const int nq, const double kzh, stream_t stream) {
    using params_t  = typename kernel_traits_<TYPE>::params_t;
    params_t* const P = (params_t*)pbuf;
    auto op = m_backend_lambda(const int q) {
        P[q] = kernel_traits_<TYPE>::params(kq[q], kzh);
    };
    parallel_for(op, nq, "2unb1spe green params", stream);
}

/**
 * @brief fill `R[q * nnx + i] = I(nx_lo + i; k_q, kzh)` on the device, for one kernel
 *
 * The kernel is a template parameter so that the dispatch is resolved at compile time and a
 * single closed form is inlined in the device code, see @ref build_r. Everything that does not
 * depend on the index has already been evaluated into @p pbuf by @ref build_params.
 *
 * Entries beyond `P[q].n_cut` are written as a plain zero instead of being evaluated: the closed
 * form decays like @f$ \rho^n @f$ there and is below @ref LGFOneUnbounded::lgf_trunc_, i.e. under
 * the round-off the length-Q dot product of @ref slice_point commits anyway. The table therefore
 * carries an absolute floor of ~1e-17, deliberately.
 *
 * @tparam TYPE the kernel, must satisfy @ref handles
 * @param pbuf the params table filled by @ref build_params for this slice
 * @param nx_lo the first (non-negative) index of the closed-form direction held by this rank
 * @param nnx the number of indexes held by this rank, i.e. the leading dimension of R
 */
template <int TYPE>
inline void build_r_(double* r, const void* pbuf, const int nx_lo, const int nnx, const int nq,
                     stream_t stream) {
    using params_t          = typename kernel_traits_<TYPE>::params_t;
    const params_t* const P = (const params_t*)pbuf;
    auto op = m_backend_lambda(const int idx) {
        const int i = idx % nnx;
        const int q = idx / nnx;
        const int n = nx_lo + i;
        r[idx]      = (n > P[q].n_cut) ? 0.0 : kernel_traits_<TYPE>::eval(P[q], n);
    };
    parallel_for(op, nnx * nq, "2unb1spe green R", stream);
}

/**
 * @brief fill the params table `P[q]` of one slice on the device
 *
 * The dispatch on the kernel happens here, on the host, once per slice.
 *
 * @param typeGreen the kernel, must satisfy @ref handles
 * @param kzh the dimensionless (and non-negative) wavenumber of the spectral direction
 */
inline void build_params(void* pbuf, const double* kq, const int typeGreen, const int nq,
                         const double kzh, stream_t stream) {
    switch (typeGreen) {
        case LGF_2:
            build_params_<LGF_2>(pbuf, kq, nq, kzh, stream);
            break;
        case LGF_4:
            build_params_<LGF_4>(pbuf, kq, nq, kzh, stream);
            break;
        case LGF_6:
            build_params_<LGF_6>(pbuf, kq, nq, kzh, stream);
            break;
        case LGF_8:
            build_params_<LGF_8>(pbuf, kq, nq, kzh, stream);
            break;
        case MEHR_4L:
            build_params_<MEHR_4L>(pbuf, kq, nq, kzh, stream);
            break;
        case MEHR_6L:
            build_params_<MEHR_6L>(pbuf, kq, nq, kzh, stream);
            break;
        case MEHR_4F:
            build_params_<MEHR_4F>(pbuf, kq, nq, kzh, stream);
            break;
        case MEHR_6F:
            build_params_<MEHR_6F>(pbuf, kq, nq, kzh, stream);
            break;
        default:
            FLUPS_CHECK(false, "[LGFTwoUnbounded] unsupported Green kernel %d", typeGreen);
    }
}

/**
 * @brief fill `R[q * nnx + i] = I(nx_lo + i; k_q, kzh)` on the device
 *
 * The dispatch on the kernel happens here, on the host, once per slice: @ref build_r_ is then
 * instantiated for a single closed form.
 *
 * @param typeGreen the kernel, must satisfy @ref handles
 * @param pbuf the params table filled by @ref build_params for this slice
 */
inline void build_r(double* r, const void* pbuf, const int typeGreen, const int nx_lo, const int nnx,
                    const int nq, stream_t stream) {
    switch (typeGreen) {
        case LGF_2:
            build_r_<LGF_2>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case LGF_4:
            build_r_<LGF_4>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case LGF_6:
            build_r_<LGF_6>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case LGF_8:
            build_r_<LGF_8>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case MEHR_4L:
            build_r_<MEHR_4L>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case MEHR_6L:
            build_r_<MEHR_6L>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case MEHR_4F:
            build_r_<MEHR_4F>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        case MEHR_6F:
            build_r_<MEHR_6F>(r, pbuf, nx_lo, nnx, nq, stream);
            break;
        default:
            FLUPS_CHECK(false, "[LGFTwoUnbounded] unsupported Green kernel %d", typeGreen);
    }
}

};  // namespace LGFTwoUnbounded

#endif  // LGF_TWO_UNBOUNDED_HPP_
