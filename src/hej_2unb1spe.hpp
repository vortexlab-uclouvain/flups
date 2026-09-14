/**
 * @file hej_2unb1spe.hpp
 * @brief exact screened slice kernels for the regularised (HEJ) Green's functions,
 *        2 unbounded + 1 spectral directions
 *
 * The HEJ kernels mollify the *continuum* Poisson equation: the Green's function is the free-space
 * one convolved with a regularisation of width @f$ \sigma @f$. A discrete Fourier transform along
 * the spectral direction alone leaves, for each mode @f$ k_z @f$, a two-dimensional *screened*
 * problem whose regularised Green's function @f$ G_{k_z}(x,y) @f$ is what the Hockney-Eastwood
 * mirror + FFT needs in physical space.
 *
 * That two-dimensional kernel is one Fourier transform away from the one-unbounded/two-spectral
 * kernel FLUPS already ships (@ref hej_2_1unb2spe_ and friends):
 *
 * @f[ G_{k_z}(x, y) = \frac{1}{\pi} \int_0^{\infty}
 *        \cos(k_y\, y)\; G_1\!\left(x;\; \sqrt{k_y^2 + k_z^2},\; \sigma\right) dk_y @f]
 *
 * because @f$ G_1 @f$ depends on the two transverse wavenumbers only through their magnitude (the
 * mollifier is isotropic). Nothing new has to be derived: the integrand is an existing, validated
 * closed form, and the construction is the exact analogue of @ref lgf_2unb1spe.hpp with the
 * Brillouin zone @f$[0,\pi]@f$ replaced by the semi-infinite line.
 *
 * Two differences with the lattice case are worth noting:
 *   - the integrand decays like @f$ e^{-(k\sigma)^2/2} @f$ (and faster still away from the origin,
 *     like @f$ e^{-k r} @f$), so the range can be truncated at @f$ k_{max} = S/\sigma @f$ with
 *     @f$ S \simeq 10 @f$ and nothing measurable is lost, see @ref quad_kmax;
 *   - the kernel is a *continuum* one, so the two unbounded directions need not share a grid
 *     spacing: @f$ x = n_x h_x @f$ and @f$ y = n_y h_y @f$ enter separately.
 *
 * The evaluation is cached in the same two stages as the lattice construction:
 *   -# @f$ Cw[q, n_y] = w_q \cos(n_y h_y k_q) @f$, built once for the whole Green's function,
 *   -# @f$ R[q, n_x] = G_1(n_x h_x; \sqrt{k_q^2 + k_z^2}, \sigma) @f$, rebuilt for each slice,
 * and, exactly as in the lattice case, the whole slice is assembled at once as the rank-@f$Q@f$
 * product @f$ T = R\, Cw^T / \pi @f$ -- a single `dgemm` (@ref backend_dgemm_nt) on the
 * column-major @f$ (n_{nx} \times Q) @f$ and @f$ (n_{ny} \times Q) @f$ tables -- followed by a
 * scatter of @f$ T[n_x + n_{nx} n_y] @f$ through the folded indexes of the Green's function.
 *
 * As in the lattice case both tables only span the interval of @f$ |n| @f$ the calling rank
 * actually reads: they are built from a first index @f$ n_{lo} @f$ and a count, and the entry
 * @f$ n @f$ lives at @f$ n - n_{lo} @f$. The quadrature rule stays a global property of the grid.
 *
 * Unlike @ref lgf_2unb1spe.hpp the values here are already in the sign convention FLUPS stores
 * (@f$ G < 0 @f$), because the @f$ *\_1unb2spe\_ @f$ kernels are. No sign flip is applied by the
 * caller. The @f$ r = 0 @f$ point needs no special case: a regularised kernel is finite there.
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
#ifndef HEJ_TWO_UNBOUNDED_HPP_
#define HEJ_TWO_UNBOUNDED_HPP_

#include <algorithm>
#include <cmath>

#include "backend.hpp"
#include "defines.hpp"
#include "green_kernels.hpp"
#include "lgf_2unb1spe.hpp"

namespace HEJTwoUnbounded {

//==============================================================================
// which kernels are supported
//==============================================================================
/**
 * @brief true if the exact screened-slice construction is available for this kernel
 *
 * HEJ_0 is excluded: it is a Fourier cut-off (spectral-like) kernel with no 1unb2spe closed form
 * of this shape, and it does not represent an unbounded direction in the first place.
 */
inline bool handles(const GreenType typeGreen) {
    return (typeGreen == HEJ_2 || typeGreen == HEJ_4 || typeGreen == HEJ_6 ||
            typeGreen == HEJ_8 || typeGreen == HEJ_10);
}

//==============================================================================
// the inner (closed-form) integrand
//==============================================================================
/**
 * @brief the one-unbounded/two-spectral regularised kernel @f$ G_1(r; k, \sigma) @f$
 *
 * This forwards to the very same expressions used by the `hej_*_1unb2spe_` Green kernels, so the
 * construction inherits whatever they are worth and nothing is re-derived here.
 *
 * @param typeGreen the kernel, must satisfy @ref handles
 * @param r the (non-negative) position along the closed-form direction
 * @param k the magnitude of the transverse wavenumber, > 0
 * @param sig the regularisation length
 */
m_device double slice_inner(const int typeGreen, const double r, const double k, const double sig) {
    const double p[3] = {r, k, sig};
    switch (typeGreen) {
        case HEJ_2:
            return hej_2_1unb2spe_(p, NULL);
        case HEJ_4:
            return hej_4_1unb2spe_(p, NULL);
        case HEJ_6:
            return hej_6_1unb2spe_(p, NULL);
        case HEJ_8:
            return hej_8_1unb2spe_(p, NULL);
        case HEJ_10:
            return hej_10_1unb2spe_(p, NULL);
        default:
            assert(false && "[HEJTwoUnbounded] unsupported Green kernel");
            return 0.0;
    }
}

//==============================================================================
// the quadrature rule
//==============================================================================
/**
 * @brief the truncation of the semi-infinite range, in wavenumber units
 *
 * At @f$ r = 0 @f$ (the slowest-decaying point of a slice) the integrand behaves like
 * @f$ e^{-(k\sigma)^2/2} @f$, so cutting at @f$ k\sigma = S @f$ leaves a relative tail of order
 * @f$ e^{-S^2/2} @f$. @f$ S = 9 @f$ gives @f$ 2\cdot10^{-18} @f$, below what the double-precision
 * FFT that follows can carry. Away from the origin the decay is @f$ e^{-kr} @f$ on top of that.
 *
 * @param sig the regularisation length
 */
inline double quad_kmax(const double sig) {
    constexpr double S = 9.0;
    return S / sig;
}

/**
 * @brief the number of Gauss-Legendre nodes for a slice construction on [0, kmax]
 *
 * Two features have to be resolved, exactly as in the lattice case:
 *   - the outer oscillation @f$ \cos(n_y h_y k) @f$, whose fastest period over the range needs
 *     about @f$ 2 k_{max} n_{y,max} h_y / \pi @f$ nodes;
 *   - the near-singularity of width @f$ \sim k_z @f$ that @f$ G_1 \sim -1/(2k) @f$ develops at
 *     @f$ k \to 0 @f$; the smallest spectral wavenumber of the grid therefore sets a floor of
 *     about @f$ 2 k_{max} / k_{z,min} @f$ nodes.
 *
 * The result is rounded up so that it splits evenly over the panels.
 *
 * @param ny_max the half-width of the table in the quadrature direction
 * @param hy the grid spacing of the quadrature direction
 * @param kmax the truncation of the range, see @ref quad_kmax
 * @param sig the regularisation length
 * @param kzmin the smallest non-zero spectral wavenumber of the grid, 0 to ignore
 */
inline int quad_size(const int ny_max, const double hy, const double kmax, const double sig, const double kzmin) {
    // the outer oscillation: cos(n_y h_y k) runs over k_max n_y h_y / (2 pi) periods on the range
    const int q_osc = (int)ceil(2.0 * kmax * ((double)ny_max) * hy * c_1opi);
    // the integrand itself varies on the scale of the regularisation, i.e. over k_max sig units of
    // the dimensionless k sig; four nodes per unit is ample for an analytic integrand
    const int q_reg = (int)ceil(4.0 * kmax * sig);
    int       q0    = q_osc + q_reg + 64;
    if (kzmin > 0.0) {
        q0 = std::max(q0, (int)ceil(2.0 * kmax / kzmin) + q_reg + 64);
    }
    const int np = LGFTwoUnbounded::quad_npanel(q0);
    return ((q0 + np - 1) / np) * np;
}

/**
 * @brief the composite Gauss-Legendre rule of @p nq nodes on [0, kmax]
 *
 * @param nq the total number of nodes, must be a multiple of the panel count
 * @param kmax the upper end of the range
 * @param[out] k the nodes, of size nq
 * @param[out] w the weights, of size nq (they sum to kmax)
 */
inline void gauss_legendre_0kmax(const int nq, const double kmax, double* k, double* w) {
    // the [0, pi] rule of the lattice construction, rescaled to [0, kmax]
    LGFTwoUnbounded::gauss_legendre_0pi(nq, k, w);
    const double fact = kmax * c_1opi;
    for (int i = 0; i < nq; ++i) {
        k[i] *= fact;
        w[i] *= fact;
    }
}

//==============================================================================
// the cached tables (device buffers, filled by parallel_for)
//==============================================================================
/**
 * @brief the value of @f$ k r @f$ beyond which a table entry is not worth evaluating
 *
 * @f$ e^{-40} \simeq 4\cdot 10^{-18} @f$ relative to the @f$ r = 0 @f$ value of the same slice.
 * The `s * rho > 100` cut the `hej_*_1unb2spe_` kernels already carry stays, this one only skips
 * the call altogether.
 */
m_constant constexpr double far_cut_ = 40.0;

/**
 * @brief fill `Cw[q * nny + j] = w_q cos((ny_lo + j) hy k_q)` on the device
 *
 * The layout is that of a column-major @f$ (n_{ny} \times Q) @f$ matrix of leading dimension
 * @p nny, i.e. the second operand of the slice `dgemm` (@ref backend_dgemm_nt).
 *
 * @param ny_lo the first (non-negative) index of the quadrature direction held by this rank
 * @param nny the number of indexes held by this rank, i.e. the leading dimension of Cw
 */
inline void build_cw(double* cw, const double* kq, const double* wq, const int ny_lo, const int nny,
                     const int nq, const double hy, stream_t stream) {
    auto op = m_backend_lambda(const int idx) {
        const int j = idx % nny;
        const int q = idx / nny;
        cw[idx]     = wq[q] * cos(((double)(ny_lo + j)) * hy * kq[q]);
    };
    parallel_for(op, nny * nq, "2unb1spe hej green Cw", stream);
}

/**
 * @brief fill `R[q * nnx + i] = G1((nx_lo + i) hx; sqrt(k_q^2 + kz^2), sig)` on the device,
 *        for one kernel
 *
 * The kernel is a template parameter so that the `switch` of @ref slice_inner is resolved at
 * compile time and a single closed form is inlined in the device code, see @ref build_r.
 *
 * Entries with @f$ k r > @f$ @ref far_cut_ are written as a plain zero instead of going through
 * the erf/exp calls: the screened kernel decays like @f$ e^{-k r} @f$, so they are below
 * @f$ 4 \cdot 10^{-18} @f$ of the maximum of the very same table (reached at @f$ r = 0 @f$), i.e.
 * under the round-off the length-Q contraction of the slice `dgemm` commits anyway. The table
 * therefore carries an absolute floor of ~1e-17, deliberately -- the same choice as the lattice
 * side, see @ref LGFTwoUnbounded::build_r_.
 *
 * @tparam TYPE the kernel, must satisfy @ref handles
 * @param nx_lo the first (non-negative) index of the closed-form direction held by this rank
 * @param nnx the number of indexes held by this rank, i.e. the leading dimension of R
 * @param kz the (non-negative) wavenumber of the spectral direction
 * @param hx the grid spacing of the closed-form direction
 * @param sig the regularisation length
 */
template <int TYPE>
inline void build_r_(double* r, const double* kq, const int nx_lo, const int nnx, const int nq,
                     const double kz, const double hx, const double sig, stream_t stream) {
    auto op = m_backend_lambda(const int idx) {
        const int    nx  = nx_lo + (idx % nnx);
        const int    q   = idx / nnx;
        const double kq_ = kq[q];
        const double k   = sqrt(kq_ * kq_ + kz * kz);
        const double rx  = ((double)nx) * hx;
        r[idx]           = (k * rx > far_cut_) ? 0.0 : slice_inner(TYPE, rx, k, sig);
    };
    parallel_for(op, nnx * nq, "2unb1spe hej green R", stream);
}

/**
 * @brief fill `R[q * nnx + i] = G1((nx_lo + i) hx; sqrt(k_q^2 + kz^2), sig)` on the device
 *
 * The dispatch on the kernel happens here, on the host, once per slice: @ref build_r_ is then
 * instantiated for a single closed form. The layout is that of a column-major
 * @f$ (n_{nx} \times Q) @f$ matrix of leading dimension @p nnx, i.e. the first operand of the
 * slice `dgemm` (@ref backend_dgemm_nt).
 *
 * @param typeGreen the kernel, must satisfy @ref handles
 */
inline void build_r(double* r, const double* kq, const int typeGreen, const int nx_lo, const int nnx,
                    const int nq, const double kz, const double hx, const double sig, stream_t stream) {
    switch (typeGreen) {
        case HEJ_2:
            build_r_<HEJ_2>(r, kq, nx_lo, nnx, nq, kz, hx, sig, stream);
            break;
        case HEJ_4:
            build_r_<HEJ_4>(r, kq, nx_lo, nnx, nq, kz, hx, sig, stream);
            break;
        case HEJ_6:
            build_r_<HEJ_6>(r, kq, nx_lo, nnx, nq, kz, hx, sig, stream);
            break;
        case HEJ_8:
            build_r_<HEJ_8>(r, kq, nx_lo, nnx, nq, kz, hx, sig, stream);
            break;
        case HEJ_10:
            build_r_<HEJ_10>(r, kq, nx_lo, nnx, nq, kz, hx, sig, stream);
            break;
        default:
            FLUPS_CHECK(false, "[HEJTwoUnbounded] unsupported Green kernel %d", typeGreen);
    }
}

};  // namespace HEJTwoUnbounded

#endif  // HEJ_TWO_UNBOUNDED_HPP_
