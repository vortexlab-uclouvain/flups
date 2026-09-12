#ifndef LGF_ONE_UNBOUNDED_HPP_
#define LGF_ONE_UNBOUNDED_HPP_

// Pick which <complex> header we should use, depending on the compiler
// NVCC can't compile <complex>, and clang can't compile <cuda/std/complex>
#if defined(__NVCC__)
#include <cuda/std/complex>
#else
#include <complex>
#endif

#include <climits>
#include <cmath>
#include <array>
#include <tuple>

#include "backend.hpp"
#include "defines.hpp"

namespace LGFOneUnbounded {

    // ----------------------------------------------------------------------------------------
    // Types
    // ----------------------------------------------------------------------------------------
    using real_t = double;                          //!< Note: floating point literals in this file are double precision.
#if __NVCC__
    using complex_t = cuda::std::complex<real_t>;   //!< complex precision matches real precision 
#else
    using complex_t = std::complex<real_t>;         //!< complex precision matches real precision 
#endif


    // ----------------------------------------------------------------------------------------
    // Utility functions
    // ----------------------------------------------------------------------------------------
    m_constant static const complex_t im = complex_t(0.0, 1.0); //!< constant for i = sqrt(-1)

    /**
     * @brief a square root with complex return type that accepts real or complex arguments
     */
    template <typename T>
    m_device complex_t complex_sqrt(T z) { 
        return sqrt(static_cast<complex_t>(z)); 
    }

    /**
     * @brief the principal cube root, with complex return type, of a real or complex argument
     *
     * This is mathematically equivalent to pow(z, 1./3.), but it does not go through
     * std::pow/std::polar: polar() carries a __glibcxx_assert that expands to a
     * host-only function, which the device compiler rejects inside device code.
     */
    template <typename T>
    m_device complex_t complex_cbrt(T z) {
        const complex_t zc  = static_cast<complex_t>(z);
        const real_t    rho = cbrt(abs(zc));
        const real_t    th  = arg(zc) / 3.0;
        return complex_t(rho * cos(th), rho * sin(th));
    }

    /**
     * @brief evaluate a polynomial with real coefficients via the Horner scheme
     * 
     * @param a an array of length n containing coefficients a_{0} to a_{n-1} 
     * @param n the number of coefficients
     * @param x the evaluation point
     */
    template <typename T>
    m_device T horner(const real_t* a, int n, T x) {
        T out = a[n - 1];
        for (int i = 2; i <= n; i++) {
            out = x*out + a[n - i];
        }
        return out;
    }

    /**
     * @brief evaluate log(x) via series expansion about x = 1
     * 
     * @param x_minus_one the value (x - 1) used in the power series
     */
    m_constant constexpr int lgf_log_series_terms_ = 12;
    m_device real_t log_series(real_t x_minus_one) {
        real_t out = 0.0;
        for (int i = 0; i < lgf_log_series_terms_; i++) {
            out += 1.0 / (lgf_log_series_terms_ - i);
            out *= -x_minus_one;
        }
        return -out;
    }

    /**
     * @brief evaluate (x - 1/x) via series expansion about x = 1
     * 
     * @param x_minus_one the value (x - 1) used in the power series
     */
    m_constant constexpr int lgf_value_minus_inverse_series_terms_ = 12;
    m_device real_t value_minus_inverse_series(real_t x_minus_one) {
        real_t out = 0.0;
        for (int i = 0; i < lgf_value_minus_inverse_series_terms_; i++) {
            out += 1.0;
            out *= -x_minus_one;
        }
        return -out + x_minus_one;
    }

    /**
     * @brief evaluate the root r = lambda - sqrt(lambda^2 - 1), |r| <= 1, of the Mehrstellen stencils, with
     * lambda = -a_0 / (2 a_1), and also r - 1 free of cancellation
     *
     * For a_1 < 0 (always for meh6, small wavenumbers for meh4) lambda >= 1 and r -> 1 as k -> 0. Forming
     * lambda first loses ~ eps / (lambda - 1) in relative precision, amplified by n in r^n (1e-10 at k h ~ 1e-3),
     * so delta = lambda - 1 = -(a_0 + 2 a_1) / (2 a_1) is formed from the y's directly: for both stencils
     * a_0 + 2 a_1 = 4 (y_0 + y_1) - 8/3 y_0 y_1 identically. With s = delta + sqrt(delta (delta + 2)),
     * r - 1 = delta - sqrt(delta (delta + 2)) is evaluated as -2 delta / s and r = 1 / (1 + s), both
     * exact for every delta >= 0. Recovering r as 1 + (r - 1) instead would cancel in the opposite
     * limit: meh4 has a_1 -> 0^- on the curve cos k_1 + cos k_2 = -1 (an endpoint of the quadrature
     * interval at k_z = pi/2), where delta -> +inf and r -> 0, and 1 + (r - 1) then leaves only
     * eps / r relative digits in r (1e-9 at the closest Q = 64 Gauss-Legendre node, which is what
     * put the MEHR_4L / MEHR_4F N = 16 table at 6e-13 on the k_z = pi/2 slice).
     * For a_1 >= 0 (meh4 at large wavenumbers) lambda <= -3 and r = lambda + sqrt(lambda^2 - 1) in (-1/3, 0],
     * with r = 0 for a_1 = 0 (infinite lambda); no cancellation there.
     * Use @ref meh_rpow for r^p and (r - 1) (r + 1) for r^2 - 1.
     */
    m_device void meh_root(const real_t a[2], const real_t y[2], real_t* r, real_t* rm1) {
        if (a[1] < 0.0) {
            const real_t delta = -(2.0 * (y[0] + y[1]) - 4. / 3. * y[0] * y[1]) / a[1];
            const real_t sdel  = delta + sqrt(delta * (delta + 2.0));
            *rm1               = -2.0 * delta / sdel;
            *r                 = 1.0 / (1.0 + sdel);
        } else {
            const real_t lambda = -a[0] / (2.0 * a[1]);
            *r                  = -1.0 / (fabs(lambda) + sqrt(lambda * lambda - 1.0));
            *rm1                = *r - 1.0;
        }
    }

    /**
     * @brief r^p from the output of @ref meh_root: through log1p(r - 1) when r is close to 1
     */
    m_device real_t meh_rpow(const real_t r, const real_t rm1, const int p) {
        return (r > 0.5) ? exp(p * log1p(rm1)) : pow(r, p);
    }

    // ----------------------------------------------------------------------------------------
    // The params / eval split
    // ----------------------------------------------------------------------------------------
    //
    // Every closed form below is a sum of a few terms of the shape (amplitude) x (root)^n: the
    // roots, the prefactors and the branch selection depend on the wavenumber only, while n only
    // enters through r^n (and sin(n theta + phi) for a conjugate pair). Each kernel is therefore
    // split in two:
    //   - `<kernel>_params(c)` returns a small POD struct, evaluated ONCE per wavenumber,
    //   - `<kernel>_eval(params, n)` is the n-dependent part alone.
    // The public `<kernel>(n, c)` entry points are thin wrappers over the pair, so the callers
    // that need a single value are unchanged.
    //
    // The split is a pure hoist: the expressions and their evaluation order are exactly those of
    // the original monolithic forms. In particular r^n keeps going through `pow` / `exp(n log r)`
    // and is NEVER replaced by a recurrence r^n = r * r^(n-1): the n eps growth of a recurrence is
    // what the `meh_root` / `meh_rpow` cancellation-free forms were written to avoid in the first
    // place (it is the mechanism behind the MEHR_6F 1e-12 floor).
    //
    /**
     * @brief which closed form applies to a given wavenumber
     *
     * Each kernel selects one of these once per wavenumber (in its `*_params` function) and the
     * matching arm of its `*_eval` function is then the only n-dependent work left.
     */
    enum Branch : int {
        kInvalid = 0,  //!< c < 0, i.e. an argument error
        kZero,         //!< c = 0 (or k1 = k2 = 0): the operator is singular, the LGF is not decaying
        kSmall,        //!< the root closest to 1 goes through its series expansion
        kReal,         //!< the roots that matter are real and distinct
        kRepeat,       //!< two roots collide, the residue expansion applies
        kImag,         //!< the roots that matter are a conjugate pair
        kGeneric,      //!< the single generic branch of a kernel (lgf6, Mehrstellen)
        kDegenerate    //!< a_1 = 0, the stencil loses its outermost point (meh4 only)
    };

    /**
     * @brief the relative level below which a table entry is not worth evaluating
     */
    m_constant constexpr real_t lgf_trunc_ = 1.0e-17;

    /**
     * @brief the largest index worth evaluating for a contribution decaying like @f$ \rho^n @f$
     *
     * The caller (@ref LGFTwoUnbounded::build_r) writes a plain zero beyond that index instead of
     * calling `eval`. This puts an absolute floor of @ref lgf_trunc_ on the table, of the same
     * class as the round-off the length-Q dot product that follows already commits. Only the
     * decay rate is taken into account, not the amplitude of the contribution: the two are
     * anti-correlated here (an amplitude blows up only when a root approaches the unit circle,
     * i.e. exactly where nothing is truncated).
     *
     * @param log_rho the logarithm of the dominant decay rate; >= 0 (or NaN) means no decay
     * @return the last index to evaluate, at least 1 so that the n = 0 and |n| = 1 special cases
     *         of the Mehrstellen forms always survive
     */
    m_device int decay_cut(const real_t log_rho) {
        if (!(log_rho < 0.0)) {
            return INT_MAX;  // no decay at all, or not a number: never truncate
        }
        const real_t nc = ceil(log(lgf_trunc_) / log_rho);
        return (nc < (real_t)INT_MAX) ? (int)fmax(nc, 1.0) : INT_MAX;
    }

    /**
     * @brief the n-independent data of an oscillatory contribution, see @ref oscillatory_form
     */
    struct OscParams {
        real_t amp;  //!< |K|
        real_t rho;  //!< |r|, i.e. the decay rate
        real_t th;   //!< arg(r)
        real_t phi;  //!< arg(K)
    };

    /**
     * @brief the n-independent part of @ref oscillatory_form
     *
     * @param r one of the complex roots
     * @param pp the complex factor p'(r, c)
     * @param w the width of the stencil
     */
    m_device OscParams oscillatory_params(complex_t r, complex_t pp, int w) {
        const complex_t K = complex_t(0.0, 2.0)*pow(r, w - 1)/pp;
        OscParams o;
        o.amp = abs(K);
        o.rho = abs(r);
        o.th  = arg(r);
        o.phi = arg(K);
        return o;
    }

    /**
     * @brief the n-dependent part of @ref oscillatory_form
     */
    m_device real_t oscillatory_eval(const OscParams& o, int n) {
        return o.amp * pow(o.rho, abs(n)) * sin(o.th*abs(n) + o.phi);
    }

    /**
     * @brief return the oscillatory contribution to the LGF from a conjugate root pair
     * 
     * @param r one of the complex roots
     * @param pp the complex factor p'(r, c)
     * @param w the width of the stencil
     * @param n the evaluation index
     */
    m_device const real_t oscillatory_form(complex_t r, complex_t pp, int w, int n) {
        return oscillatory_eval(oscillatory_params(r, pp, w), n);
    }

    // ----------------------------------------------------------------------------------------
    // Reusable constants
    // ----------------------------------------------------------------------------------------
    // large lambda expansion of r(lambda)

    // root of q(lambda) nearest 1
    m_constant constexpr int lgf_small_c_terms_ = 8;
    m_constant static const real_t lgf4_small_c_expansion_[8] = {0, 0.5, 1./24., 1./144., 5./3456., 7./20736., 7./82944., 11./497664.};
    m_constant static const real_t lgf6_small_c_expansion_[8] = {0, 0.5, 1./24., 1./720., -1./1152., -149./518400., -259./6220800., 163./62208000.};
    m_constant static const real_t lgf8_small_c_expansion_[8] = {0, 0.5, 1./24., 1./720., 1./40320., 577./3628800., 389./6220800., 34987./3048192000.};

    // p'(z, c), template to allow real_t or complex
    template <typename T> 
    m_device T lgf4_ppzc_(T z, real_t c) { return -4./3. + z*((5. + 2.*c) + z*(-4. + z/3.)); }
    template <typename T> 
    m_device T lgf6_ppzc_(T z, real_t c) { return 3./20. + z*(-3. + z*((49./6. + 3.*c) + z*(-6. + z*(3./4. - z/15.)))); }
    template <typename T> 
    m_device T lgf8_ppzc_(T z, real_t c) { return -8./315. + z*(2./5. + z*(-24./5. + z*(205./18. + 4.*c + z*(-8. + z*(6./5. + z*(-8./45. + z/70.)))))); } 

    // LGF4 repeated root residues
    m_constant static const real_t lgf4_c_star_ = 3.0;
    m_constant const real_t sq15 = 3.872983346207417;   // cannot evaluate using sqrt(15.) directly
    m_constant const real_t lgf4_repeat_g0_[2]  = {4./(5.*sq15), 1./5.};
    m_constant const real_t lgf4_repeat_g1_[4]  = {-14./(75.*sq15), -2./25., -4./(25.*sq15), -1./150.};
    m_constant const real_t lgf4_repeat_g2_[6]  = {901./(18750.*sq15), 74./3125., 253./(3750.*sq15), 23./3750., 1./(250.*sq15), 1./15000.};
    m_constant const real_t lgf4_repeat_g3_[8]  = {-37313./(2812500.*sq15), -11267./1640625., -3167./(140625.*sq15), -1./375., -31./(11250.*sq15), -31./281250., -1./(28125.*sq15), -1./3150000.};
    m_constant const real_t lgf4_repeat_g4_[10] = {3885631./(1012500000.*sq15), 1194829./590625000., 1118351./(157500000.*sq15), 12571./13125000., 13751./(11250000.*sq15), 127./1875000., 247./(6750000.*sq15), 13./15750000., 1./(6300000.*sq15), 1./1134000000.};
    
    // LGF8 repeated root residues
    m_constant static const real_t lgf8_c_star_ = 3.204471924659902;

    // The constants a_jk, u_jk, and vj_k appearing in the repeated root expansion of the order eight stencil. 
    // The naming convention allows access to constants via lgf8_a_j_[k], lgf8_u_j_[k], and lgf8_v_j_[k]
    m_constant const real_t lgf8_a_0_[2]  = {-0.1932304295173025, -9.940958723693284e-2};
    m_constant const real_t lgf8_a_1_[4]  = {5.528645180523051e-2, 6.381749748906844e-2, 1.920895723993029e-2, 1.647044339102893e-3};
    m_constant const real_t lgf8_a_2_[6]  = {-1.246195963786745e-2, -2.200316066032414e-2, -1.156690496243602e-2, -2.513956535483380e-3, -2.386943138091703e-4, -8.186599895557273e-6};
    m_constant const real_t lgf8_a_3_[8]  = {2.586595339831205e-3, 5.946351585490818e-3, 4.392324394119772e-3, 1.460419122054276e-3, 2.498408821419240e-4, 2.281049273902400e-5, 1.054600142736557e-6, 1.937682182122062e-8};
    m_constant const real_t lgf8_a_4_[10] = {-6.143387386724980e-4, -1.461770070530813e-3, -1.338631620906462e-3, -5.893614275050890e-4, -1.419718637870522e-4, -2.006233681173800e-5, -1.705558618042543e-6, -8.568226736629197e-8, -2.340119751996688e-9, -2.675335915571021e-11};

    m_constant const real_t lgf8_u_0_[1] = {0.0015031218476240784};
    m_constant const real_t lgf8_u_1_[2] = {-0.013198947960911374, -0.00309497626043345};
    m_constant const real_t lgf8_u_2_[3] = {0.002223792419320104, 0.0003211726758356002, -3.489943855291188e-6};
    m_constant const real_t lgf8_u_3_[4] = {8.218756082933073e-5, 0.00012939734109457704, 2.076196709989474e-5, 7.970739601754506e-7};
    m_constant const real_t lgf8_u_4_[5] = {-0.00011284539233578474, -4.455888056601662e-5, -4.596415525190604e-6, -1.2515613302441402e-7, 7.499072573978986e-10};

    m_constant const real_t lgf8_v_0_[1] = {-0.07869060869097218};
    m_constant const real_t lgf8_v_1_[2] = {0.005694070346472886, -0.00011828157312623746};
    m_constant const real_t lgf8_v_2_[3] = {0.001748394270280815, 0.0007853940979982054, 6.08419435502199e-5};
    m_constant const real_t lgf8_v_3_[4] = {-0.000574565978306078, -0.000152486888729598, -8.030910109571855e-6, 6.101311014540467e-8};
    m_constant const real_t lgf8_v_4_[5] = {4.706634682161246e-5, -8.921284296708267e-6, -4.142946014756904e-6, -3.42884013573887e-7, -7.828815619884859e-9};

    // ----------------------------------------------------------------------------------------
    // Symbols
    // ----------------------------------------------------------------------------------------
    m_device real_t lgf2_symbol(real_t k) { return 4.0 * pow(sin(k/2.), 2); }
    m_device real_t lgf4_symbol(real_t k) { return 4.0 * (4./3.*pow(sin(k/2.), 2) - 1./12.*pow(sin(k), 2)); }
    m_device real_t lgf6_symbol(real_t k) { return 4.0 * (3./2.*pow(sin(k/2.), 2) - 3./20.*pow(sin(k), 2) + 1./90.*pow(sin(1.5*k), 2)); }
    m_device real_t lgf8_symbol(real_t k) { return 4.0 * (8./5.*pow(sin(k/2.), 2) - 1./5. *pow(sin(k), 2) + 8./315.*pow(sin(1.5*k), 2) - 1./560.*pow(sin(2.*k), 2)); }

    // ----------------------------------------------------------------------------------------
    // LGF2
    // ----------------------------------------------------------------------------------------
    m_constant const real_t lgf2_small_cutoff_ = 1.0e-3;

    /**
     * @brief the n-independent data of @ref lgf2
     */
    struct Lgf2Params {
        int    branch;  //!< which closed form applies, see @ref Branch
        int    n_cut;   //!< the last index worth evaluating, see @ref decay_cut
        real_t r;       //!< kReal: the root
        real_t logr;    //!< kSmall: log(r), through the series
        real_t den;     //!< r - 1/r
    };

    /**
     * @brief the root and the prefactor of @ref lgf2, for one value of c
     */
    m_device Lgf2Params lgf2_params(real_t c) {
        Lgf2Params p = {kInvalid, INT_MAX, 0.0, 0.0, 0.0};
        if (c < 0) {
            assert(false && "[LGFOneUnbounded] Argument error: c must be nonnegative");
        } else if (c == 0.0) {
            p.branch = kZero;
        } else if (c < lgf2_small_cutoff_) {
            const real_t rm1 = c/2. - sqrt(c*c/4. + c);
            p.branch = kSmall;
            p.logr   = log_series(rm1);
            p.den    = value_minus_inverse_series(rm1);
            p.n_cut  = decay_cut(p.logr);
        } else {
            const real_t r = 1 + c/2. - sqrt(c*c/4. + c);
            p.branch = kReal;
            p.r      = r;
            p.den    = r - 1./r;
            p.n_cut  = decay_cut(log(fabs(r)));
        };
        return p;
    }

    /**
     * @brief the n-dependent part of @ref lgf2
     */
    m_device real_t lgf2_eval(const Lgf2Params& p, int n) {
        if (p.branch == kZero) {
            return -0.5 * abs(n);
        } else if (p.branch == kSmall) {
            return -exp(abs(n) * p.logr) / p.den;
        } else if (p.branch == kReal) {
            return -pow(p.r, abs(n)) / p.den;
        };
        return 0.0;
    }

    m_device real_t lgf2(int n, real_t c) {
        return lgf2_eval(lgf2_params(c), n);
    }

    m_device real_t lgf2(int n, real_t k1, real_t k2) {
        return lgf2(n, lgf2_symbol(k1) + lgf2_symbol(k2));
    }

    // ----------------------------------------------------------------------------------------
    // LGF4
    // ----------------------------------------------------------------------------------------

    m_constant constexpr real_t lgf4_small_cutoff_ = 1e-3;
    m_constant constexpr real_t lgf4_repeat_cutoff_ = 1e-5;

    /**
     * @brief the n-independent data of @ref lgf4
     */
    struct Lgf4Params {
        int       branch;  //!< which closed form applies, see @ref Branch
        int       n_cut;   //!< the last index worth evaluating, see @ref decay_cut
        real_t    logr1;   //!< kSmall: log(r1), through the series
        real_t    pp1;     //!< kSmall: p'(r1, c)
        real_t    r1;      //!< kReal: the larger real root
        real_t    r2;      //!< kSmall/kReal: the smaller real root
        real_t    pp2;     //!< kSmall: p'(r2, c)
        real_t    d1, d2;  //!< kReal: the two denominators
        real_t    dr;      //!< kReal: r1 - r2
        real_t    dc;      //!< kRepeat: c - c*
        OscParams osc;     //!< kImag: the conjugate pair
    };

    /**
     * @brief the n-independent data of lgf4(n, c) for c = 0
     */
    m_device void lgf4_c_zero_params_(Lgf4Params* p) {
        p->branch = kZero;
        p->n_cut  = INT_MAX;  // the -|n|/2 term does not decay
    }

    /** 
     * @brief evaluate lgf4(n, c) for c = 0 
     */
    m_device real_t lgf4_c_zero_eval_(const Lgf4Params& p, int n) {
        const real_t r = 7.179676972449123e-2;
        const real_t K = -7.216878364870408e-2;
        return -0.5*abs(n) + K*(pow(r, abs(n)) - 1.0);
    }    

    /** 
     * @brief the n-independent data of lgf4(n, c) for 0 < c < lgf4_small_c_cutoff_
     */
    m_device void lgf4_c_small_params_(Lgf4Params* p, real_t c) {
        // larger real root
        const real_t lambdam1 = horner(lgf4_small_c_expansion_, lgf_small_c_terms_, c);
        const real_t rm1 = lambdam1 - sqrt(lambdam1*lambdam1 + 2*lambdam1);
        const real_t r1 = rm1 + 1.0;
        // smaller real root
        const real_t lambda2 = 4.0 + sqrt(9.0 - 3.0*c);
        const real_t r2 = lambda2 - sqrt(lambda2*lambda2 - 1.0);
        // the prefactors
        const real_t r1_minus_r1_inv = value_minus_inverse_series(rm1);
        const real_t logr1 = log_series(rm1);
        p->branch = kSmall;
        p->logr1  = logr1;
        p->pp1    = 1./12. * r1_minus_r1_inv * (r1 - r2) * (r1 - 1.0/r2);
        p->r2     = r2;
        p->pp2    = lgf4_ppzc_(r2, c);
        p->n_cut  = decay_cut(fmax(logr1, log(fabs(r2))));
    }

    /**
     * @brief evaluate lgf4(n, c) for 0 < c < lgf4_small_c_cutoff_
     */
    m_device real_t lgf4_c_small_eval_(const Lgf4Params& p, int n) {
        return exp((abs(n) + 1)*p.logr1) / p.pp1 + pow(p.r2, abs(n) + 1) / p.pp2;
    }

    /**
     * @brief the n-independent data of lgf4(n, c) for lgf4_small_cutoff_ < c < (lgf4_c_star_ - lgf4_repeat_cutoff_)
     */
    m_device void lgf4_c_real_params_(Lgf4Params* p, real_t c) {
        const real_t lambda1 = 4.0 + sqrt(9.0 - 3.0*c);
        const real_t lambda2 = 4.0 - sqrt(9.0 - 3.0*c);
        const real_t r1 = lambda1 - sqrt(lambda1*lambda1 - 1.0);
        const real_t r2 = lambda2 - sqrt(lambda2*lambda2 - 1.0);
        p->branch = kReal;
        p->r1     = r1;
        p->r2     = r2;
        p->d1     = (r1 - 1.0/r1)*(r1 - 1.0/r2);
        p->d2     = (r2 - 1.0/r1)*(r2 - 1.0/r2);
        p->dr     = r1 - r2;
        p->n_cut  = decay_cut(fmax(log(fabs(r1)), log(fabs(r2))));
    }

    /**
     * @brief evaluate lgf(n, c) for lgf4_small_cutoff_ < c < (lgf4_c_star_ - lgf4_repeat_cutoff_)
     */
    m_device real_t lgf4_c_real_eval_(const Lgf4Params& p, int n) {
        const real_t f1 = 12*pow(p.r1, abs(n) + 1)/p.d1;
        const real_t f2 = 12*pow(p.r2, abs(n) + 1)/p.d2;
        return (f1 - f2) / p.dr;
    }

    /**
     * @brief the n-independent data of lgf4(n, c) for abs(c - lgf4_c_star_) < lgf4_repeat_cutoff_
     */
    m_device void lgf4_c_repeat_params_(Lgf4Params* p, real_t c) {
        p->branch = kRepeat;
        p->dc     = c - lgf4_c_star_;
        // the residue expansion carries a polynomial in n on top of r^n, which the cut ignores
        p->n_cut  = decay_cut(log(4.0 - sqrt(15.0)));
    }

    /**
     * @brief evaluate lgf(n, c) for abs(c - lgf4_c_star_) < lgf4_repeat_cutoff_
     */
    m_device real_t lgf4_c_repeat_eval_(const Lgf4Params& p, int n) {
        const real_t r = 4.0 - sqrt(15.0);
        const real_t nr = static_cast<real_t>(abs(n));
        const real_t lgf4_g[5] = {
            horner(lgf4_repeat_g0_, 2,  nr),
            horner(lgf4_repeat_g1_, 4,  nr),
            horner(lgf4_repeat_g2_, 6,  nr),
            horner(lgf4_repeat_g3_, 8,  nr),
            horner(lgf4_repeat_g4_, 10, nr)
        };
        return pow(r, abs(n)) * horner(lgf4_g, 5, p.dc);
    }

    /**
     * @brief the n-independent data of lgf4(n, c) for (lgf4_c_star_ + lgf4_repeat_cutoff_) < c
     */
    m_device void lgf4_c_imag_params_(Lgf4Params* p, real_t c) {
        const complex_t lambda = 4.0 - complex_sqrt(9.0 - 3.0*c);
        const complex_t r = lambda - sqrt(lambda*lambda - 1.0);
        const complex_t pp = im/6.*imag(r)*(r - 1./r)*(r - 1./conj(r)); // lgf4_ppzc_(r, c);
        p->branch = kImag;
        p->osc    = oscillatory_params(r, pp, 2);
        p->n_cut  = decay_cut(log(p->osc.rho));
    }

    /**
     * @brief evaluate lgf4(n, c) for (lgf4_c_star_ + lgf4_repeat_cutoff_) < c
     */
    m_device real_t lgf4_c_imag_eval_(const Lgf4Params& p, int n) {
        return oscillatory_eval(p.osc, n);
    }

    /**
     * @brief the roots and the prefactors of @ref lgf4, for one value of c
     */
    m_device Lgf4Params lgf4_params(real_t c) {
        Lgf4Params p = {};
        p.branch = kInvalid;
        p.n_cut  = INT_MAX;
        if (c < 0) {
            assert(false && "[LGFOneUnbounded] Argument error: c must be nonnegative");
        } else if (c == 0.0) {
            lgf4_c_zero_params_(&p);
        } else if (c < lgf4_small_cutoff_) {
            lgf4_c_small_params_(&p, c);
        } else if (c < lgf4_c_star_ - lgf4_repeat_cutoff_) {
            lgf4_c_real_params_(&p, c);
        } else if (c < lgf4_c_star_ + lgf4_repeat_cutoff_) {
            lgf4_c_repeat_params_(&p, c);
        } else {
            lgf4_c_imag_params_(&p, c);
        };
        return p;
    }

    /**
     * @brief the n-dependent part of @ref lgf4
     */
    m_device real_t lgf4_eval(const Lgf4Params& p, int n) {
        switch (p.branch) {
            case kZero:   return lgf4_c_zero_eval_(p, n);
            case kSmall:  return lgf4_c_small_eval_(p, n);
            case kReal:   return lgf4_c_real_eval_(p, n);
            case kRepeat: return lgf4_c_repeat_eval_(p, n);
            case kImag:   return lgf4_c_imag_eval_(p, n);
            default:      return 0.0;
        }
    }

    m_device real_t lgf4(int n, real_t c) {
        return lgf4_eval(lgf4_params(c), n);
    }

    m_device real_t lgf4(int n, real_t k1, real_t k2) {
        return lgf4(n, lgf4_symbol(k1) + lgf4_symbol(k2));
    }

    // ----------------------------------------------------------------------------------------
    // LGF6
    // ----------------------------------------------------------------------------------------

    /**
     * @brief return the roots of p(z, c) for the LGF6 stencil
     * 
     * Because r2 and r3 are always a conjugate pair, r3 is omitted from the return.
     */
    m_device std::tuple<real_t, complex_t> lgf6_roots_(real_t c) {
        // real root
        const real_t xi = cbrt(360.*c - 775. + 60.*sqrt(36.*c*c - 155.*c + 405.));
        const real_t lambda1 = 0.25 * (9. - 95./xi + xi);
        const real_t r1 = lambda1 - sqrt(lambda1*lambda1 - 1.0);
        // complex root
        const complex_t cbrt_one = complex_t(-0.5, sqrt(3.)/2.);
        const complex_t lambda2 = 0.25 * (9. - 95./(cbrt_one * xi) + (cbrt_one * xi));
        const complex_t r2 = lambda2 - sqrt(lambda2*lambda2 - 1.0);
        // return as a tuple so that each is typed correctly
        return std::make_tuple(r1, r2);
    }

    m_constant constexpr real_t lgf6_small_cutoff_ = 1e-3;

    /**
     * @brief the n-independent data of @ref lgf6
     */
    struct Lgf6Params {
        int       branch;  //!< which closed form applies, see @ref Branch
        int       n_cut;   //!< the last index worth evaluating, see @ref decay_cut
        real_t    logr1;   //!< kSmall: log(r1), through the series
        real_t    r1;      //!< kGeneric: the real root
        real_t    pp1;     //!< the prefactor of the real root
        OscParams osc;     //!< the conjugate pair
    };

    /**
     * @brief the n-independent data of lgf6(n, c) for c = 0
     */
    m_device void lgf6_c_zero_params_(Lgf6Params* p) {
        p->branch = kZero;
        p->n_cut  = INT_MAX;  // the -|n|/2 term does not decay
    }

    /**
     * @brief evaluate lgf6(n, c) for c = 0
     */
    m_device real_t lgf6_c_zero_eval_(const Lgf6Params& p, int n) {
        const real_t a = +1.035017642406219e-1;
        const real_t pr = +9.54295097909627e-2;
        const real_t t = -9.958762687216868e-1;
        const real_t f = -2.168575565241744e+0;
        return -0.5*abs(n) + a*pow(pr, abs(n))*sin(abs(n)*t + f) - a*sin(f);
    }

    /**
     * @brief the n-independent data of lgf6(n, c) for 0 < c < lgf6_small_cutoff_
     */
    m_device void lgf6_c_small_params_(Lgf6Params* p, real_t c) {
        // larger real root
        const real_t lambdam1 = horner(lgf6_small_c_expansion_, lgf_small_c_terms_, c);
        const real_t rm1 = lambdam1 - sqrt(lambdam1*lambdam1 + 2*lambdam1);
        const real_t r1 = rm1 + 1.0;
        // conjugate root pair
        const auto roots = lgf6_roots_(c);
        const complex_t r2 = std::get<1>(roots);
        const complex_t pp2 = lgf6_ppzc_(r2, c);
        // the prefactors
        const real_t r1_minus_r1_inv = value_minus_inverse_series(rm1);
        const real_t logr1 = log_series(rm1);
        p->branch = kSmall;
        p->logr1  = logr1;
        p->pp1    = -1./90. * r1_minus_r1_inv * norm((r1 - r2) * (r1 - 1.0/r2));
        p->osc    = oscillatory_params(r2, pp2, 3);
        p->n_cut  = decay_cut(fmax(logr1, log(p->osc.rho)));
    }

    /**
     * @brief evaluate lgf6(n, c) for 0 < c < lgf6_small_cutoff_
     */
    m_device real_t lgf6_c_small_eval_(const Lgf6Params& p, int n) {
        return exp((abs(n) + 2)*p.logr1) / p.pp1 + oscillatory_eval(p.osc, n);
    }

    /**
     * @brief the n-independent data of lgf6(n, c) for lgf6_small_cutoff_ < c
     */
    m_device void lgf6_c_generic_params_(Lgf6Params* p, real_t c) {
        const auto roots = lgf6_roots_(c);
        const real_t r1 = std::get<0>(roots);
        const real_t pp1 = lgf6_ppzc_(r1, c);
        const complex_t r2 = std::get<1>(roots);
        const complex_t pp2 = lgf6_ppzc_(r2, c);
        p->branch = kGeneric;
        p->r1     = r1;
        p->pp1    = pp1;
        p->osc    = oscillatory_params(r2, pp2, 3);
        p->n_cut  = decay_cut(fmax(log(fabs(r1)), log(p->osc.rho)));
    }

    /**
     * @brief evaluate lgf6(n, c) for lgf6_small_cutoff_ < c
     */
    m_device real_t lgf6_c_generic_eval_(const Lgf6Params& p, int n) {
        return pow(p.r1, abs(n) + 2)/p.pp1 + oscillatory_eval(p.osc, n);
    }

    /**
     * @brief the roots and the prefactors of @ref lgf6, for one value of c
     */
    m_device Lgf6Params lgf6_params(real_t c) {
        Lgf6Params p = {};
        p.branch = kInvalid;
        p.n_cut  = INT_MAX;
        if (c < 0) {
            assert(false && "[LGFOneUnbounded] Argument error: c must be nonnegative");
        } else if (c == 0.0) {
            lgf6_c_zero_params_(&p);
        } else if (c < lgf6_small_cutoff_) {
            lgf6_c_small_params_(&p, c);
        } else {
            lgf6_c_generic_params_(&p, c);
        };
        return p;
    }

    /**
     * @brief the n-dependent part of @ref lgf6
     */
    m_device real_t lgf6_eval(const Lgf6Params& p, int n) {
        switch (p.branch) {
            case kZero:    return lgf6_c_zero_eval_(p, n);
            case kSmall:   return lgf6_c_small_eval_(p, n);
            case kGeneric: return lgf6_c_generic_eval_(p, n);
            default:       return 0.0;
        }
    }

    m_device real_t lgf6(int n, real_t c) {
        return lgf6_eval(lgf6_params(c), n);
    }

    m_device real_t lgf6(int n, real_t k1, real_t k2) {
        return lgf6(n, lgf6_symbol(k1) + lgf6_symbol(k2));
    }
    
    // ----------------------------------------------------------------------------------------
    // LGF8
    // ----------------------------------------------------------------------------------------

    /**
     * @brief return the roots of p(z, c) for LGF8. 
     * 
     * Because r3 and r4 are always a conjugate pair, r4 is omitted
     */
    m_device std::array<complex_t, 3> lgf8_roots_(real_t c) {
        const real_t tmp1 = 1968941520. + c*(-879221700. + c*(223915104. - c*44089920));
        const complex_t xi = complex_cbrt(1520225. - 273420.*c + 35.*complex_sqrt(tmp1));
        const real_t del = real(xi + (-4655. + 3780.*c)/xi);
        const complex_t eta = complex_sqrt(del/9. - 434./81.);
        const complex_t lambda1 = 0.5 * (32./9. + eta - complex_sqrt(-868./81. + 81088./(729.*eta) - del/9.));
        const complex_t lambda2 = 0.5 * (32./9. + eta + complex_sqrt(-868./81. + 81088./(729.*eta) - del/9.));
        const complex_t lambda3 = 0.5 * (32./9. - eta - complex_sqrt(-868./81. - 81088./(729.*eta) - del/9.));
        return {
            lambda1 - sqrt(lambda1 + 1.0)*sqrt(lambda1 - 1.0),
            lambda2 - sqrt(lambda2 + 1.0)*sqrt(lambda2 - 1.0),
            lambda3 - sqrt(lambda3 + 1.0)*sqrt(lambda3 - 1.0)
        };
    }

    m_constant constexpr real_t lgf8_small_cutoff_ = 1e-3;
    m_constant constexpr real_t lgf8_repeat_cutoff_ = 1e-5;

    /**
     * @brief the n-independent data of @ref lgf8
     */
    struct Lgf8Params {
        int       branch;      //!< which closed form applies, see @ref Branch
        int       n_cut;       //!< the last index worth evaluating, see @ref decay_cut
        real_t    logr1;       //!< kSmall: log(r1), through the series
        real_t    r1, r2;      //!< the real roots
        real_t    pp1, pp2;    //!< kSmall: p'(r1, c) and p'(r2, c)
        real_t    d1, d2;      //!< kReal: the two denominators
        real_t    dr;          //!< kReal: r2 - r1
        real_t    dc;          //!< kRepeat: c - c*
        OscParams osc1, osc3;  //!< the conjugate pairs (osc1 is kImag only)
    };

    /**
     * @brief the n-independent data of lgf8(n, c) for c = 0
     */
    m_device void lgf8_c_zero_params_(Lgf8Params* p) {
        p->branch = kZero;
        p->n_cut  = INT_MAX;  // the -|n|/2 term does not decay
    }

    /**
     * @brief evaluate lgf8(n, c) for c = 0
     */
    m_device real_t lgf8_c_zero_eval_(const Lgf8Params& p, int n) {
        const real_t rho1 = +1.2185576475724652e-01;
        const real_t th1  = +1.4667413886166037e+00;
        const real_t a1   = +8.2123351081853858e-02;
        const real_t phi1 = -6.5666289805823130e-01;
        const real_t r2   = +9.6191214171019512e-02;
        const real_t a2   = +4.0692417408429382e-02;
        return -0.5*abs(n) + a1*pow(rho1, abs(n))*sin(abs(n)*th1 + phi1) - a1*sin(phi1) + a2*(1.0 - pow(r2, abs(n)));
    }

    /**
     * @brief the n-independent data of lgf8(n, c) for 0 < c < lgf8_small_cutoff_
     */
    m_device void lgf8_c_small_params_(Lgf8Params* p, real_t c) {
        // larger real root
        const real_t lambdam1 = horner(lgf8_small_c_expansion_, lgf_small_c_terms_, c);
        const real_t rm1 = lambdam1 - sqrt(lambdam1*lambdam1 + 2*lambdam1);
        const real_t r1 = rm1 + 1.0;
        // second real root and comlex pair
        const auto roots = lgf8_roots_(c);
        const real_t r2 = real(roots[1]);
        const real_t pp2 = lgf8_ppzc_(r2, c);
        const complex_t r3 = roots[2];
        const complex_t pp3 = lgf8_ppzc_(r3, c);
        // the prefactors
        const real_t r1_minus_r1_inv = value_minus_inverse_series(rm1);
        const real_t logr1 = log_series(rm1);
        p->branch = kSmall;
        p->logr1  = logr1;
        p->pp1    = 1./560. * r1_minus_r1_inv * (r1 - r2)*(r1 - 1/r2)*norm((r1 - r3)*(r1 - 1.0/r3));
        p->r2     = r2;
        p->pp2    = pp2;
        p->osc3   = oscillatory_params(r3, pp3, 4);
        p->n_cut  = decay_cut(fmax(fmax(logr1, log(fabs(r2))), log(p->osc3.rho)));
    }

    /**
     * @brief evaluate lgf8(n, c) for 0 < c < lgf8_small_cutoff_
     */
    m_device real_t lgf8_c_small_eval_(const Lgf8Params& p, int n) {
        return exp((abs(n) + 3)*p.logr1)/p.pp1 + pow(p.r2, abs(n) + 3)/p.pp2 + oscillatory_eval(p.osc3, n);
    }

    /**
     * @brief the n-independent data of lgf8(n, c) for lgf8_small_cutoff_ < c < (lgf8_c_star_ - lgf8_repeat_cutoff_)
     */
    m_device void lgf8_c_real_params_(Lgf8Params* p, real_t c) {
        const auto roots = lgf8_roots_(c);
        const real_t r1 = real(roots[0]);
        const real_t r2 = real(roots[1]);
        const complex_t r3 = roots[2];
        const complex_t pp3 = lgf8_ppzc_(r3, c);
        p->branch = kReal;
        p->r1     = r1;
        p->r2     = r2;
        p->d1     = (r1 - 1./r1)*(r1 - 1./r2)*norm((r1 - r3)*(r1 - 1./r3));
        p->d2     = (r2 - 1./r1)*(r2 - 1./r2)*norm((r2 - r3)*(r2 - 1./r3));
        p->dr     = r2 - r1;
        p->osc3   = oscillatory_params(r3, pp3, 4);
        p->n_cut  = decay_cut(fmax(fmax(log(fabs(r1)), log(fabs(r2))), log(p->osc3.rho)));
    }

    /**
     * @brief evaluate lgf8(n, c) for lgf8_small_cutoff_ < c < (lgf8_c_star_ - lgf8_repeat_cutoff_)
     */
    m_device real_t lgf8_c_real_eval_(const Lgf8Params& p, int n) {
        const real_t f1 = 560. * pow(p.r1, abs(n) + 3) / p.d1;
        const real_t f2 = 560. * pow(p.r2, abs(n) + 3) / p.d2;
        return (f2 - f1) / p.dr + oscillatory_eval(p.osc3, n);
    }

    /**
     * @brief the n-independent data of lgf8(n, c) for abs(c - lgf8_c_star_) < lgf8_repeat_cutoff_
     */
    m_device void lgf8_c_repeat_params_(Lgf8Params* p, real_t c) {
        p->branch = kRepeat;
        p->dc     = c - lgf8_c_star_;
        // the residue expansion carries a polynomial in n on top of r^n, which the cut ignores
        p->n_cut  = decay_cut(log(1.4016094393066413e-01));
    }

    /**
     * @brief evaluate lgf8(n, c) for abs(c - lgf8_c_star_) < lgf8_repeat_cutoff_
     */
    m_device real_t lgf8_c_repeat_eval_(const Lgf8Params& p, int n) {
        const real_t r0 = 1.4016094393066413e-01;
        const real_t rho1 = 1.2718191072648583e-01;
        const real_t theta1 = 1.5912866598657263e+00;
        const real_t r0n = pow(r0, abs(n));
        const real_t p1n = pow(rho1, abs(n));
        const real_t cnt1 = cos(abs(n)*theta1);
        const real_t snt1 = sin(abs(n)*theta1);
        const real_t nr = abs(n);
        const real_t lgf8_g[5] = {
            r0n * horner(lgf8_a_0_, 2,  nr) + p1n*(cnt1*horner(lgf8_u_0_, 1, nr) + snt1*horner(lgf8_v_0_, 1, nr)),
            r0n * horner(lgf8_a_1_, 4,  nr) + p1n*(cnt1*horner(lgf8_u_1_, 2, nr) + snt1*horner(lgf8_v_1_, 2, nr)),
            r0n * horner(lgf8_a_2_, 6,  nr) + p1n*(cnt1*horner(lgf8_u_2_, 3, nr) + snt1*horner(lgf8_v_2_, 3, nr)),
            r0n * horner(lgf8_a_3_, 8,  nr) + p1n*(cnt1*horner(lgf8_u_3_, 4, nr) + snt1*horner(lgf8_v_3_, 4, nr)),
            r0n * horner(lgf8_a_4_, 10, nr) + p1n*(cnt1*horner(lgf8_u_4_, 5, nr) + snt1*horner(lgf8_v_4_, 5, nr))
        };
        return -horner(lgf8_g, 5, p.dc);
    }

    /**
     * @brief the n-independent data of lgf8(n, c) for (lgf8_c_star_ + lgf8_repeat_cutoff_) < c
     */
    m_device void lgf8_c_imag_params_(Lgf8Params* p, real_t c) {
        const auto roots = lgf8_roots_(c);
        const complex_t r1 = roots[1];
        const complex_t r3 = roots[2];
        const complex_t cr1 = conj(r1);
        const complex_t cr3 = conj(r3);
        const complex_t pp1 = im/280.*imag(r1)*(r1 - 1./r1)*(r1 - 1./cr1)*(r1 - r3)*(r1 - 1./r3)*(r1 - cr3)*(r1 - 1./cr3); // lgf8_ppzc_(r1, c);
        const complex_t pp3 = lgf8_ppzc_(r3, c);
        p->branch = kImag;
        p->osc1   = oscillatory_params(r1, pp1, 4);
        p->osc3   = oscillatory_params(r3, pp3, 4);
        p->n_cut  = decay_cut(fmax(log(p->osc1.rho), log(p->osc3.rho)));
    }

    /**
     * @brief evaluate lgf8(n, c) for (lgf8_c_star_ + lgf8_repeat_cutoff_) < c
     */
    m_device real_t lgf8_c_imag_eval_(const Lgf8Params& p, int n) {
        return oscillatory_eval(p.osc1, n) + oscillatory_eval(p.osc3, n);
    }

    /**
     * @brief the roots and the prefactors of @ref lgf8, for one value of c
     */
    m_device Lgf8Params lgf8_params(real_t c) {
        Lgf8Params p = {};
        p.branch = kInvalid;
        p.n_cut  = INT_MAX;
        if (c < 0) {
            assert(false && "[LGFOneUnbounded] Argument error: c must be nonnegative");
        } else if (c == 0.0) {
            lgf8_c_zero_params_(&p);
        } else if (c < lgf8_small_cutoff_) {
            lgf8_c_small_params_(&p, c);
        } else if (c < lgf8_c_star_ - lgf8_repeat_cutoff_) {
            lgf8_c_real_params_(&p, c);
        } else if (c < lgf8_c_star_ + lgf8_repeat_cutoff_) {
            lgf8_c_repeat_params_(&p, c);
        } else {
            lgf8_c_imag_params_(&p, c);
        };
        return p;
    }

    /**
     * @brief the n-dependent part of @ref lgf8
     */
    m_device real_t lgf8_eval(const Lgf8Params& p, int n) {
        switch (p.branch) {
            case kZero:   return lgf8_c_zero_eval_(p, n);
            case kSmall:  return lgf8_c_small_eval_(p, n);
            case kReal:   return lgf8_c_real_eval_(p, n);
            case kRepeat: return lgf8_c_repeat_eval_(p, n);
            case kImag:   return lgf8_c_imag_eval_(p, n);
            default:      return 0.0;
        }
    }

    m_device real_t lgf8(int n, real_t c) {
        return lgf8_eval(lgf8_params(c), n);
    }

    m_device real_t lgf8(int n, real_t k1, real_t k2) {
        return lgf8(n, lgf8_symbol(k1) + lgf8_symbol(k2));
    }

    // ----------------------------------------------------------------------------------------
    // Mehrstellen: the shared params machinery
    // ----------------------------------------------------------------------------------------
    /**
     * @brief the n-independent data of any of the four Mehrstellen closed forms
     *
     * All four have the same shape, `num * r^(|n| - shift) / den` plus a couple of corrections at
     * |n| <= 1, so one struct covers them; the unused members simply stay at zero.
     */
    struct MehParams {
        int    branch;      //!< kZero (k1 = k2 = 0), kDegenerate (a1 = 0, meh4 only) or kGeneric
        int    n_cut;       //!< the last index worth evaluating, see @ref decay_cut
        bool   use_log;     //!< r^p goes through exp(p log1p(r - 1)), see @ref meh_rpow
        real_t r;           //!< the root of the stencil, see @ref meh_root
        real_t logr;        //!< log1p(r - 1), only when @ref MehParams::use_log
        real_t num;         //!< the n-independent numerator: r/a1 (left) or p(r) (full)
        real_t den;         //!< the n-independent denominator
        real_t v0;          //!< the value at n = 0 when it is a special case (meh4_full)
        real_t c1;          //!< the correction at |n| = 1 (meh6_full)
        real_t c0num;       //!< the numerator of the correction at n = 0 (meh6_full)
        real_t c0den;       //!< its denominator, kept apart so that the order of operations holds
        real_t a0;          //!< kDegenerate: a_0
        real_t b0a0, b1a0;  //!< kDegenerate: b_0/a_0 and b_1/a_0
    };

    /**
     * @brief @ref meh_rpow with the branch and the logarithm resolved once and for all
     */
    m_device real_t meh_rpow_(const MehParams& p, const int q) {
        return p.use_log ? exp(q * p.logr) : pow(p.r, q);
    }

    /**
     * @brief fill the root-dependent part of @ref MehParams, i.e. @ref meh_root and @ref meh_rpow
     *
     * @param[out] p the params being built
     * @param[out] rm1 the value r - 1, needed by the callers to form their denominators
     */
    m_device void meh_root_params_(MehParams* p, real_t* rm1, const real_t a[2], const real_t y[2]) {
        real_t r;
        meh_root(a, y, &r, rm1);
        p->r       = r;
        p->use_log = (r > 0.5);
        p->logr    = p->use_log ? log1p(*rm1) : 0.0;
        p->n_cut   = decay_cut(log(fabs(r)));
    }

    // ----------------------------------------------------------------------------------------
    // Mehrstellen 4
    // ----------------------------------------------------------------------------------------
    /**
     * @brief evaluate the coefficients a_j(y) for the meh4 stencil
     */
    m_device void meh4_left_coeffs(real_t a[2], const real_t y[2]) {
        a[0] = -8./3.*(y[0]*y[1] - y[0] - y[1]) + 2.0;
        a[1] = 2./3.*(y[0] + y[1]) - 1.0;
    }

    /**
     * @brief evaluate the coefficients b_j(y) for the meh4 stencil
     */
    m_device void meh4_right_coeffs(real_t b[2], const real_t y[2]) {
        b[0] = 5./6. - 1./3.*(y[0] + y[1]);
        b[1] = +1./12.;
    }

    /**
     * @brief the n-independent data of @ref meh4_left
     */
    m_device MehParams meh4_left_params(const real_t k1, const real_t k2) {
        MehParams p = {};
        p.branch = kZero;
        p.n_cut  = INT_MAX;  // the -|n|/2 term does not decay
        if (k1 == 0 && k2 == 0) {
            return p;
        }
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        real_t a[2];
        meh4_left_coeffs(a, y);
        real_t rm1;
        meh_root_params_(&p, &rm1, a, y);
        if (a[1] != 0.0) {
            p.branch = kGeneric;
            p.num    = p.r / a[1];
            p.den    = rm1 * (p.r + 1.0);
        } else {
            p.branch = kDegenerate;
            p.n_cut  = 1;  // only n = 0 is non-zero
            p.a0     = a[0];
        }
        return p;
    }

    /**
     * @brief the n-dependent part of @ref meh4_left
     */
    m_device real_t meh4_left_eval(const MehParams& p, const int n) {
        if (p.branch == kZero) {
            return -0.5 * abs(n);
        }
        const real_t base = p.num * meh_rpow_(p, abs(n)) / p.den;
        return (p.branch == kGeneric) ? base : (n == 0) / p.a0;
    }

    /**
     * @brief evaluate the LGF G(n, k1, k2) for the left hand side of the meh4 stencil
     */
    m_device real_t meh4_left(const int n, const real_t k1, const real_t k2) {
        return meh4_left_eval(meh4_left_params(k1, k2), n);
    };
    
    /**
     * @brief the n-independent data of @ref meh4_full
     */
    m_device MehParams meh4_full_params(const real_t k1, const real_t k2) {
        MehParams p = {};
        p.branch = kZero;
        p.n_cut  = INT_MAX;  // the -|n|/2 term does not decay
        if (k1 == 0 && k2 == 0) {
            return p;
        }; 
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        real_t a[2], b[2];
        meh4_left_coeffs(a, y);
        meh4_right_coeffs(b, y);
        real_t rm1;
        meh_root_params_(&p, &rm1, a, y);
        const real_t r = p.r;
        p.num = b[1] + b[0]*r + b[1]*r*r;  // p(r)
        if (a[1] != 0.0) {
            const real_t r2m1 = rm1 * (r + 1.0);  // r^2 - 1
            p.branch = kGeneric;
            p.den    = a[1]*r2m1;
            p.v0     = (r/a[1])*(b[0] + 2.*b[1]*r)/r2m1;
        } else {
            p.branch = kDegenerate;
            p.n_cut  = 1;  // only n = 0 and n = 1 are non-zero
            p.b0a0   = b[0] / a[0];
            p.b1a0   = b[1] / a[0];
        }
        return p;
    }

    /**
     * @brief the n-dependent part of @ref meh4_full
     */
    m_device real_t meh4_full_eval(const MehParams& p, const int n) {
        if (p.branch == kZero) {
            return -0.5 * abs(n) - (n == 0) * 1./12.;
        }
        if (p.branch == kGeneric) {
            const real_t base = p.num * meh_rpow_(p, abs(n)) / p.den;
            return (n == 0) ? p.v0 : base;
        } else {
            return (n == 0) * p.b0a0 + (n == 1) * p.b1a0;
        }
    }

    /**
     * @brief evaluate the LGF G(n, k1, k2) for the full meh4 stencil
     */
    m_device real_t meh4_full(const int n, const real_t k1, const real_t k2) {
        return meh4_full_eval(meh4_full_params(k1, k2), n);
    }

    // ----------------------------------------------------------------------------------------
    // Mehrstellen 4
    // ----------------------------------------------------------------------------------------
    /**
     * @brief evaluate the coefficients a_j(y) for the meh6 stencil
     */
    m_device void meh6_left_coeffs(real_t a[2], const real_t y[2]) {
        a[0] = -8./5.*y[0]*y[1] + 8./3.*(y[0] + y[1]) + 2.0;
        a[1] = -8./15.*y[0]*y[1] + 2./3.*(y[0] + y[1]) - 1.0;
    }

    /**
     * @brief evaluate the coefficients b_j(y) for the meh6 stencil
     */
    m_device void meh6_right_coeffs(real_t b[3], const real_t y[2]) {
        b[0] = -1./15.*(pow(y[0], 2) + pow(y[1], 2)) + 8./45.*y[0]*y[1] - 11./45.*(y[0] + y[1]) + 97./120.;
        b[1] = -2./45.*(y[0] + y[1]) + 1./10.;
        b[2] = -1./240.;
    }

    /**
     * @brief the n-independent data of @ref meh6_left
     */
    m_device MehParams meh6_left_params(const real_t k1, const real_t k2) {
        MehParams p = {};
        p.branch = kZero;
        p.n_cut  = INT_MAX;  // the -|n|/2 term does not decay
        if (k1 == 0 && k2 == 0) {
            return p;
        }; 
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        real_t a[2];
        meh6_left_coeffs(a, y);
        real_t rm1;
        meh_root_params_(&p, &rm1, a, y);
        p.branch = kGeneric;
        p.num    = p.r / a[1];
        p.den    = rm1 * (p.r + 1.0);
        return p;
    }

    /**
     * @brief the n-dependent part of @ref meh6_left
     */
    m_device real_t meh6_left_eval(const MehParams& p, const int n) {
        if (p.branch == kZero) {
            return -0.5 * abs(n);
        };
        return p.num * meh_rpow_(p, abs(n)) / p.den;
    }

    /**
     * @brief evaluate the LGF G(n, k1, k2) for the left hand side of the meh6 stencil
     */
    m_device real_t meh6_left(const int n, const real_t k1, const real_t k2) {
        return meh6_left_eval(meh6_left_params(k1, k2), n);
    };

    /**
     * @brief the n-independent data of @ref meh6_full
     */
    m_device MehParams meh6_full_params(const real_t k1, const real_t k2) {
        MehParams p = {};
        p.branch = kZero;
        p.n_cut  = INT_MAX;  // the -|n|/2 term does not decay
        if (k1 == 0 && k2 == 0) {
            return p;
        }; 
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        real_t a[2], b[3];
        meh6_left_coeffs(a, y);
        meh6_right_coeffs(b, y);
        real_t rm1;
        meh_root_params_(&p, &rm1, a, y);
        const real_t r = p.r;
        p.branch = kGeneric;
        p.num    = b[0]*pow(r, 2) + b[1]*(r + pow(r, 3)) + b[2]*(1.0 + pow(r, 4));
        p.den    = a[1] * rm1 * (r + 1.0);
        p.c1     = b[2] / a[1];
        p.c0num  = b[1]*a[1] - b[2]*a[0];
        p.c0den  = a[1]*a[1];
        // the decaying term carries r^(|n| - 1), i.e. one power less than @ref decay_cut assumes
        p.n_cut  = (p.n_cut < INT_MAX) ? (p.n_cut + 1) : INT_MAX;
        return p;
    }

    /**
     * @brief the n-dependent part of @ref meh6_full
     */
    m_device real_t meh6_full_eval(const MehParams& p, const int n) {
        if (p.branch == kZero) {
            return -0.5 * abs(n) - (abs(n) == 1) * (-1./240.) - (n == 0) * (11./120.);
        };
        const real_t base = p.num * meh_rpow_(p, abs(n) - 1) / p.den;
        return base + (abs(n) == 1) * p.c1 + (n == 0) * p.c0num / p.c0den;
    }

    /**
     * @brief evaluate the LGF G(n, k1, k2) for the full meh6 stencil
     */
    m_device real_t meh6_full(const int n, const real_t k1, const real_t k2) {
        return meh6_full_eval(meh6_full_params(k1, k2), n);
    }
};

#endif