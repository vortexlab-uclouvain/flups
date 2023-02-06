#ifndef LGF_1UNB2SPE_HPP
#define LGF_1UNB2SPE_HPP

#include <cmath>
#include <complex>
#include <array>
#include <tuple>

namespace LGFOneUnbounded {

    // ----------------------------------------------------------------------------------------
    // Types
    // ----------------------------------------------------------------------------------------

    using real_t = double;                  //!< Note: floating point literals in this file are double precision.
    using complex_t = std::complex<real_t>; //!< complex precision matches real precision 

    // ----------------------------------------------------------------------------------------
    // Utility functions
    // ----------------------------------------------------------------------------------------
    static const complex_t im = complex_t(0.0, 1.0); //!< constant for i = sqrt(-1)

    /**
     * @brief a square root with complex return type that accepts real or complex arguments
     */
    template <typename T> inline complex_t complex_sqrt(T z) { 
        return sqrt(static_cast<complex_t>(z)); 
    }

    /**
     * @brief evaluate a polynomial with real coefficients via the Horner scheme
     * 
     * @param a an array of length n containing coefficients a_{0} to a_{n-1} 
     * @param n the number of coefficients
     * @param x the evaluation point
     */
    template<typename T> inline T horner(const real_t* a, int n, T x) {
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
    constexpr int lgf_log_series_terms_ = 12;
    inline real_t log_series(real_t x_minus_one) {
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
    constexpr int lgf_value_minus_inverse_series_terms_ = 12;
    inline real_t value_minus_inverse_series(real_t x_minus_one) {
        real_t out = 0.0;
        for (int i = 0; i < lgf_value_minus_inverse_series_terms_; i++) {
            out += 1.0;
            out *= -x_minus_one;
        }
        return -out + x_minus_one;
    }

    /**
     * @brief return the oscillatory contribution to the LGF from a conjugate root pair
     * 
     * @param r one of the complex roots
     * @param pp the complex factor p'(r, c)
     * @param w the width of the stencil
     * @param n the evaluation index
     */
    const real_t oscillatory_form(complex_t r, complex_t pp, int w, int n) {
        const complex_t K = complex_t(0.0, -2.0)*pow(r, w - 1)/pp;
        return abs(K) * pow(abs(r), abs(n)) * sin(arg(r)*abs(n) + arg(K));
    }

    // ----------------------------------------------------------------------------------------
    // Reusable constants
    // ----------------------------------------------------------------------------------------

    // root of q(lambda) nearest 1
    constexpr int lgf_small_c_terms_ = 8;
    static const real_t lgf4_small_c_expansion_[8] = {0, 0.5, 1./24., 1./144., 5./3456., 7./20736., 7./82944., 11./497664.};
    static const real_t lgf6_small_c_expansion_[8] = {0, 0.5, 1./24., 1./720., -1./1152., -149./518400., -259./6220800., 163./62208000.};

    // p'(z, c), template to allow real_t or complex
    template <typename T> inline T lgf4_ppzc_(T z, real_t c) { return 4./3. + z*((-5. - 2.*c) + z*(4. - z/3.)); }
    template <typename T> inline T lgf6_ppzc_(T z, real_t c) { return -3./20. + z*(3. + z*((-49./6. - 3.*c) + z*(6. + z*(-3./4. + z/15.)))); }

    // LGF4 repeated root residues
    static const real_t lgf4_c_star_ = 3.0;
    const real_t sq15 = sqrt(15.);
    const real_t lgf4_repeat_g0_[2]  = {-4./(5.*sq15), -1./5.};
    const real_t lgf4_repeat_g1_[4]  = {14./(75.*sq15), 2./25., 4./(25.*sq15), 1./150.};
    const real_t lgf4_repeat_g2_[6]  = {-901./(18750.*sq15), -74./3125., -253./(3750.*sq15), -23./3750., -1./(250.*sq15), -1./15000.};
    const real_t lgf4_repeat_g3_[8]  = {37313./(2812500.*sq15), 11267./1640625., 3167./(140625.*sq15), 1./375., 31./(11250.*sq15), 31./281250., 1./(28125.*sq15), 1./3150000.};
    const real_t lgf4_repeat_g4_[10] = {-3885631./(1012500000.*sq15), -1194829./590625000., -1118351./(157500000.*sq15), -12571./13125000., -13751./(11250000.*sq15), -127./1875000., -247./(6750000.*sq15), -13./15750000., -1./(6300000.*sq15), -1./1134000000.};
    
    // ----------------------------------------------------------------------------------------
    // Symbols
    // ----------------------------------------------------------------------------------------
    inline real_t lgf2_symbol(real_t k) { return 4.0 * pow(sin(k/2.), 2); }
    inline real_t lgf4_symbol(real_t k) { return 4.0 * (4./3.*pow(sin(k/2.), 2) - 1./12.*pow(sin(k), 2)); }
    inline real_t lgf6_symbol(real_t k) { return 4.0 * (3./2.*pow(sin(k/2.), 2) - 3./20.*pow(sin(k), 2) + 1./90.*pow(sin(1.5*k), 2)); }

    // ----------------------------------------------------------------------------------------
    // LGF2
    // ----------------------------------------------------------------------------------------
    const real_t lgf2_small_cutoff_ = 1.0e-3;
    real_t lgf2(int n, real_t c) {
        if (c < 0) {
            printf("[ExpandLGF] Argument error: c must be nonnegative");
            abort();
        } else if (c == 0.0) {
            return -0.5 * abs(n);
        } else if (c < lgf2_small_cutoff_) {
            const real_t rm1 = c/2. - sqrt(c*c/4. + c);
            return -exp(abs(n) * log_series(rm1)) / value_minus_inverse_series(rm1);
        } else {
            const real_t r = 1 + c/2. - sqrt(c*c/4. + c);
            return -pow(r, abs(n)) / (r - 1./r); 
        };
    }

    // ----------------------------------------------------------------------------------------
    // LGF4
    // ----------------------------------------------------------------------------------------

    /** 
     * @brief evaluate lgf4(n, c) for c = 0 
     */
    inline real_t lgf4_c_zero_(int n) {
        const real_t r = 7.179676972449123e-2;
        const real_t K = 7.216878364870408e-2;
        return -0.5*abs(n) + K*(1.0 - pow(r, abs(n)));
    }    

    /** 
     * @brief evaluate lgf4(n, c) for 0 < c < lgf4_small_c_cutoff_
     */
    inline real_t lgf4_c_small_(int n, real_t c) {
        // larger real root
        const real_t lambdam1 = horner(lgf4_small_c_expansion_, lgf_small_c_terms_, c);
        const real_t rm1 = lambdam1 - sqrt(lambdam1*lambdam1 + 2*lambdam1);
        const real_t r1 = rm1 + 1.0;
        // smaller real root
        const real_t lambda2 = 4.0 + sqrt(9.0 - 3.0*c);
        const real_t r2 = lambda2 - sqrt(lambda2*lambda2 - 1.0);
        // assemble LGF value
        const real_t r1_minus_r1_inv = value_minus_inverse_series(rm1);
        const real_t logr1 = log_series(rm1);
        const real_t pp1 = -1./12. * r1_minus_r1_inv * (r1 - r2) * (r1 - 1.0/r2);
        const real_t pp2 = lgf4_ppzc_(r2, c);
        return -exp((abs(n) + 1)*logr1) / pp1 - pow(r2, abs(n) + 1) / pp2;
    }

    /**
     * @brief evaluate lgf(n, c) for lgf4_small_cutoff_ < c < (lgf4_c_star_ - lgf4_repeat_cutoff_)
     */
    inline real_t lgf4_c_real_(int n, real_t c) {
        const real_t lambda1 = 4.0 + sqrt(9.0 - 3.0*c);
        const real_t lambda2 = 4.0 - sqrt(9.0 - 3.0*c);
        const real_t r1 = lambda1 - sqrt(lambda1*lambda1 - 1.0);
        const real_t r2 = lambda2 - sqrt(lambda2*lambda2 - 1.0);
        const real_t f1 = 12*pow(r1, abs(n) + 1)/((r1 - 1.0/r1)*(r1 - 1.0/r2));
        const real_t f2 = 12*pow(r2, abs(n) + 1)/((r2 - 1.0/r1)*(r2 - 1.0/r2));
        return (f1 - f2) / (r1 - r2);
    }

    /**
     * @brief evaluate lgf(n, c) for abs(c - lgf4_c_star_) < lgf4_repeat_cutoff_
     */
    inline real_t lgf4_c_repeat_(int n, real_t c) {
        const real_t r = 4.0 - sqrt(15.0);
        const real_t nr = static_cast<real_t>(abs(n));
        const real_t lgf4_g[5] = {
            horner(lgf4_repeat_g0_, 2,  nr),
            horner(lgf4_repeat_g1_, 4,  nr),
            horner(lgf4_repeat_g2_, 6,  nr),
            horner(lgf4_repeat_g3_, 8,  nr),
            horner(lgf4_repeat_g4_, 10, nr)
        };
        const real_t dc = c - lgf4_c_star_;
        return -pow(r, abs(n)) * horner(lgf4_g, 5, dc);
    }

    /**
     * @brief evaluate lgf4(n, c) for (lgf4_c_star_ + lgf4_repeat_cutoff_) < c
     */
    inline real_t lgf4_c_imag_(int n, real_t c) {
        const complex_t lambda = 4.0 - complex_sqrt(9.0 - 3.0*c);
        const complex_t r = lambda - sqrt(lambda*lambda - 1.0);
        const complex_t pp = -im/6.*imag(r)*(r - 1./r)*(r - 1./conj(r)); // lgf4_ppzc_(r, c);
        return oscillatory_form(r, pp, 2, n);
    }

    constexpr real_t lgf4_small_cutoff_ = 1e-3;
    constexpr real_t lgf4_repeat_cutoff_ = 1e-5;

    real_t lgf4(int n, real_t c) {
        if (c < 0) {
            printf("[ExpandLGF] Argument error: c must be nonnegative");
            abort();
        } else if (c == 0.0) {
            return lgf4_c_zero_(n);
        } else if (c < lgf4_small_cutoff_) {
            return lgf4_c_small_(n, c);
        } else if (c < lgf4_c_star_ - lgf4_repeat_cutoff_) {
            return lgf4_c_real_(n, c);
        } else if (c < lgf4_c_star_ + lgf4_repeat_cutoff_) {
            return lgf4_c_repeat_(n, c);
        } else {
            return lgf4_c_imag_(n, c);
        };
    }

    // ----------------------------------------------------------------------------------------
    // LGF6
    // ----------------------------------------------------------------------------------------

    /**
     * @brief return the roots of p(z, c) for the LGF6 stencil
     * 
     * Because r2 and r3 are always a conjugate pair, r3 is omitted from the return.
     */
    std::tuple<real_t, complex_t> lgf6_roots_(real_t c) {
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

    /**
     * @brief evaluate lgf6(n, c) for c = 0
     */
    inline real_t lgf6_c_zero_(int n) {
        const real_t a = +1.035017642406219e-1;
        const real_t p = +9.54295097909627e-2; // TODO
        const real_t t = -9.958762687216868e-1;
        const real_t f = -2.168575565241744e+0;
        return -0.5*abs(n) + a*pow(p, abs(n))*sin(abs(n)*t + f) - a*sin(f);
    }

    /**
     * @brief evaluate lgf6(n, c) for 0 < c < lgf6_small_cutoff_
     */
    inline real_t lgf6_c_small_(int n, real_t c) {
        // larger real root
        const real_t lambdam1 = horner(lgf6_small_c_expansion_, lgf_small_c_terms_, c);
        const real_t rm1 = lambdam1 - sqrt(lambdam1*lambdam1 + 2*lambdam1);
        const real_t r1 = rm1 + 1.0;
        // conjugate root pair
        const auto roots = lgf6_roots_(c);
        const complex_t r2 = std::get<1>(roots);
        const complex_t pp2 = lgf6_ppzc_(r2, c);
        // assemble LGF value
        const real_t r1_minus_r1_inv = value_minus_inverse_series(rm1);
        const real_t logr1 = log_series(rm1);
        const real_t pp1 = 1./90. * r1_minus_r1_inv * norm((r1 - r2) * (r1 - 1.0/r2));
        return -exp((abs(n) + 2)*logr1) / pp1 + oscillatory_form(r2, pp2, 3, n);
    }

    /**
     * @brief evaluate lgf6(n, c) for lgf6_small_cutoff_ < c
     */
    inline real_t lgf6_c_generic_(int n, real_t c) {
        const auto roots = lgf6_roots_(c);
        const real_t r1 = std::get<0>(roots);
        const real_t pp1 = lgf6_ppzc_(r1, c);
        const complex_t r2 = std::get<1>(roots);
        const complex_t pp2 = lgf6_ppzc_(r2, c);
        return -pow(r1, abs(n) + 2)/pp1 + oscillatory_form(r2, pp2, 3, n);
    }

    constexpr real_t lgf6_small_cutoff_ = 1e-3;

    real_t lgf6(int n, real_t c) {
        if (c < 0) {
            printf("[ExpandLGF] Argument error: c must be nonnegative");
            abort();
        } else if (c == 0.0) {
            return lgf6_c_zero_(n);
        } else if (c < lgf6_small_cutoff_) {
            return lgf6_c_small_(n, c);
        } else {
            return lgf6_c_generic_(n, c);
        };
    }
}; // namespace LGFOneUnbounded

#endif