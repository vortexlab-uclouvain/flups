#ifndef LGF_1UNB2SPE_HPP
#define LGF_1UNB2SPE_HPP

#include <cmath>
#include <complex>
#include <array>
#include <tuple>

#include "defines.hpp"

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
    // large lambda expansion of r(lambda)
    constexpr int large_lambda_terms_ = 10;
    static const real_t large_lambda_expansion_[] = {0, 1./2., 0, 1./8., 0, 1./16., 0, 5./128., 0, 7./256.};

    // root of q(lambda) nearest 1
    constexpr int lgf_small_c_terms_ = 8;
    static const real_t lgf4_small_c_expansion_[8] = {0, 0.5, 1./24., 1./144., 5./3456., 7./20736., 7./82944., 11./497664.};
    static const real_t lgf6_small_c_expansion_[8] = {0, 0.5, 1./24., 1./720., -1./1152., -149./518400., -259./6220800., 163./62208000.};
    static const real_t lgf8_small_c_expansion_[8] = {0, 0.5, 1./24., 1./720., 1./40320., 577./3628800., 389./6220800., 34987./3048192000.};

    // p'(z, c), template to allow real_t or complex
    template <typename T> inline T lgf4_ppzc_(T z, real_t c) { return 4./3. + z*((-5. - 2.*c) + z*(4. - z/3.)); }
    template <typename T> inline T lgf6_ppzc_(T z, real_t c) { return -3./20. + z*(3. + z*((-49./6. - 3.*c) + z*(6. + z*(-3./4. + z/15.)))); }
    template <typename T> inline T lgf8_ppzc_(T z, real_t c) { return 8./315. + z*(-2./5. + z*(24./5. + z*(-205./18. - 4.*c + z*(8. + z*(-6./5. + z*(8./45. - z/70.)))))); } 

    // LGF4 repeated root residues
    static const real_t lgf4_c_star_ = 3.0;
    const real_t sq15 = sqrt(15.);
    const real_t lgf4_repeat_g0_[2]  = {-4./(5.*sq15), -1./5.};
    const real_t lgf4_repeat_g1_[4]  = {14./(75.*sq15), 2./25., 4./(25.*sq15), 1./150.};
    const real_t lgf4_repeat_g2_[6]  = {-901./(18750.*sq15), -74./3125., -253./(3750.*sq15), -23./3750., -1./(250.*sq15), -1./15000.};
    const real_t lgf4_repeat_g3_[8]  = {37313./(2812500.*sq15), 11267./1640625., 3167./(140625.*sq15), 1./375., 31./(11250.*sq15), 31./281250., 1./(28125.*sq15), 1./3150000.};
    const real_t lgf4_repeat_g4_[10] = {-3885631./(1012500000.*sq15), -1194829./590625000., -1118351./(157500000.*sq15), -12571./13125000., -13751./(11250000.*sq15), -127./1875000., -247./(6750000.*sq15), -13./15750000., -1./(6300000.*sq15), -1./1134000000.};
    
    // LGF8 repeated root residues
    static const real_t lgf8_c_star_ = 3.204471924659902;

    const real_t lgf8_a_0_[2]  = {-0.1932304295173025, -9.940958723693284e-2};
    const real_t lgf8_a_1_[4]  = {5.528645180523051e-2, 6.381749748906844e-2, 1.920895723993029e-2, 1.647044339102893e-3};
    const real_t lgf8_a_2_[6]  = {-1.246195963786745e-2, -2.200316066032414e-2, -1.156690496243602e-2, -2.513956535483380e-3, -2.386943138091703e-4, -8.186599895557273e-6};
    const real_t lgf8_a_3_[8]  = {2.586595339831205e-3, 5.946351585490818e-3, 4.392324394119772e-3, 1.460419122054276e-3, 2.498408821419240e-4, 2.281049273902400e-5, 1.054600142736557e-6, 1.937682182122062e-8};
    const real_t lgf8_a_4_[10] = {-6.143387386724980e-4, -1.461770070530813e-3, -1.338631620906462e-3, -5.893614275050890e-4, -1.419718637870522e-4, -2.006233681173800e-5, -1.705558618042543e-6, -8.568226736629197e-8, -2.340119751996688e-9, -2.675335915571021e-11};

    const real_t lgf8_u_0_[1] = {0.0015031218476240784};
    const real_t lgf8_u_1_[2] = {-0.013198947960911374, -0.00309497626043345};
    const real_t lgf8_u_2_[3] = {0.002223792419320104, 0.0003211726758356002, -3.489943855291188e-6};
    const real_t lgf8_u_3_[4] = {8.218756082933073e-5, 0.00012939734109457704, 2.076196709989474e-5, 7.970739601754506e-7};
    const real_t lgf8_u_4_[5] = {-0.00011284539233578474, -4.455888056601662e-5, -4.596415525190604e-6, -1.2515613302441402e-7, 7.499072573978986e-10};

    const real_t lgf8_v_0_[1] = {-0.07869060869097218};
    const real_t lgf8_v_1_[2] = {0.005694070346472886, -0.00011828157312623746};
    const real_t lgf8_v_2_[3] = {0.001748394270280815, 0.0007853940979982054, 6.08419435502199e-5};
    const real_t lgf8_v_3_[4] = {-0.000574565978306078, -0.000152486888729598, -8.030910109571855e-6, 6.101311014540467e-8};
    const real_t lgf8_v_4_[5] = {4.706634682161246e-5, -8.921284296708267e-6, -4.142946014756904e-6, -3.42884013573887e-7, -7.828815619884859e-9};
    
    // ----------------------------------------------------------------------------------------
    // Symbols
    // ----------------------------------------------------------------------------------------
    inline real_t lgf2_symbol(real_t k) { return 4.0 * pow(sin(k/2.), 2); }
    inline real_t lgf4_symbol(real_t k) { return 4.0 * (4./3.*pow(sin(k/2.), 2) - 1./12.*pow(sin(k), 2)); }
    inline real_t lgf6_symbol(real_t k) { return 4.0 * (3./2.*pow(sin(k/2.), 2) - 3./20.*pow(sin(k), 2) + 1./90.*pow(sin(1.5*k), 2)); }
    inline real_t lgf8_symbol(real_t k) { return 4.0 * (8./5.*pow(sin(k/2.), 2) - 1./5. *pow(sin(k), 2) + 8./315.*pow(sin(1.5*k), 2) - 1./560.*pow(sin(2.*k), 2)); }

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


    // ----------------------------------------------------------------------------------------
    // LGF8
    // ----------------------------------------------------------------------------------------

    /**
     * @brief return the roots of p(z, c) for LGF8. 
     * 
     * Because r3 and r4 are always a conjugate pair, r4 is omitted
     */
    std::array<complex_t, 3> lgf8_roots_(real_t c) {
        const real_t tmp1 = 1968941520. + c*(-879221700. + c*(223915104. - c*44089920));
        const complex_t xi = pow(1520225. - 273420.*c + 35.*complex_sqrt(tmp1), 1./3.);
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

    /**
     * @brief evaluate lgf8(n, c) for c = 0
     */
    inline real_t lgf8_c_zero_(int n) {
        const real_t rho1 = +1.2185576475724652e-01;
        const real_t th1  = +1.4667413886166037e+00;
        const real_t a1   = +8.2123351081853858e-02;
        const real_t phi1 = -6.5666289805823130e-01;
        const real_t r2   = +9.6191214171019512e-02;
        const real_t a2   = +4.0692417408429382e-02;
        return -0.5*abs(n) + a1*pow(rho1, abs(n))*sin(abs(n)*th1 + phi1) - a1*sin(phi1) + a2*(1.0 - pow(r2, abs(n)));
    }

    /**
     * @brief evaluate lgf8(n, c) for 0 < c < lgf8_small_cutoff_
     */
    inline real_t lgf8_c_small_(int n, real_t c) {
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
        // assemble LGF value
        const real_t r1_minus_r1_inv = value_minus_inverse_series(rm1);
        const real_t logr1 = log_series(rm1);
        const real_t pp1 = -1./560. * r1_minus_r1_inv * (r1 - r2)*(r1 - 1/r2)*norm((r1 - r3)*(r1 - 1.0/r3));
        return -exp((abs(n) + 3)*logr1)/pp1 - pow(r2, abs(n) + 3)/pp2 + oscillatory_form(r3, pp3, 4, n);
    }

    /**
     * @brief evaluate lgf8(n, c) for lgf8_small_cutoff_ < c < (lgf8_c_star_ - lgf8_repeat_cutoff_) 
     */
    inline real_t lgf8_c_real_(int n, real_t c) {
        const auto roots = lgf8_roots_(c);
        const real_t r1 = real(roots[0]);
        const real_t r2 = real(roots[1]);
        const complex_t r3 = roots[2];
        const complex_t pp3 = lgf8_ppzc_(r3, c);
        const real_t f1 = 560. * pow(r1, abs(n) + 3) / ((r1 - 1./r1)*(r1 - 1./r2)*norm((r1 - r3)*(r1 - 1./r3)));
        const real_t f2 = 560. * pow(r2, abs(n) + 3) / ((r2 - 1./r1)*(r2 - 1./r2)*norm((r2 - r3)*(r2 - 1./r3)));
        return (f2 - f1) / (r2 - r1) + oscillatory_form(r3, pp3, 4, n);
    }

    /**
     * @brief evaluate lgf8(n, c) for abs(c - lgf8_c_star_) < lgf8_repeat_cutoff_
     */
    inline real_t lgf8_c_repeat_(int n, real_t c) {
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
        const real_t dc = c - lgf8_c_star_;
        return -horner(lgf8_g, 5, dc);
    }

    /**
     * @brief evaluate lgf8(n, c) for (lgf8_c_star_ + lgf8_repeat_cutoff_) < c 
     */
    inline real_t lgf8_c_imag_(int n, real_t c) {
        const auto roots = lgf8_roots_(c);
        const complex_t r1 = roots[1];
        const complex_t r3 = roots[2];
        const complex_t cr1 = conj(r1);
        const complex_t cr3 = conj(r3);
        const complex_t pp1 = -im/280.*imag(r1)*(r1 - 1./r1)*(r1 - 1./cr1)*(r1 - r3)*(r1 - 1./r3)*(r1 - cr3)*(r1 - 1./cr3); // lgf8_ppzc_(r1, c);
        const complex_t pp3 = lgf8_ppzc_(r3, c);
        return oscillatory_form(r1, pp1, 4, n) + oscillatory_form(r3, pp3, 4, n); 
    }

    constexpr real_t lgf8_small_cutoff_ = 1e-3;
    constexpr real_t lgf8_repeat_cutoff_ = 1e-5;

    real_t lgf8(int n, real_t c) {
        if (c < 0) {
            printf("[ExpandLGF] Argument error: c must be nonnegative");
            abort();
        } else if (c == 0.0) {
            return lgf8_c_zero_(n);
        } else if (c < lgf4_small_cutoff_) {
            return lgf8_c_small_(n, c);
        } else if (c < lgf8_c_star_ - lgf8_repeat_cutoff_) {
            return lgf8_c_real_(n, c);
        } else if (c < lgf8_c_star_ + lgf8_repeat_cutoff_) {
            return lgf8_c_repeat_(n, c);
        } else {
            return lgf8_c_imag_(n, c);
        };
    }

    // ----------------------------------------------------------------------------------------
    // Mehrstellen Stencils 4 and 6
    // ----------------------------------------------------------------------------------------
    static const real_t meh_large_lambda_tol_ = 1.0e2;
    real_t mehrstellen_r(real_t lambda) {
        real_t r = 0.0;
        if (std::isinf(lambda)) { // infinite lambda
            r = 0.0;
        } else if (fabs(lambda) > meh_large_lambda_tol_) { // large lambda
            r = horner(large_lambda_expansion_, large_lambda_terms_, 1.0 / lambda);
        } else { // moderate lambda
            r = real(lambda - complex_sqrt(lambda - 1.)*complex_sqrt(lambda + 1.));
        }
        return r;
    }

    real_t meh4_left(const int n, const real_t k1, const real_t k2) {
        if (k1 == 0 && k2 == 0) {
            return -0.5 * abs(n);
        }
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        const real_t q0 = +8./3.*(y[0]*y[1] - y[0] - y[1]) - 2.0;
        const real_t q1 = -4./3.*(y[0] + y[1]) + 2.0;
        const real_t lambda = -q0 / q1;
        const real_t r = mehrstellen_r(lambda);
        const real_t base = 2. * (r / q1) * pow(r, abs(n)) / (1.0 - r*r);
        return std::isfinite(lambda) ? base : -(n == 0) / q0; 
    };

    real_t meh6_left(const int n, const real_t k1, const real_t k2) {
        if (k1 == 0 && k2 == 0) {
            return -0.5 * abs(n);
        }; 
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        const real_t q0 = 8./5.*y[0]*y[1] - 8./3.*(y[0] + y[1]) - 2.0;
        const real_t q1 = 16./15.*y[0]*y[1] - 4./3.*(y[0] + y[1]) + 2.0;
        const real_t lambda = -q0 / q1;
        const real_t r = mehrstellen_r(lambda);
        const real_t base = 2. * (r / q1) * pow(r, abs(n)) / (1.0 - r*r);
        return std::isfinite(lambda) ? base : -(n == 0) / q0; 
    };
    
    real_t meh4_full(const int n, const real_t k1, const real_t k2) {
        if (k1 == 0 && k2 == 0) {
            return -0.5 * abs(n) - (n == 0) * 1./12.;
        }; 
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        const real_t pl1 = +8./3.*(y[0]*y[1] - y[0] - y[1]) - 2.0;
        const real_t pl0 = -2./3.*(y[0] + y[1]) + 1.0;
        const real_t pr1 = -1./3.*(y[0] + y[1]) + 5./6.;
        const real_t pr0 = 1./12.;
        const real_t lambda = -pl1 / (2*pl0);
        const real_t r = mehrstellen_r(lambda);
        if (std::isfinite(lambda)) {
            const real_t base = pow(r, abs(n)) * (pr1*r + pr0*(1 + r*r)) / (pl0*(r*r - 1));
            return (n == 0) ? -(r/pl0)*(pr1 + 2.*pr0*r)/(r*r - 1.) : -base;
        } else {
            return (n == 0) * (-pr1 / pl1) + (n == 1) * (-pr0 / pl1); 
        }
    }

    real_t meh6_full(const int n, const real_t k1, const real_t k2) {
        if (k1 == 0 && k2 == 0) {
            return -0.5 * abs(n) - (abs(n) == 1) * (-1./240.) - (n == 0) * (11./120.);
        }; 
        // deal with the base portion
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        const real_t ql0 = 8./5.*y[0]*y[1] - 8./3.*(y[0] + y[1]) - 2.0;
        const real_t ql1 = 16./15.*y[0]*y[1] - 4./3.*(y[0] + y[1]) + 2.0;
        const real_t lambda = -ql0 / ql1;
        const real_t qr0 = -1./15.*(pow(y[0], 2) + pow(y[1], 2)) + 8./45.*y[0]*y[1] - 11./45.*(y[0] + y[1]) + 49./60.;
        const real_t qr1 = -4./45.*(y[0] + y[1]) + 1./5.;
        const real_t qr2 = -1./60.;
        const real_t qr = std::isinf(lambda) ? 0.0 : qr0 + lambda*(qr1 + lambda*qr2);
        const real_t r = mehrstellen_r(lambda);
        const real_t base = 2. * qr * pow(r, abs(n) + 1) / (ql1 * (1.0 - r*r));
        // corrections near the origin
        const real_t pl0 = 8./5.*y[0]*y[1] - 8./3.*(y[0] + y[1]) - 2.0;
        const real_t pl1 = 8./15.*y[0]*y[1] - 2./3.*(y[0] + y[1]) + 1.0;
        const real_t pr1 = -2./45.*(y[0] + y[1]) + 1./10.;
        const real_t pr2 = -1./240.;
        return base - (abs(n) == 1) * (pr2 / pl1) - (n == 0) * (pr1*pl1 - pr2*pl0) / (pl1*pl1);
    }

    void meh4_left_coeffs(real_t coeffs[2], const real_t k1, const real_t k2) {
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        coeffs[0] = +8./3.*(y[0]*y[1] - y[0] - y[1]) - 2.0;
        coeffs[1] = -2./3.*(y[0] + y[1]) + 1.0;
    }

    void meh4_right_coeffs(real_t coeffs[2], const real_t k1, const real_t k2) {
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        coeffs[0] = -1./3.*(y[0] + y[1]) + 5./6.;
        coeffs[1] = 1./12.;
    }

    void meh6_left_coeffs(real_t coeffs[2], const real_t k1, const real_t k2) {
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        coeffs[0] = 8./5.*y[0]*y[1] - 8./3.*(y[0] + y[1]) - 2.0;
        coeffs[1] = 8./15.*y[0]*y[1] - 2./3.*(y[0] + y[1]) + 1.0;
    }

    void meh6_right_coeffs(real_t coeffs[3], const real_t k1, const real_t k2) {
        const real_t y[2] = {pow(sin(k1/2), 2), pow(sin(k2/2), 2)};
        coeffs[0] = -1./15.*(pow(y[0], 2) + pow(y[1], 2)) + 8./45.*y[0]*y[1] - 11./45.*(y[0] + y[1]) + 97./120.;
        coeffs[1] = -2./45.*(y[0] + y[1]) + 1./10.;
        coeffs[2] = -1./240.;
    }

}; // namespace LGFOneUnbounded

#endif