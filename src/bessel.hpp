/**
 * @file bessel.hpp
 * @author Denis-Gabriel Caprace
 * @brief implementation of the modified bessel functions (with nu=0).
 * @date 2019-09-05
 * 
 * References:
 * - Abramowitz and Stegun, "Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables", 1964; §9.4 "Bessel functions"
 * - Press et al., "Numerical Recipes", 3rd edition, Cambridge University Press, 2007; §6.5.1 "Modified Bessel Functions of Integer Order", pp. 279
 */

#include <cmath>
#include <iostream>


#define APPROX

#ifdef APPROX
static inline double besseli0(double x) {
    const double t = fabs(x / 3.75) * fabs(x / 3.75);

    if (t < 1.) {
        return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.360768e-1 + t * 0.45813e-2)))));
    } else {
        const double ax = abs(x);
        const double y  = 3.75 / ax;
        return (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
    }
}

static inline double besseli1(double x) {
    const double t = fabs(x / 3.75) * fabs(x / 3.75);
    const double ax = abs(x);

    if (t < 1.) {
        return ax * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934 + t * (0.2658733e-1 + t * (0.301532e-2 + t * 0.32411e-3))))));
    } else {
        const double y = 3.75 / ax;
        return (exp(ax) / sqrt(ax)) * (0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * (0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2))))))));
    }
}

static double besselk0(double x) {
    if (x <= 2.0) {
        const double y = x * x / 4.0;
        return (-log(x / 2.0) * besseli0(x)) + (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))));
    } else {
        const double y = 2.0 / x;
        return (exp(-x) / sqrt(x)) * (1.25331414 + y * (-0.7832358e-1 + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
    }
}

static inline double besselk1(double x) {
    if (x <= 2.0) {
        const double y = x * x / 4.0;
        return (log(x / 2.0) * besseli1(x)) + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
    } else {
        const double y = 2.0 / x;
        return (exp(-x) / sqrt(x)) * (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
}

#else

static inline double poly(const double *coef, const int n, const double x) {
    double ans = coef[n];
    for (int i=n-1;i>=0;i--) ans = ans*x+coef[i];
    return ans;
}

const double c_i0p[14] = {9.999999999999997e-1, 2.466405579426905e-1, 1.478980363444585e-2, 3.826993559940360e-4, 5.395676869878828e-6, 4.700912200921704e-8, 2.733894920915608e-10, 1.115830108455192e-12, 3.301093025084127e-15, 7.209167098020555e-18, 1.166898488777214e-20, 1.378948246502109e-23, 1.124884061857506e-26, 5.498556929587117e-30};
const double c_i0q[5]  = {4.463598170691436e-1, 1.702205745042606e-3, 2.792125684538934e-6, 2.369902034785866e-9, 8.965900179621208e-13};
const double c_i0pp[5] = {1.192273748120670e-1, 1.947452015979746e-1, 7.629241821600588e-2, 8.474903580801549e-3, 2.023821945835647e-4};
const double c_i0qq[6] = {2.962898424533095e-1, 4.866115913196384e-1, 1.938352806477617e-1, 2.261671093400046e-2, 6.450448095075585e-4, 1.529835782400450e-6};

static inline double besseli0(double x) {
    const double ax = abs(x);
    if (ax < 15.0) {
        const double y = x * x;
        return poly(c_i0p, 13, y) / poly(c_i0q, 4, 225. - y);
    } else {
        const double z = 1.0 - 15.0 / ax;
        return exp(ax) * poly(c_i0pp, 4, z) / (poly(c_i0qq, 5, z) * sqrt(ax));
    }
}

const double c_i1p[14] = {5.000000000000000e-1, 6.090824836578078e-2, 2.407288574545340e-3, 4.622311145544158e-5, 5.161743818147913e-7, 3.712362374847555e-9, 1.833983433811517e-11, 6.493125133990706e-14, 1.693074927497696e-16, 3.299609473102338e-19, 4.813071975603122e-22, 5.164275442089090e-25, 3.846870021788629e-28, 1.712948291408736e-31};
const double c_i1q[5]  = {4.665973211630446e-1, 1.677754477613006e-3, 2.583049634689725e-6, 2.045930934253556e-9, 7.166133240195285e-13};
const double c_i1pp[5] = {1.286515211317124e-1, 1.930915272916783e-1, 6.965689298161343e-2, 7.345978783504595e-3, 1.963602129240502e-4};
const double c_i1qq[6] = {3.309385098860755e-1, 4.878218424097628e-1, 1.663088501568696e-1, 1.473541892809522e-2, 1.964131438571051e-4, -1.034524660214173e-6};

static inline double besseli1(const double x) {
    const double ax = abs(x);
    if (ax < 15.0) {
        const double y = x * x;
        return x * poly(c_i1p, 13, y) / poly(c_i1q, 4, 225. - y);
    } else {
        const double z   = 1.0 - 15.0 / ax;
        double       ans = exp(ax) * poly(c_i1pp, 4, z) / (poly(c_i1qq, 5, z) * sqrt(ax));
        return x > 0.0 ? ans : -ans;
    }
}

const double c_k0pi[5] = {1.0, 2.346487949187396e-1, 1.187082088663404e-2, 2.150707366040937e-4, 1.425433617130587e-6};
const double c_k0qi[3] = {9.847324170755358e-1, 1.518396076767770e-2, 8.362215678646257e-5};
const double c_k0p[5]  = {1.159315156584126e-1, 2.770731240515333e-1, 2.066458134619875e-2, 4.574734709978264e-4, 3.454715527986737e-6};
const double c_k0q[3]  = {9.836249671709183e-1, 1.627693622304549e-2, 9.809660603621949e-5};
const double c_k0pp[8] = {1.253314137315499, 1.475731032429900e1, 6.123767403223466e1, 1.121012633939949e2, 9.285288485892228e1, 3.198289277679660e1, 3.595376024148513, 6.160228690102976e-2};
const double c_k0qq[8] = {1.0, 1.189963006673403e1, 5.027773590829784e1, 9.496513373427093e1, 8.318077493230258e1, 3.181399777449301e1, 4.443672926432041, 1.408295601966600e-1};

static inline double besselk0(const double x) {
    if (x <= 1.0) {
        const double z = x * x;
        const double term = poly(c_k0pi, 4, z) * log(x) / poly(c_k0qi, 2, 1. - z);
        return poly(c_k0p, 4, z) / poly(c_k0q, 2, 1. - z) - term;
    } else {
        const double z = 1. / x;
        return exp(-x) * poly(c_k0pp, 7, z) / (poly(c_k0qq, 7, z) * sqrt(x));
    }
}
const double c_k1pi[5] = {0.5, 5.598072040178741e-2, 1.818666382168295e-3, 2.397509908859959e-5, 1.239567816344855e-7};
const double c_k1qi[3] = {9.870202601341150e-1, 1.292092053534579e-2, 5.881933053917096e-5};
const double c_k1p[5]  = {-3.079657578292062e-1, -8.109417631822442e-2, -3.477550948593604e-3, -5.385594871975406e-5, -3.110372465429008e-7};
const double c_k1q[3]  = {9.861813171751389e-1, 1.375094061153160e-2, 6.774221332947002e-5};
const double c_k1pp[8] = {1.253314137315502, 1.457171340220454e1, 6.063161173098803e1, 1.147386690867892e2, 1.040442011439181e2, 4.356596656837691e1, 7.265230396353690, 3.144418558991021e-1};
const double c_k1qq[8] = {1.0, 1.125154514806458e1, 4.427488496597630e1, 7.616113213117645e1, 5.863377227890893e1, 1.850303673841586e1, 1.857244676566022, 2.538540887654872e-2};

static inline double besselk1(const double x) {
    if (x <= 1.0) {
        double const z    = x * x;
        double const term = poly(c_k1pi, 4, z) * log(x) / poly(c_k1qi, 2, 1. - z);
        return x * (poly(c_k1p, 4, z) / poly(c_k1q, 2, 1. - z) + term) + 1. / x;
    } else {
        double const z = 1.0 / x;
        return exp(-x) * poly(c_k1pp, 7, z) / (poly(c_k1qq, 7, z) * sqrt(x));
    }
}

#endif