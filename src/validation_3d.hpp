/**
 * @file validation.hpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-19
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#ifndef VALIDATION_3D_HPP
#define VALIDATION_3D_HPP

#include "Solver.hpp"
#include "defines.hpp"
#include "expint.hpp"

#include "hdf5_io.hpp"
#include "mpi.h"

#define MANUFACTURED_SOLUTION


struct DomainDescr {
    int          nglob[3]   = {64, 64, 64};
    int          nproc[3]   = {2, 2, 1};
    double       L[3]       = {1.0, 1.0, 1.0};
    FLUPS::BoundaryType mybc[3][2] = {{FLUPS::UNB, FLUPS::UNB}, {FLUPS::UNB, FLUPS::UNB}, {FLUPS::UNB, FLUPS::UNB}};
    // BoundaryType mybc[3][2] = {{UNB, UNB}, {UNB, UNB}, {UNB, UNB}};
};

/**
 * @name validation of the solver using a gaussian blob
 * 
 */
/**@{ */
void validation_3d(const DomainDescr myCase, const FLUPS::SolverType type, const FLUPS::GreenType typeGreen);
void validation_3d(const DomainDescr myCase, const FLUPS::SolverType type, const FLUPS::GreenType typeGreen, const int nSolve);
/**@} */


#ifdef MANUFACTURED_SOLUTION


struct manuParams {
    double       freq    = 1; //an integer or 0.5
    double       sign[2] = {0., 0.};
    double       sigma   = 0.15;
    double       center  = 0.5;
};

typedef double (*manuF)(const double, const double, const manuParams);


static inline double fOddOdd(const double x, const double L, const manuParams params) {
    return sin((c_2pi / L * params.freq) * x);
}
static inline double d2dx2_fOddOdd(const double x, const double L, const manuParams params) {
    return -(c_2pi / L * params.freq) * (c_2pi / L * params.freq) * sin((c_2pi / L * params.freq) * x);
}

static inline double fEvenEven(const double x, const double L, const manuParams params) {
    return cos((M_PI / L * params.freq) * x);
}
static inline double d2dx2_fEvenEven(const double x, const double L, const manuParams params) {
    return -(M_PI / L * params.freq) * (M_PI / L * params.freq) * cos((M_PI / L * params.freq) * x);
}

static inline double fOddEven(const double x, const double L, const manuParams params) {
    return sin((M_PI / L * (params.freq+.5)) * x);
}
static inline double d2dx2_fOddEven(const double x, const double L, const manuParams params) {
    return -(M_PI / L * (params.freq+.5)) * (M_PI / L * (params.freq+.5)) * sin((M_PI / L * (params.freq+.5)) * x);
}

static inline double fEvenOdd(const double x, const double L, const manuParams params) {
    return cos((M_PI / L * (params.freq+.5)) * x);
}
static inline double d2dx2_fEvenOdd(const double x, const double L, const manuParams params) {
    return -(M_PI / L * (params.freq+.5)) * (M_PI / L * (params.freq+.5)) * cos((M_PI / L * (params.freq+.5)) * x);
}

static inline double fUnb(const double x, const double L, const manuParams params) {
    const double x0 = (x -       params.center  * L) / params.sigma;
    const double x1 = (x +       params.center  * L) / params.sigma;
    const double x2 = (x - (2. - params.center) * L) / params.sigma;
    return                   exp(-x0*x0) + \
            params.sign[0] * exp(-x1*x1) + \
            params.sign[1] * exp(-x2*x2) ;
}
static inline double d2dx2_fUnb(const double x, const double L, const manuParams params) {
    const double x0 = (x -       params.center  * L) / params.sigma;
    const double x1 = (x +       params.center  * L) / params.sigma;
    const double x2 = (x - (2. - params.center) * L) / params.sigma;
    return                   -2. / (params.sigma * params.sigma) * exp(-x0*x0 ) * (1. - 2. * ( x0 * x0)) + \
            params.sign[0] * -2. / (params.sigma * params.sigma) * exp(-x1*x1 ) * (1. - 2. * ( x1 * x1)) + \
            params.sign[1] * -2. / (params.sigma * params.sigma) * exp(-x2*x2 ) * (1. - 2. * ( x2 * x2)) ;
}


static const inline double c_C = 10.;
static const inline double c_sigma = .5;

static inline double fUnbSpietz(const double x, const double L, const manuParams params) {
    const double x0 = (x -       params.center  * L) / c_sigma;
    const double x1 = (x +       params.center  * L) / c_sigma;
    const double x2 = (x - (2. - params.center) * L) / c_sigma;
    return    fabs(x0)>=1. ?  0.0 : exp(c_C * (1. - 1. / (1. - x0 * x0))); //+ \   
            // params.sign[0] * exp(c_C * (1 - 1 / (1 - x1 * x1))) + \
            // params.sign[1] * exp(c_C * (1 - 1 / (1 - x2 * x2)));
}
static inline double d2dx2_fUnbSpietz(const double x, const double L, const manuParams params) {
    const double x0sq = pow((x -       params.center  * L) / c_sigma,2);
    const double x1sq = pow((x +       params.center  * L) / c_sigma,2);
    const double x2sq = pow((x - (2. - params.center) * L) / c_sigma,2);
    return    fabs(x0sq)>=1. ?  0.0 :              exp(c_C * (1. - 1. / (1. - x0sq))) * (2. * c_C * (2. * (c_C - 1.) * x0sq + 3. * x0sq * x0sq - 1.)) / pow(x0sq - 1., 4) / (c_sigma*c_sigma);// +
        //    params.sign[0] * exp(c_C * (1 - 1 / (1 - x1sq))) * (2 * c_C * (2 * (c_C - 1) * x1sq + 3 * x1sq * x1sq - 1)) / pow(1 - x1sq, 4) +
        //    params.sign[1] * exp(c_C * (1 - 1 / (1 - x2sq))) * (2 * c_C * (2 * (c_C - 1) * x2sq + 3 * x2sq * x2sq - 1)) / pow(1 - x2sq, 4);
}


static inline double fCst(const double x, const double L, const manuParams params) {
    return 1.0;
}
static inline double fZero(const double x, const double L, const manuParams params) {
    return 0.0;
}

#endif


#endif