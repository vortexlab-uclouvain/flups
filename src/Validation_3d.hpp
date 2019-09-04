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

#include "FFTW_Solver.hpp"
#include "defines.hpp"
#include "expint.hpp"

#include "hdf5_io.hpp"
#include "mpi.h"

#define MANUFACTURED_SOLUTION

/**
 * @brief everything required to characterize the computational domain and initial condition 
 * 
 */
struct DomainDescr {
    int          nglob[3]   = {64, 64, 64};
    int          nproc[3]   = {2, 2, 1};
    double       L[3]       = {1.0, 1.0, 1.0};
    double       sigma      = 0.1;              // smoothing length scale of the blob
    double       center[3]  = {0.5, 0.5, 0.5};  //center of the blob (in % of L)
    BoundaryType mybc[3][2] = {{UNB, UNB}, {UNB, UNB}, {UNB, UNB}};
};

/**
 * @name validation of the solver using a gaussian blob
 * 
 */
/**@{ */
void validation_3d(const DomainDescr myCase, const SolverType type, const GreenType typeGreen);
void validation_3d(const DomainDescr myCase, const SolverType type, const GreenType typeGreen, const int nSolve);
/**@} */


#ifdef MANUFACTURED_SOLUTION


struct manuParams {
    double       freq    = 1; //an integer or 0.5
    double       sign[2] = {0, 0};
    double       sigma   = 0.15;
};

typedef double (*manuF)(const double, const double, const manuParams);


static inline double fOddOdd(const double x, const double L, const manuParams params) {
    return sin((c_2pi / L * params.freq) * x);
}
static inline double d2dx2_fOddOdd(const double x, const double L, const manuParams params) {
    return -(c_2pi / L * params.freq) * (c_2pi / L * params.freq) * sin((c_2pi / L * params.freq) * x);
}

static inline double fEvenEven(const double x, const double L, const manuParams params) {
    return cos((c_2pi / L * params.freq) * x);
}
static inline double d2dx2_fEvenEven(const double x, const double L, const manuParams params) {
    return -(c_2pi / L * params.freq) * (c_2pi / L * params.freq) * cos((c_2pi / L * params.freq) * x);
}

static inline double fUnb(const double x, const double L, const manuParams params) {
    return                   exp(-(x - .5 * L) * (x - .5 * L) / (params.sigma * params.sigma)) + \
            params.sign[0] * exp(-(x + .5 * L) * (x + .5 * L) / (params.sigma * params.sigma)) + \
            params.sign[1] * exp(-(x -1.5 * L) * (x -1.5 * L) / (params.sigma * params.sigma)) ;
}
static inline double d2dx2_fUnb(const double x, const double L, const manuParams params) {
    return                   -2. / (params.sigma * params.sigma) * exp(-(x - .5*L) * (x - .5*L) / (params.sigma * params.sigma)) * (1. - 2. * ((x - .5*L) * (x - .5*L) / (params.sigma * params.sigma))) + \
            params.sign[0] * -2. / (params.sigma * params.sigma) * exp(-(x + .5*L) * (x + .5*L) / (params.sigma * params.sigma)) * (1. - 2. * ((x + .5*L) * (x + .5*L) / (params.sigma * params.sigma))) + \
            params.sign[1] * -2. / (params.sigma * params.sigma) * exp(-(x -1.5*L) * (x -1.5*L) / (params.sigma * params.sigma)) * (1. - 2. * ((x -1.5*L) * (x -1.5*L) / (params.sigma * params.sigma))) ;
}

#endif


#endif