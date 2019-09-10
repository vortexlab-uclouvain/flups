/**
 * @file defines.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <cassert>
#include <cmath>
#include <iostream>

#include <execinfo.h>
#include "fftw3.h"
#include "mpi.h"

#define GREEN_DIM 3

// #define DUMP_H5
#undef DUMP_H5

#define FFTW_FLAG FFTW_MEASURE

//=============================================================================
// WARNINGS
//=============================================================================
static inline void FLUPS_WARNING(char* a) {
    char tmp[512];
    sprintf(tmp, a);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};
template<typename T1>
static inline void FLUPS_WARNING(char* a, T1 b) {
    char tmp[512];
    sprintf(tmp, a, b);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2>
static inline void FLUPS_WARNING(char* a, T1 b, T2 c) {
    char tmp[512];
    sprintf(tmp, a, b, c);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_WARNING(char* a, T1 b, T2 c, T3 d) {
    char tmp[512];
    sprintf(tmp, a, b, c, d);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_WARNING(char* a, T1 b, T2 c, T3 d, T4 e) {
    char tmp[512];
    sprintf(tmp, a, b, c, d, e);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_WARNING(char* a, T1 b, T2 c, T3 d, T4 e, T5 f) {
    char tmp[512];
    sprintf(tmp, a, b, c, d, e, f);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_WARNING(char* a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g) {
    char tmp[512];
    sprintf(tmp, a, b, c, d, e, f, g);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - from %s\n", tmp, __func__);
    printf(msg_error);
    fflush(stdout);
};

//=============================================================================
// LOG / FLUPS_INFO
//=============================================================================
#ifdef VERBOSE
static inline void FLUPS_INFO_DISP(char* a) {
    char msg[1024];
    sprintf(msg, "[FLUPS] %s\n", a);
    printf(msg);
    fflush(stdout);
}
static inline void FLUPS_INFO(char* a) {
    char tmp[512];
    sprintf(tmp, a);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1>
static inline void FLUPS_INFO(char* a, T1 b) {
    char tmp[512];
    sprintf(tmp, a, b);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2>
static inline void FLUPS_INFO(char* a, T1 b, T2 c) {
    char tmp[512];
    sprintf(tmp, a, b, c);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d) {
    char tmp[512];
    sprintf(tmp, a, b, c, d);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d, T4 e) {
    char tmp[512];
    sprintf(tmp, a, b, c, d, e);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d, T4 e, T5 f) {
    char tmp[512];
    sprintf(tmp, a, b, c, d, e, f);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g) {
    char tmp[512];
    sprintf(tmp, a, b, c, d, e, f, g);
    FLUPS_INFO_DISP(tmp);
};
#define BEGIN_FUNC {FLUPS_INFO(">>> entering %s from %s at line %d", __func__, __FILE__, __LINE__);};
#else
static inline void FLUPS_INFO_DISP(char* a) {
    (void(0));
}
static inline void FLUPS_INFO(char* a) {
    (void(0));
};
template <typename T1>
static inline void FLUPS_INFO(char* a, T1 b) {
    ((void)0);
};
template <typename T1,typename T2>
static inline void FLUPS_INFO(char* a, T1 b, T2 c) {
    ((void)0);
};
template <typename T1,typename T2,typename T3>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d, T4 e) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d, T4 e, T5 f) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_INFO(char* a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g) {
    ((void)0);
};
// static inline void BEGIN_FUNC() {
//     ((void)0);
// };
#define BEGIN_FUNC { ((void)0);};
#endif

//=============================================================================
// ERRORS AND ASSERTS
//=============================================================================
static inline void FLUPS_ERROR(char* a) {
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - in %s from %s at line %d\n", a, __func__, __FILE__, __LINE__);
    printf(msg_error);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,0);
};

#ifndef NDEBUG
static inline void FLUPS_CHECK(bool a, char* b) {
    if (!(a)) {
        FLUPS_ERROR(b);
    }
};
template<typename T1>
static inline void FLUPS_CHECK(bool a, char* b, T1 c) {
    if (!(a)) {
        char msg_chk[1024];
        sprintf(msg_chk, b, c);
        FLUPS_ERROR(msg_chk);
    }
};
template<typename T1,typename T2>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d) {
    if (!(a)) {
        char msg_chk[1024];
        sprintf(msg_chk, b, c, d);
        FLUPS_ERROR(msg_chk);
    }
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e) {
    if (!(a)) {
        char msg_chk[1024];
        sprintf(msg_chk, b, c, d, e);
        FLUPS_ERROR(msg_chk);
    }
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e, T4 f) {
    if (!(a)) {
        char msg_chk[1024];
        sprintf(msg_chk, b, c, d, e, f);
        FLUPS_ERROR(msg_chk);
    }
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e, T4 f, T5 g) {
    if (!(a)) {
        char msg_chk[1024];
        sprintf(msg_chk, b, c, d, e, f, g);
        FLUPS_ERROR(msg_chk);
    }
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e, T4 f, T5 g, T6 h) {
    if (!(a)) {
        char msg_chk[1024];
        sprintf(msg_chk, b, c, d, e, f, g, h);
        FLUPS_ERROR(msg_chk);
    }
};
#else
static inline void FLUPS_CHECK(bool a, char* b) {
    ((void)0);
};
template<typename T1>
static inline void FLUPS_CHECK(bool a, char* b, T1 c) {
    (void(0));
};
template<typename T1,typename T2>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d) {
    (void(0));
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e) {
    (void(0));
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e, T4 f) {
    (void(0));
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e, T4 f, T5 g) {
    (void(0));
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_CHECK(bool a, char* b, T1 c, T2 d, T3 e, T4 f, T5 g, T6 h) {
    (void(0));
};
#endif

    //=============================================================================
    // CONSTANTS AND OTHERS
    //=============================================================================

#define GAMMA 0.5772156649015328606

#define FLUPS_FORWARD -1  // = FFTW_FORWARD
#define FLUPS_BACKWARD 1  // = FFTW_BACKWARD

#define FLUPS_ALIGNMENT 32

template <typename T>
static inline bool FLUPS_ISALIGNED(T a) {
    return ((uintptr_t)(const void*)a) % FLUPS_ALIGNMENT == 0;
}

typedef int* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_int_ptr;
typedef double* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_double_ptr;
typedef fftw_complex* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_complex_ptr;

namespace FLUPS {
/**
     * @brief the boundary condition can be EVEN, ODD, PERiodic or UNBounded
     * 
     */
enum BoundaryType {
    EVEN = 0, /**< EVEN boundary condition = zero flux  */
    ODD  = 1, /**< ODD boundary condition = zero value */
    PER  = 3, /**< PERiodic boundary conditions */
    UNB  = 4  /**< UNBounded boundary condition */
};

/**
     * @brief The type of Green's function used for the Poisson solver
     * 
     */
enum GreenType {
    CHAT_2, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2,  /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2,  /**< @brief regularized in zero, order 2, Hejlesen et al. (2015)*/
    HEJ_4,  /**< @brief regularized in zero, order 4, Hejlesen et al. (2015)*/
    HEJ_6,  /**< @brief regularized in zero, order 6, Hejlesen et al. (2015)*/
};

/**
     * @brief Type of Poisson equation solved
     * 
     */
enum SolverType {
    SRHS, /**<@brief scalar \f$ \nabla^2 f = rhs \f$ */
    VRHS, /**<@brief vectorial \f$ \nabla^2 f = rhs \f$ */
    ROT,  /**<@brief vectorial \f$ \nabla^2 f = \nabla \times rhs \f$ */
    DIV   /**<@brief scalar \f$ \nabla^2 f = \nabla \cdot rhs \f$ */
};

class Solver;
class FFTW_plan_dim;
class Topology;
class SwitchTopo;
}  // namespace FLUPS

static const double c_1opi     = 1.0 / (1.0 * M_PI);
static const double c_1o2pi    = 1.0 / (2.0 * M_PI);
static const double c_1o4pi    = 1.0 / (4.0 * M_PI);
static const double c_1osqrtpi = 1.0 / sqrt(M_PI);
static const double c_1o2      = 1.0 / 2.0;
static const double c_1o4      = 1.0 / 4.0;
static const double c_7o4      = 7. / 4.;
static const double c_19o8     = 19. / 8;
static const double c_2o3      = 2. / 3;
static const double c_1o24     = 1. / 24;

static const double c_1osqrt2 = 1.0 / M_SQRT2;

static const double c_2pi = 2.0 * M_PI;

#endif