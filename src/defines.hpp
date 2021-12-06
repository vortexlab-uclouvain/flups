/**
 * @file defines.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
#include "flups.h"


#ifdef NODE_CENTERED
#define FLUPS_CELL_CENTERED 0
#else 
#define FLUPS_CELL_CENTERED 1
#endif

//=============================================================================
// LOCATORS
//=============================================================================
#define LOCATION ("in " + std::string(__func__) + " from  " + std::string(__FILE__) + " at line " + std::to_string(__LINE__) )
#define LOC      ("in " + std::string(__func__) )

//=============================================================================
// WARNINGS
//=============================================================================
static inline void FLUPS_WARNING(std::string a, std::string loc) {
    char tmp[512];
    sprintf(tmp, "%s", a.c_str());
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};
template<typename T1>
static inline void FLUPS_WARNING(std::string a, T1 b, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};
template<typename T1,typename T2>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, T4 e, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fflush(stdout);
};

//=============================================================================
// LOG / FLUPS_INFO
//=============================================================================
#ifdef VERBOSE
static inline void FLUPS_INFO_DISP(std::string a) {
    char msg[1024];
    sprintf(msg, "[FLUPS] %s\n", a.c_str());
    printf("%s", msg);
    fflush(stdout);
}
static inline void FLUPS_INFO(std::string a) {
    char tmp[512];
    sprintf(tmp, "%s", a.c_str());
    FLUPS_INFO_DISP(tmp);
};
template <typename T1>
static inline void FLUPS_INFO(std::string a, T1 b) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, T7 h) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g, h);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, T7 h, T8 i) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g, h, i);
    FLUPS_INFO_DISP(tmp);
};
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, T7 h, T8 i, T9 j) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g, h, i, j);
    FLUPS_INFO_DISP(tmp);
};
#else
static inline void FLUPS_INFO_DISP(std::string a) {
    (void(0));
}
static inline void FLUPS_INFO(std::string a) {
    (void(0));
};
template <typename T1>
static inline void FLUPS_INFO(std::string a, T1 b) {
    ((void)0);
};
template <typename T1,typename T2>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c) {
    ((void)0);
};
template <typename T1,typename T2,typename T3>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6,typename T7>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, T7 h) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6,typename T7,typename T8>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, T7 h, T8 i) {
    ((void)0);
};
template <typename T1,typename T2,typename T3,typename T4,typename T5,typename T6,typename T7,typename T8, typename T9>
static inline void FLUPS_INFO(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, T7 h, T8 i, T9 j) {
    ((void)0);
};
#endif

//=============================================================================
// VERBOSITY LEVELS
//=============================================================================
#if VERBOSE>=2
    #ifdef PROF
        #define BEGIN_FUNC double T0 = MPI_Wtime();\
                            FLUPS_INFO(">>> entering %s from %s at line %d", __func__, __FILE__, __LINE__);
        #define END_FUNC double T1 = MPI_Wtime();\
                            FLUPS_INFO(">>> leaving %s from %s at line %d after %lf [s]", __func__, __FILE__, __LINE__,(T1)-(T0));
    #else
        #define BEGIN_FUNC {FLUPS_INFO(">>> entering %s from %s at line %d", __func__, __FILE__, __LINE__);};
        #define END_FUNC {FLUPS_INFO(">>> leaving %s from %s at line %d", __func__, __FILE__, __LINE__);};
    #endif
#else
    #define BEGIN_FUNC { ((void)0);};
    #define END_FUNC { ((void)0);};
#endif


#if VERBOSE>=1
    #define FLUPS_INFO_1(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_1(...) ((void)0);
#endif
#if VERBOSE>=2
    #define FLUPS_INFO_2(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_2(...) ((void)0);
#endif
#if VERBOSE>=3
    #define FLUPS_INFO_3(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_3(...) ((void)0);
#endif
#if VERBOSE>=4
    #define FLUPS_INFO_4(...) FLUPS_INFO(__VA_ARGS__)
#else
    #define FLUPS_INFO_4(...) ((void)0);
#endif

//=============================================================================
// ERRORS AND ASSERTS
//=============================================================================
//=============================================================================
// WARNINGS
//=============================================================================
static inline void FLUPS_ERROR(std::string a, std::string loc) {
    char tmp[512];
    sprintf(tmp, "%s", a.c_str());
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};
template<typename T1>
static inline void FLUPS_ERROR(std::string a, T1 b, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};
template<typename T1,typename T2>
static inline void FLUPS_ERROR(std::string a, T1 b, T2 c, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_ERROR(std::string a, T1 b, T2 c, T3 d, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_ERROR(std::string a, T1 b, T2 c, T3 d, T4 e, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_ERROR(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_ERROR(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf("%s", msg_error);
    fprintf(stderr, "%s", msg_error);
    fflush(stderr);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD,1);
};

#ifndef NDEBUG
static inline void FLUPS_CHECK(bool a, std::string b, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,loc);
    }
};
template<typename T1>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,c,loc);
    }
};
template<typename T1,typename T2>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,c,d,loc);
    }
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,c,d,e,loc);
    }
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, T4 f, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,c,d,e,f,loc);
    }
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, T4 f, T5 g, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,c,d,e,f,g,loc);
    }
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, T4 f, T5 g, T6 h, std::string loc) {
    if (!(a)) {
        FLUPS_ERROR(b,c,d,e,f,g,h,loc);
    }
};
#else
static inline void FLUPS_CHECK(bool a, std::string b, std::string loc) {
    ((void)0);
};
template<typename T1>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, std::string loc) {
    (void(0));
};
template<typename T1,typename T2>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, std::string loc) {
    (void(0));
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, std::string loc) {
    (void(0));
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, T4 f, std::string loc) {
    (void(0));
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, T4 f, T5 g, std::string loc) {
    (void(0));
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_CHECK(bool a, std::string b, T1 c, T2 d, T3 e, T4 f, T5 g, T6 h, std::string loc) {
    (void(0));
};
#endif

    //=============================================================================
    // CONSTANTS AND OTHERS
    //=============================================================================

template <typename T>
static inline bool FLUPS_ISALIGNED(T a) {
    return ((uintptr_t)(const void*)a) % FLUPS_ALIGNMENT == 0;
}

template <typename T>
static inline int FLUPS_CMPT_ALIGNMENT(T a) {
    return ((uintptr_t)(const void*)a) % FLUPS_ALIGNMENT;
}

typedef int* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_int_ptr;
typedef double* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_double_ptr;
typedef fftw_complex* __restrict __attribute__((aligned(FLUPS_ALIGNMENT))) opt_complex_ptr;


//=============================================================================
// MEMORY ALLOCATION AND FREE
//=============================================================================
static inline void* flups_mem_malloc(size_t size) {
#if defined(__INTEL_COMPILER)
    return _mm_malloc(size, FLUPS_ALIGNMENT);
#elif defined(__GNUC__)
    void* data;
    posix_memalign(&data, FLUPS_ALIGNMENT, size);
    return data;
#endif
}

static inline void flups_mem_free(void* data) {
#if defined(__INTEL_COMPILER)
    _mm_free(data);
#elif defined(__GNUC__)
    free(data);    
#endif
}

#if defined(__INTEL_COMPILER)
    #define FLUPS_ASSUME_ALIGNED(a,b) __assume_aligned(a,b)
#elif defined(__GNUC__)
    #define FLUPS_ASSUME_ALIGNED(a,b) __builtin_assume_aligned(a,b)
#endif

typedef enum FLUPS_BoundaryType BoundaryType;
typedef enum FLUPS_GreenType    GreenType;

static const double c_1opi     = 1.0 / (1.0 * M_PI);
static const double c_1o2pi    = 1.0 / (2.0 * M_PI);
static const double c_1o2pi2   = 1.0 / (2.0 * M_PI*M_PI);
static const double c_1o4pi    = 1.0 / (4.0 * M_PI);
static const double c_1osqrtpi = 1.0 / sqrt(M_PI);
static const double c_sqrtpi   = sqrt(M_PI);
static const double c_1o2      = 1.0 / 2.0;
static const double c_1o4      = 1.0 / 4.0;
static const double c_7o4      = 7. / 4.;
static const double c_19o8     = 19. / 8;
static const double c_2o3      = 2. / 3;
static const double c_1o12     = 1. / 12.;
static const double c_1o24     = 1. / 24;
static const double c_5o16     = 5. / 16.;
static const double c_1o16     = 1. / 16.;
static const double c_187o64   = 187. / 64.;
static const double c_233o192  = 233. / 192.;
static const double c_29o192   = 29. / 192.;
static const double c_13o192   = 13. / 192.;
static const double c_1o192    = 1. / 192.;
static const double c_11o12    = 11. / 12.;
static const double c_7o24     = 7. / 24.;
static const double c_25o24    = 25. / 24.;
static const double c_1o48     = 1. / 48.;
static const double c_23o48    = 23. / 48.;
static const double c_1o384    = 1. / 384.;
static const double c_17o384   = 17. / 384.;
static const double c_23o384   = 23. / 384.;
static const double c_11o32    = 11. / 32.;
static const double c_1o96     = 1. / 96.;
static const double c_35o128   = 35. / 128.;
static const double c_47o128   = 47. / 128.;
static const double c_93o256   = 93. / 256.;
static const double c_47o256   = 47. / 256.;
static const double c_1o256    = 1. / 256.;
static const double c_73o768   = 73. / 768.;
static const double c_11o768   = 11. / 768.;
static const double c_23o768   = 23. / 768.;
static const double c_1o768    = 1. / 768.;


static const double c_1osqrt2 = 1.0 / M_SQRT2;

static const double c_2pi = 2.0 * M_PI;

#endif