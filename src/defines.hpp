/**
 * @file defines.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
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

#define GREEN_DIM 3

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
    sprintf(tmp, a.c_str());
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};
template<typename T1>
static inline void FLUPS_WARNING(std::string a, T1 b, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, T4 e, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};
template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
static inline void FLUPS_WARNING(std::string a, T1 b, T2 c, T3 d, T4 e, T5 f, T6 g, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str(), b, c, d, e, f, g);
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - WARNING] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fflush(stdout);
};

//=============================================================================
// LOG / FLUPS_INFO
//=============================================================================
#ifdef VERBOSE
static inline void FLUPS_INFO_DISP(std::string a) {
    char msg[1024];
    sprintf(msg, "[FLUPS] %s\n", a.c_str());
    printf(msg);
    fflush(stdout);
}
static inline void FLUPS_INFO(std::string a) {
    char tmp[512];
    sprintf(tmp, a.c_str());
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
#define BEGIN_FUNC {FLUPS_INFO(">>> entering %s from %s at line %d", __func__, __FILE__, __LINE__);};
#define END_FUNC {FLUPS_INFO(">>> leaving %s from %s at line %d", __func__, __FILE__, __LINE__);};
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
#define BEGIN_FUNC { ((void)0);};
#define END_FUNC { ((void)0);};
#endif

//=============================================================================
// ERRORS AND ASSERTS
//=============================================================================
//=============================================================================
// WARNINGS
//=============================================================================
static inline void FLUPS_ERROR(std::string a, std::string loc) {
    char tmp[512];
    sprintf(tmp, a.c_str());
    char msg_error[1024];
    sprintf(msg_error, "[FLUPS - ERROR] %s - %s\n", tmp, loc.c_str());
    printf(msg_error);
    fprintf(stderr,msg_error);
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
    printf(msg_error);
    fprintf(stderr,msg_error);
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
    printf(msg_error);
    fprintf(stderr,msg_error);
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
    printf(msg_error);
    fprintf(stderr,msg_error);
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
    printf(msg_error);
    fprintf(stderr,msg_error);
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
    printf(msg_error);
    fprintf(stderr,msg_error);
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
    printf(msg_error);
    fprintf(stderr,msg_error);
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

#define GAMMA 0.5772156649015328606

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

typedef enum FLUPS_BoundaryType BoundaryType;
typedef enum FLUPS_GreenType    GreenType;
typedef enum FLUPS_SolverType   SolverType;

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