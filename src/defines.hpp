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

#include <cmath>
#include <iostream>
#include <cassert>

#include <execinfo.h>
#include "fftw3.h"

#define DIM 3

#define FFTW_FLAG FFTW_MEASURE

//-----------------------------------------------------------------------------
// need to add __LINE__ ; __FILE__ ; __func__
#ifdef VERBOSE
    #define INFO(a)         printf(a);
    #define INFO2(a,b)      printf(a,b);
    #define INFO3(a,b,c)    printf(a,b,c);
#else
    #define INFO(a)    ((void)0);
    #define INFO2(a,b) ((void)0);
    #define INFO3(a,b,c) ((void)0);
#endif

#ifdef VERBOSE
    #define INFOLOG(a)        printf(a);
    #define INFOLOG2(a,b)     printf(a,b);
    #define INFOLOG3(a,b,c)   printf(a,b,c);
    #define INFOLOG4(a,b,c,d) printf(a,b,c,d);

    #define BEGIN_FUNC ((void)0); 
    // #define BEGIN_FUNC {INFOLOG4(">>entering %s from %s at line %d\n",__func__,__FILE__,__LINE__);}
#else
    #define INFOLOG(a)    ((void)0);
    #define INFOLOG2(a,b)    ((void)0);
    #define INFOLOG3(a,b,c)    ((void)0);
    #define INFOLOG4(a,b,c,d)    ((void)0);
    #define BEGIN_FUNC  ((void)0);
#endif

//-----------------------------------------------------------------------------
#ifndef NDEBUG
    // #define UP_ERROR(a) {\
    //                     printf("ERROR - %s\n",a);\
    //                     void* callstack[128];\
    //                     int i, frames = backtrace(callstack, 128);\
    //                     char** strs = backtrace_symbols(callstack, frames);\
    //                     for (i = 0; i < frames; ++i) printf("   %s\n", strs[i]);\
    //                     free(strs);\
    //                 }

    #define UP_ERROR(a) { char msg_error[1024]; sprintf(msg_error,"ERROR - %s - in %s from %s at line %d\n",a,__func__,__FILE__,__LINE__); printf(msg_error);}

    #define UP_CHECK0(a,b)        if(!(a)){ UP_ERROR(b);}
    #define UP_CHECK1(a,b,c)      if(!(a)){ char msg_chk[1024]; sprintf(msg_chk,b,c); UP_ERROR(msg_chk);}
    #define UP_CHECK2(a,b,c,d)    if(!(a)){ char msg_chk[1024]; sprintf(msg_chk,b,c,d); UP_ERROR(msg_chk);}
    #define UP_CHECK3(a,b,c,d,e)  if(!(a)){ char msg_chk[1024]; sprintf(msg_chk,b,c,d,e); UP_ERROR(msg_chk);}
#else
    #define UP_ERROR(a) ((void)0);

    #define UP_CHECK0(a,b)        ((void)0);
    #define UP_CHECK1(a,b,c)      ((void)0);
    #define UP_CHECK2(a,b,c,d)    ((void)0);
    #define UP_CHECK3(a,b,c,d,e)  ((void)0);
#endif

#define GAMMA 0.5772156649015328606


#define UP_ALIGNMENT 16
#define UP_ISALIGNED(a) {((uintptr_t)(const void*) a)%UP_ALIGNMENT == 0 ;}

typedef double* __restrict __attribute__((aligned (UP_ALIGNMENT))) opt_double_ptr;
typedef fftw_complex* __restrict __attribute__((aligned (UP_ALIGNMENT))) opt_complex_ptr;

static const double c_1opi  = 1.0/(1.0*M_PI);
static const double c_1o2pi = 1.0/(2.0*M_PI);
static const double c_1o4pi = 1.0/(4.0*M_PI);
static const double c_1o2   = 1.0/2.0;
static const double c_1o4   = 1.0/4.0;

static const double c_1osqrt2 = 1.0/M_SQRT2;

static const double c_2pi = 2.0*M_PI;


#endif