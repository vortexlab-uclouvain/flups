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

#define DIM 2

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

    #define BEGIN_FUNC ((void)0); //INFOLOG4(">>entering %s from %s at line %d\n",__func__,__FILE__,__LINE__);
#else
    #define INFOLOG(a)    ((void)0);
    #define INFOLOG2(a,b)    ((void)0);
    #define INFOLOG3(a,b,c)    ((void)0);
    #define INFOLOG4(a,b,c,d)    ((void)0);
    #define BEGIN_FUNC  ((void)0);
#endif

//-----------------------------------------------------------------------------
#ifndef NDEBUG
    #define UP_ERROR(a) {\
                        printf("ERROR - %s\n",a);\
                        void* callstack[128];\
                        int i, frames = backtrace(callstack, 128);\
                        char** strs = backtrace_symbols(callstack, frames);\
                        for (i = 0; i < frames; ++i) printf("   %s\n", strs[i]);\
                        free(strs);\
                    }

    #define UP_CHECK0(a,b)        if(!(a)){ UP_ERROR(b);}
    #define UP_CHECK1(a,b,c)      if(!(a)){ char msg[256]; sprintf(msg,b,c); UP_ERROR(msg);}
    #define UP_CHECK2(a,b,c,d)    if(!(a)){ char msg[256]; sprintf(msg,b,c,d); UP_ERROR(msg);}
    #define UP_CHECK3(a,b,c,d,e)  if(!(a)){ char msg[256]; sprintf(msg,b,c,d,e); UP_ERROR(msg);}
#else
    #define ERROR(a) ((void)0);
#endif

#define GAMMA 0.5772156649015328606



static const double c_1opi  = 1.0/(1.0*M_PI);
static const double c_1o2pi = 1.0/(2.0*M_PI);
static const double c_1o4pi = 1.0/(4.0*M_PI);
static const double c_1o2   = 1.0/2.0;
static const double c_1o4   = 1.0/4.0;


static const double c_2pi = 2.0*M_PI;


#endif