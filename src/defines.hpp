#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <cmath>
#include <iostream>
#include <cassert>

#define DIM 2


#define C2PI 2.0*M_PI

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
    #define INFO3(a,b) ((void)0);
#endif

#ifdef VERBOSE
    #define INFOLOG(a)    printf(a);
#else
    #define INFOLOG(a)    ((void)0);
#endif

//-----------------------------------------------------------------------------
#ifndef NDEBUG
    #define ERROR(a)  printf(a);
#else
    #define ERROR(a) ((void)0);
#endif

#endif