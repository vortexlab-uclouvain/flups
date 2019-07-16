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
    #define ERROR(a)  printf(a);
#else
    #define ERROR(a) ((void)0);
#endif




#endif