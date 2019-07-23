/**
 * @file green_functions_3d.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "green_functions_3d.hpp"

void Green_3D_3dirunbounded_0dirspectral_chat2(const int size[3],const double hfact[3],double* green){
    UP_CHECK0(false,"non implemented function");
}


void Green_3D_3dirunbounded_0dirspectral_hej2 (const int size[3],const double hfact[3],double* green,const double alpha){
    BEGIN_FUNC

    // assert that the green spacing is not 0.0 everywhere
    assert(hfact[0] != 0.0);
    assert(hfact[1] != 0.0);
    assert(hfact[0] == hfact[1]);
    assert(hfact[2] == 0.0); // since we are 2d


    const double eps     = alpha*hfact[0];
    const double oo2eps2 = 1.0/(2.0*eps*eps);

    for(int i2=0; i2<size[2]; i2++){
        for(int i1=0; i1<size[1]; i1++){
            for(int i0=0; i0<size[0]; i0++){
                const double x = i0*hfact[0];
                const double y = i1*hfact[1];
                const double z = i2*hfact[2];

                const double r = sqrt(x*x+y*y+z*z);

                const size_t id = i0 + size[0]*(i1 + i2*size[1]);

                green[id] = - c_1o4pi * 1.0/r * (erf(r*c_1osqrt2/eps));
            }
        }
    }
    // reset the value in 0.0
    green[0] = - M_SQRT2/(4.0* eps *sqrt(M_PI*M_PI*M_PI));

}

void Green_3D_3dirunbounded_0dirspectral_hej4 (const int size[3],const double hfact[3],double* green,const double alpha){
    BEGIN_FUNC

    UP_CHECK0(false,"non implemented function");
    
}