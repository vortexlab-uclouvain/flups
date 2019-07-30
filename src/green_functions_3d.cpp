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

void Green_3D_3dirunbounded_0dirspectral_chat2(const Topology* topo,const double hfact[3],double* green){
    /**
     * @todo get the 3dir unbounded green function
     */

    UP_CHECK0(false,"non implemented function");
}


void Green_3D_3dirunbounded_0dirspectral_hej2 (const Topology* topo,const double hfact[3],double* green,const double alpha){
    BEGIN_FUNC

    // assert that the green spacing is not 0.0 everywhere
    UP_CHECK0(hfact[0] != 0.0,"grid spacing cannot be 0");
    UP_CHECK0(hfact[1] != 0.0,"grid spacing cannot be 0");
    UP_CHECK0(hfact[2] != 0.0,"grid spacing cannot be 0");

    int ax0 = topo->axis();
    int ax1 = (ax0+1)%3;
    int ax2 = (ax0+2)%3;

    const double eps     = alpha*hfact[0];
    const double oo2eps2 = 1.0/(2.0*eps*eps);

    for(int i2=0; i2<topo->nloc(ax2); i2++){
        for(int i1=0; i1<topo->nloc(ax1); i1++){
            for(int i0=0; i0<topo->nloc(ax0); i0++){

                const double x0 = i0*hfact[ax0];
                const double x1 = i1*hfact[ax1];
                const double x2 = i2*hfact[ax2];

                const double r = sqrt(x0*x0+x1*x1+x2*x2);

                const size_t id = i0 + topo->nloc(ax0)*(i1 + i2*topo->nloc(ax1));

                green[id] = - c_1o4pi * 1.0/r * (erf(r*c_1osqrt2/eps));
            }
        }
    }
    // reset the value in 0.0
    green[0] = - M_SQRT2/(4.0* eps *sqrt(M_PI*M_PI*M_PI));
}

void Green_3D_3dirunbounded_0dirspectral_hej4 (const Topology* topo,const double hfact[3],double* green,const double alpha){
    BEGIN_FUNC

    /**
     * @todo get the 3dir unbounded green function
     */

    UP_CHECK0(false,"non implemented function");
    
}