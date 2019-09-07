/**
 * @file green_functions.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-17
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "expint.hpp"
#include "green_functions_2d.hpp"

/**
 * @brief 2D Green function for 3 dir unbounded form Chatelain et al.
 * 
 * The Green function is given by
 * 
 * \f$ \frac{1}{4 \pi} log\left( r^2 \right) \f$
 * 
 * and its value in 0 is given by
 * 
 * \f$ \frac{2}{\pi^{3/2} \, r_{eq}} log\left( 1+\sqrt{2} \right) \f$
 * 
 * with \f$ r_{eq} = \sqrt{\frac{h_X h_Y}{\pi}} \f$
 * 
 * see Chatelain et al. (JCP, 2010) for more details
 * 
 * @param size 
 * @param kfact 
 * @param green 
 */
void Green_2D_3dirunbounded_0dirspectral_chat2(const int size[3],const double hfact[3],double* green){
    BEGIN_FUNC
    // assert that the green spacing is not 0.0 everywhere
    assert(hfact[0] != 0.0);
    assert(hfact[1] != 0.0);
    assert(hfact[2] == 0.0); // since we are 2d

    INFOLOG("assigning Chatelain second order to Green function\n");

    for(int i2=0; i2<size[2]; i2++){
        for(int i1=0; i1<size[1]; i1++){
            for(int i0=0; i0<size[0]; i0++){
                const double x = i0*hfact[0];
                const double y = i1*hfact[1];
                const double r2 = x*x+y*y;

                const size_t id = i0 + i1*size[0] + i2*size[0]*size[1];

                green[id] = c_1o4pi * log(r2);
            }
        }
    }
    // reset the value in 0.0
    const double req = sqrt(c_1opi*hfact[0]*hfact[1]);
    green[0] = - 2.0*log(1.0+M_SQRT2)/(sqrt(M_PI*M_PI*M_PI) * req) ;

    // printf("G(0,0) = %f\n",green[0]);
    // printf("G(1,0) = %f\n",green[1]);
}

/**
 * @brief 2D Green function for 3 dir unbounded regularized 2nd order, Hejlesen et al.
 * 
 The Green function is given by
 * 
 * \f$ \frac{1}{4 \pi} \left[ log\left( r^2 \right) + Ei(\frac{r^2}{2 \epsilon^2}) \right] \f$
 * 
 * and its value in 0 is given by
 * 
 * \f$ - \frac{1}{2 \pi} \left[ \frac{\gamma}{2} - log\left( \sqrt{2} \epsilon \right) \right] \f$
 * 
 * with \f$ \epsilon = \alpha h \f$
 * 
 * see Hejlesen et al. (JCP, 2013) for more details
 * 
 * @param size the size of the Green's function
 * @param hfact the mesh grid spacing
 * @param green the Green memory
 * @param alpha regularization parameter
 */
void Green_2D_3dirunbounded_0dirspectral_hej2(const int size[3],const double hfact[3],double* green,const double alpha){
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
                const double r2 = x*x+y*y;

                const size_t id = i0 + size[0]*(i1 + i2*size[1]);

                green[id] = c_1o4pi * (log(r2) + expint_ei(r2*oo2eps2));
            }
        }
    }
    // reset the value in 0.0
    green[0] = - c_1o2pi*(GAMMA*0.5 - log(M_SQRT2*eps));
}



/**
 * @brief 2D Green function for 3 dir unbounded regularized 4th order, Hejlesen et al.
 * 
 The Green function is given by
 * 
 * \f$ \frac{1}{4 \pi} \left[ log\left( r^2 \right) + Ei(\frac{r^2}{2 \epsilon^2}) - exp(\frac{r^2}{2 \epsilon^2}) \right] \f$
 * 
 * and its value in 0 is given by
 * 
 * \f$ - \frac{1}{2 \pi} \left[ \frac{\gamma}{2} - log\left( \sqrt{2} \epsilon \right) + \frac{1}{2} \right] \f$
 * 
 * with \f$ \epsilon = \alpha h \f$
 * 
 * see Hejlesen et al. (JCP, 2013) for more details
 * 
 * @param size the size of the Green's function
 * @param hfact the mesh grid spacing
 * @param green the Green memory
 * @param alpha regularization parameter
 */
void Green_2D_3dirunbounded_0dirspectral_hej4(const int size[3],const double hfact[3],double* green,const double alpha){
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
                const double r2 = x*x+y*y;

                const size_t id = i0 + size[0]*(i1 + i2*size[1]);

                green[id] = c_1o4pi * (log(r2) + expint_ei(r2*oo2eps2) - exp(-r2*oo2eps2));
            }
        }
    }
    // reset the value in 0.0
    green[0] = - c_1o2pi*(GAMMA*0.5 - log(M_SQRT2*eps) + c_1o2);
}
