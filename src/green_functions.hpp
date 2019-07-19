/**
 * @file green_functions.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-17
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#ifndef GREEN_FUNCTIONS_HPP
#define GREEN_FUNCTIONS_HPP

#include "defines.hpp"

/**
 * @brief The type of Green's function used for the Poisson solver
 * 
 * For the moment, only the Green's functions from Chatelain et al. (CHAT2) support #_dospectral = true
 */
enum OrderDiff{
    CHAT_2, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2,  /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2,  /**< @brief regularized in zero, order 2, Hejlesen et al. (2015)*/
    HEJ_4,  /**< @brief regularized in zero, order 4, Hejlesen et al. (2015)*/
    DIF_2,  /**< @brief finite difference order 2 for the spectral diff */
    DIF_4,  /**< @brief finite difference order 4 for the spectral diff */
    DIF_INF /**< @brief finite difference order infinite (=spectral) for the spectral diff */
};


/**
 * @name 3 dir unbounded - 0 dir spectral 
 * 
 */
/**@{ */
void Green_2D_3dirunbounded_0dirspectral_chat2(const int size[3],const double hfact[3],double* green);
// void Green_2D_3dirunbounded_0dirspectral_lgf2 (const size_t size[3],const double hfact[3],double* green);
void Green_2D_3dirunbounded_0dirspectral_hej2 (const int size[3],const double hfact[3],double* green,const double alpha);
void Green_2D_3dirunbounded_0dirspectral_hej4(const int size[3],const double hfact[3],double* green,const double alpha);
/**@} */


/**
 * @name 2 dir unbounded - 1 dir spectral 
 * 
 */
/**@{ */
// void Green_2D_2dirunbounded_1dirspectral_chat2_diffinf(const size_t size[3],const double hfact[3],const double kfact[3],double* green);
// void Green_2D_2dirunbounded_1dirspectral_chat2_diff2  (const size_t size[3],const double hfact[3],const double kfact[3],double* green);
/**@} */

/**
 * @name 1 dir unbounded - 2 dir spectral 
 * 
 */
/**@{ */
// void Green_2D_1dirunbounded_2dirspectral_chat2_diffinf(const size_t size[3],const double hfact[3],const double kfact[3],double* green);
// void Green_2D_1dirunbounded_2dirspectral_chat2_diff2  (const size_t size[3],const double hfact[3],const double kfact[3],double* green);
/**@} */

/**
 * @name 0 dir unbounded - 3 dir spectral 
 * 
 */
/**@{ */
// void Green_2D_0dirunbounded_3dirspectral_chat2_diffinf(const size_t size[3],const double hfact[3],const double kfact[3],double* green);
// void Green_2D_0dirunbounded_3dirspectral_chat2_diff2  (const size_t size[3],const double hfact[3],const double kfact[3],double* green);
/**@} */

#endif