/**
 * @file green_functions_3d.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "defines.hpp"

/**
 * @name 3 dir unbounded - 0 dir spectral 
 * 
 */
/**@{ */
void Green_3D_3dirunbounded_0dirspectral_chat2(const int size[3],const double hfact[3],double* green);
// void Green_2D_3dirunbounded_0dirspectral_lgf2 (const size_t size[3],const double hfact[3],double* green);
void Green_3D_3dirunbounded_0dirspectral_hej2 (const int size[3],const double hfact[3],double* green,const double alpha);
void Green_3D_3dirunbounded_0dirspectral_hej4(const int size[3],const double hfact[3],double* green,const double alpha);
/**@} */
