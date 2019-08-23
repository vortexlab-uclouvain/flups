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
#include "topology.hpp"

/**
 * @brief The type of Green's function used for the Poisson solver
 * 
 */
enum GreenType{
    CHAT_2, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2,  /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2,  /**< @brief regularized in zero, order 2, Hejlesen et al. (2015)*/
    HEJ_4,  /**< @brief regularized in zero, order 4, Hejlesen et al. (2015)*/
    HEJ_6,  /**< @brief regularized in zero, order 6, Hejlesen et al. (2015)*/
};


/**
 * @name 3 dir unbounded - 0 dir spectral 
 * 
 */
/**@{ */
// void Green_3D_3dirunbounded_0dirspectral_chat2(const Topology* topo,const double hfact[3],double* green);
// void Green_2D_3dirunbounded_0dirspectral_lgf2 (const size_t size[3],const double hfact[3],double* green);
void Green_3D_3dirunbounded_0dirspectral_hej2(const Topology *topo, const double hfact[3], const int symstart[3], double *green, const double alpha);
void Green_3D_3dirunbounded_0dirspectral_hej4(const Topology *topo, const double hfact[3], const int symstart[3], double *green, const double alpha);
/**@} */
