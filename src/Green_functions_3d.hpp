/**
 * @file Green_functions_3d.hpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
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
enum GreenType {
    CHAT_2, /**< @brief quadrature in zero, order 2, Chatelain et al. (2010) */
    LGF_2,  /**< @brief Lattice Green's function, order 2, Gillis et al. (2018)*/
    HEJ_2,  /**< @brief regularized in zero, order 2, Hejlesen et al. (2015)*/
    HEJ_4,  /**< @brief regularized in zero, order 4, Hejlesen et al. (2015)*/
    HEJ_6,  /**< @brief regularized in zero, order 6, Hejlesen et al. (2015)*/
};

void cmpt_Green_3D_3dirunbounded_0dirspectral(const Topology *topo, const double hfact[3], const int symstart[3], double *green, GreenType typeGreen, const double alpha);
void cmpt_Green_3D_2dirunbounded_1dirspectral(const Topology *topo, const double hfact[3], const int symstart[3], double *green, GreenType typeGreen, const double alpha);
void cmpt_Green_3D_1dirunbounded_2dirspectral(const Topology *topo, const double hfact[3], const int symstart[3], double *green, GreenType typeGreen, const double alpha);
void cmpt_Green_3D_0dirunbounded_3dirspectral(const Topology *topo, const double hfact[3], const int symstart[3], double *green, GreenType typeGreen, const double alpha);