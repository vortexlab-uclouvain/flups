/**
 * @file green_functions_3d.hpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "defines.hpp"
#include "Topology.hpp"
#include "bessel.hpp"
#include "expint.hpp"

void cmpt_Green_3D_3dirunbounded_0dirspectral(const FLUPS::Topology *topo, const double hfact[3],                                                 const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double eps);
void cmpt_Green_3D_2dirunbounded_1dirspectral(const FLUPS::Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double eps);
void cmpt_Green_3D_1dirunbounded_2dirspectral(const FLUPS::Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double eps);
void cmpt_Green_3D_0dirunbounded_3dirspectral(const FLUPS::Topology *topo,                        const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double eps);
void cmpt_Green_3D_0dirunbounded_3dirspectral(const FLUPS::Topology *topo,                        const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double eps, const int istart0[3], const int ishift[3]);