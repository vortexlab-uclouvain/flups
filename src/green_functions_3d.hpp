/**
 * @file green_functions_3d.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

#include "defines.hpp"
#include "Topology.hpp"
#include "bessel.hpp"

void cmpt_Green_3D_3dirunbounded_0dirspectral(const FLUPS::Topology *topo, const double hfact[3],                                                 const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double alpha);
void cmpt_Green_3D_2dirunbounded_1dirspectral(const FLUPS::Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double alpha);
void cmpt_Green_3D_1dirunbounded_2dirspectral(const FLUPS::Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double alpha);
void cmpt_Green_3D_0dirunbounded_3dirspectral(const FLUPS::Topology *topo,                        const double kfact[3], const double koffset[3], const double symstart[3], double *green, FLUPS::GreenType typeGreen, const double alpha);