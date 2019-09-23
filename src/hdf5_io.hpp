/**
 * @file hdf5_io.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
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

#ifndef HDF5_IO_HPP
#define HDF5_IO_HPP

#include "defines.hpp"
#include "hdf5.h"
#include "Topology.hpp"

using namespace std;

void hdf5_dumptest();
void hdf5_dump(const FLUPS::Topology *topo, const string filename, const double *data);

void xmf_write(const FLUPS::Topology *topo, const string filename, const string attribute);
void hdf5_write(const FLUPS::Topology *topo, const string filename, const string attribute, const double *data);

#endif