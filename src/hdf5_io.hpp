/**
 * @file hdf5_io.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef HDF5_IO_HPP
#define HDF5_IO_HPP

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>

#include "defines.hpp"
#include "hdf5.h"
#include "Topology.hpp"

using namespace std;

void hdf5_dumptest();
void hdf5_dump(const Topology *topo, const string filename, const double *data);

void xmf_write(const Topology *topo, const string filename, const string attribute);
void hdf5_write(const Topology *topo, const string filename, const string attribute, const double *data);

#endif