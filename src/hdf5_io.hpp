/**
 * @file hdf5_io.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef HDF5_IO_HPP
#define HDF5_IO_HPP



#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>


#include "defines.hpp"
#include "Topology.hpp"

#if (FLUPS_HDF5)
#include "hdf5.h"
#endif

void hdf5_dumptest();
void hdf5_dump(const Topology *topo, const std::string filename, const double *data);

void xmf_write(const Topology *topo, const std::string filename, const std::string attribute);
void hdf5_write(const Topology *topo, const std::string filename, const std::string attribute, const double *data);
#endif
