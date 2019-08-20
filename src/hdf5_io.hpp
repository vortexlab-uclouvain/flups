/**
 * @file hdf5_write.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-08-20
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef HDF5_IO_HPP
#define HDF5_IO_HPP

#include "hdf5.h"
#include "defines.hpp"
#include "topology.hpp"

using namespace std;

void hdf5_write(const Topology *topo, const string filename, const string attribute, const double *data);

void xmf_write(const Topology *topo, const string filename, const string attribute);

#endif