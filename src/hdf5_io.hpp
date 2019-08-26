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

#include "defines.hpp"
#include "hdf5.h"
#include "topology.hpp"

using namespace std;

void hdf5_dumptest();
void hdf5_dump(const Topology *topo, const string filename, const double *data);

void xmf_write(const Topology *topo, const string filename, const string attribute);
void hdf5_write(const Topology *topo, const string filename, const string attribute, const double *data);

#endif