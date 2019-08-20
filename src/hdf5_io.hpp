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
#include "topology.hpp"

void hdf5_write(Topology* topo,const double* data);

void xmf_write(Topology* topo);


#endif