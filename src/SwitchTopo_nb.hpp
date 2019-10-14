/**
 * @file SwitchTopo_nb.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-28
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

#ifndef SWITCHTOPO_NB_HPP
#define SWITCHTOPO_NB_HPP

#include <cstring>
#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"
#include "Topology.hpp"
#include "Profiler.hpp"
#include "omp.h"
#include "SwitchTopo.hpp"

/**
 * @brief Takes care of the switch between to different topologies
 * 
 * Reorganize the memory between 2 different topologies, also accounting for a
 * "principal axis" which is aligned with the fast rotating index.
 * 
 * Communications are handled by packing data in blocks (i.e. chuncks). These 
 * smaller pieces are sent and recieved asynchronesouly, in a hope for reducing
 * the time required for communication.
 * 
 * The a memory shift can be specified in the switch between the topologies, when
 * there is a need for skipping some data at the left or right side of a given direction.
 * 
 */
class FLUPS::SwitchTopo_nb : public SwitchTopo {
   protected:
    int _selfBlockN=0;
    int* _iselfBlockID = NULL;
    int* _oselfBlockID = NULL;

    opt_int_ptr _i2o_destTag = NULL; /**<@brief The destination rank in the output topo of each block */
    opt_int_ptr _o2i_destTag = NULL; /**<@brief The destination rank in the output topo of each block */

    MPI_Request *_i2o_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_i2o_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */
    MPI_Request *_o2i_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_o2i_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */

   public:
    SwitchTopo_nb(const Topology *topo_input, const Topology *topo_output, const int shift[3],Profiler* prof);
    ~SwitchTopo_nb();

    void setup_buffers(opt_double_ptr _sendBuf,opt_double_ptr _recvBuf);
    void execute(double* v, const int sign) const;

    void disp() const;
};

void SwitchTopo_test();

#endif