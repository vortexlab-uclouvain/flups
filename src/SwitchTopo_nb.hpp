/**
 * @file SwitchTopo_nb.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * 
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
 * @brief Switch between to different topologies using non-blocking communications
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
class SwitchTopo_nb : public SwitchTopo {
   protected:
    int _selfBlockN=0;
    int* _iselfBlockID = NULL; /**<@brief The list of the block iD that stays on the current rank in the input topology (used while output to input)  */
    int* _oselfBlockID = NULL; /**<@brief The list of the block iD that stays on the current rank in the output topoloy (used while input to ouput) */

    int* _i2o_destTag = NULL; /**<@brief The destination rank in the output topo of each block */
    int* _o2i_destTag = NULL; /**<@brief The destination rank in the output topo of each block */

    MPI_Request *_i2o_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_i2o_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */
    MPI_Request *_o2i_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_o2i_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */

    void _init_blockInfo(const Topology* topo_in, const Topology* topo_out);
    void _free_blockInfo();

   public:
    SwitchTopo_nb(const Topology *topo_input, const Topology *topo_output, const int shift[3],Profiler* prof);
    ~SwitchTopo_nb();

    void setup_buffers(opt_double_ptr _sendBuf,opt_double_ptr _recvBuf);
    void execute(double* v, const int sign) const;
    void setup() ;
    void disp() const;
};

void SwitchTopo_test();

#endif