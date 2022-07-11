/**
 * @file SwitchTopo_a2a.hpp
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


#ifndef SWITCHTOPO_A2A_HPP
#define SWITCHTOPO_A2A_HPP

#include <cstring>
#include "mpi.h"


#include "defines.hpp"
#include "Topology.hpp"
#include "SwitchTopo.hpp"
#include "hdf5_io.hpp"


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
class SwitchTopo_a2a : public SwitchTopo {
   protected:
    bool is_all2all_ = false; /**<@brief is the call an alltoall or an alltoall_v */

    int *i2o_count_ = NULL; /**<@brief count argument of the all_to_all_v for input to output */
    int *o2i_count_ = NULL; /**<@brief count argument of the all_to_all_v for output to input */
    int *i2o_start_ = NULL; /**<@brief start argument of the all_to_all_v for input to output */
    int *o2i_start_ = NULL; /**<@brief start argument of the all_to_all_v for output to input */

    opt_double_ptr sendBufG_ = NULL; /**<@brief pointer to the globally allocated memory for the send buffers */
    opt_double_ptr recvBufG_ = NULL; /**<@brief pointer to the globally allocated memory for the recv buffers */

    void init_blockInfo_(const Topology* topo_in, const Topology* topo_out);
    void free_blockInfo_();

   public:
    SwitchTopo_a2a(const Topology *topo_input, const Topology *topo_output, const int shift[3], H3LPR::Profiler *prof);
    ~SwitchTopo_a2a();

    void setup_buffers(opt_double_ptr sendBuf, opt_double_ptr recvBuf) ;
    void execute(double* v, const int sign) const;
    void setup();
    void disp() const;
};

void SwitchTopo_a2a_test();
void SwitchTopo_a2a_test2();

#endif