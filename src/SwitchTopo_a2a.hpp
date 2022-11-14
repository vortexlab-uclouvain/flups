/**
 * @file SwitchTopo_a2a.hpp
 * @copyright Copyright © UCLouvain 2020
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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