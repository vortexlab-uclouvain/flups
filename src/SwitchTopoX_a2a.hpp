/**
 * @file SwitchTopoX_a2a.hpp
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

#ifndef SRC_SWITCHTOPOX_A2A_HPP_
#define SRC_SWITCHTOPOX_A2A_HPP_

#include "SwitchTopoX.hpp"

class SwitchTopoX_a2a : public SwitchTopoX {
    MPI_Request* i2o_rqst_ = NULL;  //!< MPI i2o requests
    MPI_Request* o2i_rqst_ = NULL;  //!< MPI o2i requests

    int *i2o_count_ = NULL; /**<@brief count argument of the all_to_all_v for input to output */
    int *i2o_disp_  = NULL; /**<@brief start argument of the all_to_all_v for input to output */
    int *o2i_count_ = NULL; /**<@brief count argument of the all_to_all_v for output to input */
    int *o2i_disp_  = NULL; /**<@brief start argument of the all_to_all_v for output to input */

   public:
    explicit SwitchTopoX_a2a(const Topology* topo_in, const Topology* topo_out, const int shift[3], H3LPR::Profiler* prof);
    ~SwitchTopoX_a2a();

    void print_info() const override;

    virtual bool need_send_buf()const override{return true;};
    virtual bool need_recv_buf()const override{return true;};

    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) override;
    virtual void execute(opt_double_ptr data, const int sign) const override;
    virtual void disp() const override;
};

#endif