/**
 * @file SwitchTopoX_isr.hpp
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
#ifndef SRC_SWITCHTOPOX_ISR_HPP_
#define SRC_SWITCHTOPOX_ISR_HPP_

#include "SwitchTopoX.hpp"

class SwitchTopoX_isr : public SwitchTopoX {
    int* i2o_send_order_ = nullptr;  //!< order in which to perform the send for input to output
    int* o2i_send_order_ = nullptr;  //!< order in which to perform the send for output to input
    int* completed_id_   = nullptr;  //!< array used by the Wait/Test in the non-blocking comms
    int* recv_order_     = nullptr;  //!< array used by the Wait/Test in the non-blocking comms

    MPI_Comm shared_comm_ = MPI_COMM_NULL;  //<! communicators with ranks on the same node

    MPI_Request* send_rqst_;  //<! storage for send requests
    MPI_Request* recv_rqst_;  //<! storage for recv requests

   public:
    explicit SwitchTopoX_isr(const Topology* topo_in, const Topology* topo_out, const int shift[3], H3LPR::Profiler* prof);
    ~SwitchTopoX_isr();

    virtual bool need_send_buf() const override { return false; };
    virtual bool need_recv_buf() const override { return true; };

    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) override;
    virtual void execute(opt_double_ptr data, const int sign) const override;
    virtual void disp() const override;
};

#endif