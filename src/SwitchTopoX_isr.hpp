/**
 * @file SwitchTopoX_isr.hpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
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