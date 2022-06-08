#ifndef SRC_SWITCHTOPOX_ISR_HPP_
#define SRC_SWITCHTOPOX_ISR_HPP_

#include "SwitchTopoX.hpp"

class SwitchTopoX_isr : public SwitchTopoX {
    int*     completed_id_ = nullptr;        //!< array used by the Wait/Test in the non-blocking comms
    MPI_Comm shared_comm_  = MPI_COMM_NULL;  //<! communicators with ranks on the same node

    int* i2o_send_order_ = nullptr;
    int* o2i_send_order_ = nullptr;

    size_t*       i2o_offset_;
    size_t*       o2i_offset_;

    MPI_Datatype* i2o_dtype_;
    MPI_Datatype* o2i_dtype_;

    MPI_Request* send_rqst_; //<! storage for send requests
    MPI_Request* recv_rqst_; //<! storage for recv requests

   public:
    explicit SwitchTopoX_isr(const Topology* topo_in, const Topology* topo_out, const int shift[3], H3LPR::Profiler* prof);
    ~SwitchTopoX_isr();

    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) override;
    virtual void execute(opt_double_ptr data, const int sign) const override;
    virtual void disp() const override;
};

#endif