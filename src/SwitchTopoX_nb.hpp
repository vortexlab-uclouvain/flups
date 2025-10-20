#ifndef SRC_SWITCHTOPOX_NB_HPP_
#define SRC_SWITCHTOPOX_NB_HPP_

#include "SwitchTopoX.hpp"

class SwitchTopoX_nb : public SwitchTopoX {
    int*     completed_id_ = nullptr;        //!< array used by the Wait/Test in the non-blocking comms
    int*     recv_order_   = nullptr;        //!< array used by the Wait/Test in the non-blocking comms
    MPI_Comm shared_comm_  = MPI_COMM_NULL;  //<! communicators with ranks on the same node

    // bool* i2o_priority_list_ = nullptr; //!< indicate if a request is in the priority list or not
    // bool* o2i_priority_list_ = nullptr; //!< indicate if a request is in the priority list or not
    int* i2o_send_order_ = nullptr;
    int* o2i_send_order_ = nullptr;

    MPI_Request* i2o_send_rqst_ = NULL;  //!< MPI send requests
    MPI_Request* i2o_recv_rqst_ = NULL;  //!< MPI recv requests
    MPI_Request* o2i_send_rqst_ = NULL;  //!< MPI send requests
    MPI_Request* o2i_recv_rqst_ = NULL;  //!< MPI recv requests

   public:
    explicit SwitchTopoX_nb(const Topology* topo_in, const Topology* topo_out, const int shift[3], H3LPR::Profiler* prof);
    ~SwitchTopoX_nb();

    virtual bool need_send_buf()const override{return true;};
    virtual bool need_recv_buf()const override{return true;};


    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) override;
    virtual void execute(opt_double_ptr data, const int sign) const override;
    virtual void disp() const override;
};

#endif