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