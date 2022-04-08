#ifndef SRC_SWITCHTOPOX_NB_HPP_
#define SRC_SWITCHTOPOX_NB_HPP_

#include "SwitchTopoX.hpp"

class SwitchTopoX_nb : public SwitchTopoX {
    MPI_Request* i2o_send_rqst_ = NULL;  //!< MPI send requests
    MPI_Request* i2o_recv_rqst_ = NULL;  //!< MPI recv requests
    MPI_Request* o2i_send_rqst_ = NULL;  //!< MPI send requests
    MPI_Request* o2i_recv_rqst_ = NULL;  //!< MPI recv requests

   public:
    explicit SwitchTopoX_nb(const int shift[3], Topology* topo_in, Topology* topo_out, H3LPR::Profiler* prof);
    ~SwitchTopoX_nb();

    virtual void setup();
    virtual void execute(opt_double_ptr data, const int sign) const override;
    virtual void disp() const override;
};

#endif