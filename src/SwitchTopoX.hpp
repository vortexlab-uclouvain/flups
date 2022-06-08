#ifndef SWITCHTOPOX_HPP_
#define SWITCHTOPOX_HPP_

#include "Topology.hpp"
#include "chunk_tools.hpp"
#include "defines.hpp"

/**
 * @brief More efficient implementation of the SwitchTopo
 *
 */
class SwitchTopoX {
   protected:
    const int i2o_shift_[3];  //!< position of the (0,0,0) of topo_in_ in topo_out_
    const int o2i_shift_[3];  //!< position of the (0,0,0) of topo_out_ in topo_in_

    const Topology *topo_in_  = NULL; /**<@brief input topology  */
    const Topology *topo_out_ = NULL; /**<@brief  output topology */

    MPI_Comm inComm_  = MPI_COMM_NULL; /**<@brief the reference input communicator */
    MPI_Comm outComm_ = MPI_COMM_NULL; /**<@brief the reference output communicator */
    MPI_Comm subcomm_ = MPI_COMM_NULL; /**<@brief the subcomm for this switchTopo */

    int i2o_nchunks_ = 0;  //!< local number of chunks in the input topology
    int o2i_nchunks_ = 0;  //!< local number of chunks in the output topology

    MemChunk *i2o_chunks_ = NULL;  //!< the local chunks of memory in the output topology
    MemChunk *o2i_chunks_ = NULL;  //!< the local chunks of memory in the output topology

    int i2o_selfcomm_ = -1; //!< Index of the self communication chunk (remains at -1 if there is no self communication)
    int o2i_selfcomm_ = -1; //!< Index of the self communication chunk (remains at -1 if there is no self communication)

    opt_double_ptr send_buf_ = NULL; /**<@brief The send buffer for MPI send */
    opt_double_ptr recv_buf_ = NULL; /**<@brief The recv buffer for MPI recv */

    fftw_plan *i2o_shuffle_ = NULL;  //!< FFTW plan to shuffle the indexes around from the input topo to the output topo
    fftw_plan *o2i_shuffle_ = NULL;  //!< FFTW plan to shuffle the indexes around from the input topo to the ouput topo

    H3LPR::Profiler *prof_         = NULL;
    int              idswitchtopo_ = -1;

   public:
    explicit SwitchTopoX(const Topology *topo_in, const Topology *topo_out, const int shift[3], H3LPR::Profiler *prof);
    virtual ~SwitchTopoX();

    virtual bool need_send_buf()const  = 0;
    virtual bool need_recv_buf()const  = 0;

    // abstract functions
    void setup();
    virtual void print_info() const;
    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData);
    virtual void execute(opt_double_ptr data, const int sign) const = 0;
    virtual void disp() const                                       = 0;
    

    size_t get_bufMemSize() const;
    size_t get_ChunkArraysMemSize(const int lda, const int nchunks, const MemChunk *chunks) const;

   protected:
    void SubCom_SplitComm();
    // void SubCom_UpdateRanks();
    // setup_subComm_(const int nBlock, const int lda, int *blockSize[3], int *destRank, int **count, int **start);
};

#endif  // SWITCHTOPOX_HPP_