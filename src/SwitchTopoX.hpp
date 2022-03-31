#ifndef SWITCHTOPOX_HPP_
#define SWITCHTOPOX_HPP_

#include <mpi.h>
#include "Topology.hpp"
#include "defines.hpp"

typedef struct
{
    const int istart[3];  //!< local start index for the memory chunk (012 indexing)
    const int isize[3];   // !< local size for the memory chunks (012 indexing)

    const int dest_rank;  //!< destination ranks
    const int tag;        //!< tag to be used

} MemChunk;

/**
 * @brief More efficient implementation of the SwitchTopo
 *
 */
class SwitchTopoX {
   protected:
    MPI_Comm inComm_  = NULL; /**<@brief the reference input communicator */
    MPI_Comm outComm_ = NULL; /**<@brief the reference output communicator */
    MPI_Comm subcomm_ = NULL; /**<@brief the subcomm for this switchTopo */

    const int shift[3];  //!< shift from global to local indexes

    const int inchunks = 0;  //!< local number of chunks in the input topology
    const int onchunks = 0;  //!< local number of chunks in the output topology

    MemChunk *ichunks;  //!< the local chunks of memory in the output topology
    MemChunk *ochunks;  //!< the local chunks of memory in the output topology

    const Topology *topo_in_  = NULL; /**<@brief input topology  */
    const Topology *topo_out_ = NULL; /**<@brief  output topology */

    opt_double_ptr *sendBuf_ = NULL; /**<@brief The send buffer for MPI send */
    opt_double_ptr *recvBuf_ = NULL; /**<@brief The recv buffer for MPI recv */

    fftw_plan *i2o_shuffle_ = NULL;  //!< FFTW plan to shuffle the indexes around from the input topo to the output topo
    fftw_plan *o2i_shuffle_ = NULL;  //!< FFTW plan to shuffle the indexes around from the input topo to the ouput topo

#ifdef PROF
    Profiler *prof_    = NULL;
    int       iswitch_ = -1;
#endif


    public:
    explicit SwitchTopoX();
    virtual ~SwitchTopoX(){};

    virtual void setup()                                                                    = 0;
    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData)            = 0;
    virtual void execute(opt_double_ptr v, const int sign) const                            = 0;
    virtual void disp() const                                                               = 0;



    /**
     * @brief return the buffer size for one proc = number of blocks * blocks memory size * lda component
     * 
     * @return size_t 
     */
    inline size_t get_bufMemSize() const {
        // the nf is the maximum between in and out
        const int nf = std::max(topo_in_->nf(),topo_out_->nf());
        // nultiply by the number of blocks
        size_t total = 0;
        for(int ib=0; ib<ichunks; ib++){
            total += get_blockMemSize(ib,nf,iBlockSize_) * ((size_t)topo_in_->lda());
        }
        for(int ib=0; ib<onBlock_; ib++){
            total += get_blockMemSize(ib,nf,oBlockSize_) * ((size_t)topo_in_->lda());
        }
        // return the total size
        return total;
    };

   protected:
    
    /**
     * @brief returns the memory size (!NOT padded!) of a MemChunk
     * 
     * @param nf 
     * @param chunk 
     * @return size_t 
     */
    inline size_t get_ChunkMemSize(const int nf, const MemChunk *chunk) {
        //----------------------------------------------------------------------
        const size_t total      = (size_t)(chunk->isize[0]) * (size_t)(chunk->isize[1]) * (size_t)(chunk->isize[2]) * (size_t)(nf);
        // size_t alignDelta = ((total * sizeof(double)) % FLUPS_ALIGNMENT == 0) ? 0 : (FLUPS_ALIGNMENT - (total * sizeof(double)) % FLUPS_ALIGNMENT) / sizeof(double);
        // return total + alignDelta;
        return total;
        //----------------------------------------------------------------------
    }
}

#endif  // SWITCHTOPOX_HPP_