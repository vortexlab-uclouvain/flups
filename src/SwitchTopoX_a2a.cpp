#include "SwitchTopoX_a2a.hpp"

void All2Allv(MemChunk *send_chunks, const int *count_send, const int *disp_send,
              MemChunk *recv_chunks, const int *count_recv, const int *disp_recv, 
              opt_double_ptr send_buf, opt_double_ptr recv_buf, MPI_Comm subcomm,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof);

SwitchTopoX_a2a::SwitchTopoX_a2a(const Topology *topo_in, const Topology *topo_out, const int shift[3], H3LPR::Profiler *prof)
    : SwitchTopoX(topo_in, topo_out, shift, prof) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // nothing special to do here
    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_a2a::setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // first setup the basic stuffs
    this->SwitchTopoX::setup_buffers(sendData, recvData);

    // Retrieve MPI information for the subcomm
    int sub_rank, sub_size;
    MPI_Comm_rank(subcomm_, &sub_rank);
    MPI_Comm_size(subcomm_, &sub_size);

    //..........................................................................
    // Allocate the arrays needed by MPI_all2allv
    i2o_count_ = reinterpret_cast<int *>(m_calloc(sub_size * sizeof(int)));
    i2o_disp_  = reinterpret_cast<int *>(m_calloc(sub_size * sizeof(int)));
    o2i_count_ = reinterpret_cast<int *>(m_calloc(sub_size * sizeof(int)));
    o2i_disp_  = reinterpret_cast<int *>(m_calloc(sub_size * sizeof(int)));

    //..........................................................................
    // this is the loop over the input topo and the associated chunks
    // Chunks are organised by rank -> First chunks in the array correspond to the 
    // first rank in the sub comm. We can fill the array of displacement in the same 
    // loop than the array of count.
    FLUPS_CHECK( i2o_nchunks_ == sub_size, "Number of chunk should be equal to the number of rank in the communication %d vs %d", cchunk->dest_rank, i2o_nchunks_, sub_size);
    int cdisp = 0;
    for (int ic = 0; ic < i2o_nchunks_; ++ic){
        MemChunk *cchunk = i2o_chunks_ + ic;
        size_t count = cchunk->size_padded * cchunk->nda;
        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        FLUPS_CHECK(cchunk->dest_rank < sub_size, "Destination rank of the chunk should be inside the subcomm %d vs %d", cchunk->dest_rank, sub_size);
        FLUPS_CHECK(cchunk->dest_rank == ic, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ic);
        i2o_count_[cchunk->dest_rank] = count;
        i2o_disp_[cchunk->dest_rank]  = cdisp;
        cdisp += count;
    }

    FLUPS_CHECK( o2i_nchunks_ == sub_size, "Number of chunk should be equal to the number of rank in the communication %d vs %d", cchunk->dest_rank, i2o_nchunks_, sub_size);
    cdisp = 0;
    for (int ic = 0; ic < o2i_nchunks_; ++ic){
        MemChunk *cchunk = o2i_chunks_ + ic;
        size_t count = cchunk->size_padded * cchunk->nda;
        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        FLUPS_CHECK(cchunk->dest_rank < sub_size, "Destination rank of the chunk should be inside the subcomm %d vs %d", cchunk->dest_rank, sub_size);
        FLUPS_CHECK(cchunk->dest_rank == ic, "Destination rank of the chunk should correspond to the rank indexing %d vs %d", cchunk->dest_rank, ic);
        o2i_count_[cchunk->dest_rank] = count;
        o2i_disp_[cchunk->dest_rank]  = cdisp;
        cdisp += count;
    }
    
    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX_a2a::~SwitchTopoX_a2a(){
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // free the request arrays
    m_free(i2o_count_);
    m_free(i2o_disp_);
    m_free(o2i_count_);
    m_free(o2i_disp_);
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Send and receive the non-blocking calls, overlaping with the shuffle execution
 *
 * @param v
 * @param sign
 */
void SwitchTopoX_a2a::execute(opt_double_ptr v, const int sign) const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    m_profStarti(prof_, "Switchtopo%d", idswitchtopo_);
    if (sign == FLUPS_FORWARD) {
        All2Allv(i2o_chunks_, i2o_count_, i2o_disp_,
                 o2i_chunks_, o2i_count_, o2i_disp_,
                 send_buf_, recv_buf_, subcomm_,
                 topo_in_, topo_out_, v, prof_);
    } else {
        All2Allv(o2i_chunks_, o2i_count_, o2i_disp_,
                 i2o_chunks_, i2o_count_, i2o_disp_,
                 recv_buf_, send_buf_, subcomm_,
                 topo_out_, topo_in_, v, prof_);
    }
    m_profStopi(prof_, "Switchtopo%d", idswitchtopo_);
    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_a2a::disp() const {
    BEGIN_FUNC;
    FLUPS_INFO("------------------------------------------");
    FLUPS_INFO("## Topo Swticher MPI");
    FLUPS_INFO("--- INPUT");
    FLUPS_INFO("  - input axis = %d", topo_in_->axis());
    FLUPS_INFO("  - input local = %d %d %d", topo_in_->nloc(0), topo_in_->nloc(1), topo_in_->nloc(2));
    FLUPS_INFO("  - input global = %d %d %d", topo_in_->nglob(0), topo_in_->nglob(1), topo_in_->nglob(2));
    FLUPS_INFO("--- OUTPUT");
    FLUPS_INFO("  - output axis = %d", topo_out_->axis());
    FLUPS_INFO("  - output local = %d %d %d", topo_out_->nloc(0), topo_out_->nloc(1), topo_out_->nloc(2));
    FLUPS_INFO("  - output global = %d %d %d", topo_out_->nglob(0), topo_out_->nglob(1), topo_out_->nglob(2));
    FLUPS_INFO("------------------------------------------");
}

/**
 * @brief process to the Send/Recv operation to go from topo_in to topo_out
 *
 * The indexing of the requests and the chunks is the same
 *
 * @param n_send_rqst
 * @param send_rqst
 * @param send_chunks
 * @param n_recv_rqst
 * @param recv_rqst
 * @param recv_chunks
 * @param topo_in
 * @param topo_out
 * @param mem
 */
void All2Allv(MemChunk *send_chunks, const int *count_send, const int *disp_send,
              MemChunk *recv_chunks, const int *count_recv, const int *disp_recv, 
              opt_double_ptr send_buf, opt_double_ptr recv_buf, MPI_Comm subcomm,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    int sub_rank, sub_size; 
    MPI_Comm_rank(subcomm, &sub_rank);
    MPI_Comm_size(subcomm, &sub_size);
    //..........................................................................
    auto set_sendbuf = [=](MemChunk *chunk) {
        FLUPS_INFO("sending request to rank %d of size %d %d %d",chunk->dest_rank,chunk->isize[0],chunk->isize[1],chunk->isize[2]);
        // copy here the chunk from the input topo to the chunk
        m_profStarti(prof, "Copy data 2 chunk");
        CopyData2Chunk(nmem_in, mem, chunk);
        m_profStopi(prof, "Copy data 2 chunk");
    };

    auto complete_recv = [=](MemChunk *chunk) {
        FLUPS_INFO("recving request from rank %d of size %d %d %d",chunk->dest_rank,chunk->isize[0],chunk->isize[1],chunk->isize[2]);
        // shuffle the data
        m_profStarti(prof, "Do Shuffle Chunk");
        DoShuffleChunk(chunk);
        m_profStopi(prof, "Do Shuffle Chunk");
        m_profStarti(prof, "Copy Chunk 2 data ");
        CopyChunk2Data(chunk, nmem_out, mem);
        m_profStopi(prof, "Copy Chunk 2 data ");
    };
    //..........................................................................
    // Prepare the send buffer
    for(int ir = 0; ir < sub_size; ++ir){
        MemChunk *cchunk = send_chunks + ir;
        FLUPS_CHECK(cchunk->dest_rank == ir, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ir);
        set_sendbuf(cchunk);
    }

    m_profStarti(prof, "All2All");
    MPI_Alltoallv(send_buf, count_send, disp_send, MPI_DOUBLE, recv_buf, count_recv, disp_recv, MPI_DOUBLE, subcomm);
    m_profStopi(prof, "All2All");

    // Copy back the recveived data 
    for(int ir = 0; ir < sub_size; ++ir){
        MemChunk *cchunk = recv_chunks + ir;
        FLUPS_CHECK(cchunk->dest_rank == ir, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ir);
        complete_recv(cchunk);
    }

    //--------------------------------------------------------------------------
    END_FUNC;
}
