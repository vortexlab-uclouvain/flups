#include "SwitchTopoX_nb.hpp"

void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof);

SwitchTopoX_nb::SwitchTopoX_nb(const Topology *topo_in, const Topology *topo_out, const int shift[3], H3LPR::Profiler *prof)
    : SwitchTopoX(topo_in, topo_out, shift, prof) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // nothing special to do here
    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_nb::setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // first setup the basic stuffs
    this->SwitchTopoX::setup_buffers(sendData, recvData);

    int sub_rank;
    MPI_Comm_rank(subcomm_, &sub_rank);

    //..........................................................................
    // once we have the MemChunks we can allocate the requests
    // i2o transfert goes from i2o_chunks to o2i_chunks
    i2o_send_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Request)));
    i2o_recv_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Request)));
    // o2i transfert goes from o2i_chunks to i2o_chunks
    o2i_send_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Request)));
    o2i_recv_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Request)));

    //..........................................................................
    // this is the loop over the input topo and the associated chunks
    for (int ir = 0; ir < i2o_nchunks_; ++ir) {
        MemChunk      *cchunk = i2o_chunks_ + ir;
        opt_double_ptr buf    = cchunk->data;
        size_t         count  = cchunk->size_padded * cchunk->nda;

        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        MPI_Send_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, sub_rank, subcomm_, i2o_send_rqst_ + ir);
        MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, sub_rank, subcomm_, o2i_recv_rqst_ + ir);
    }

    // here we go for the output topo and the associated chunks
    for (int ir = 0; ir < o2i_nchunks_; ++ir) {
        MemChunk      *cchunk = o2i_chunks_ + ir;
        opt_double_ptr buf    = cchunk->data;
        size_t         count  = cchunk->size_padded * cchunk->nda;
        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, cchunk->dest_rank, subcomm_, i2o_recv_rqst_ + ir);
        MPI_Send_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, cchunk->dest_rank, subcomm_, o2i_send_rqst_ + ir);
    }

    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX_nb::~SwitchTopoX_nb(){
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    for (int ir = 0; ir < i2o_nchunks_; ++ir) {
        MPI_Request_free(i2o_send_rqst_ + ir);
        MPI_Request_free(o2i_recv_rqst_ + ir);
    }
    // here we go for the output topo and the associated chunks
    for (int ir = 0; ir < o2i_nchunks_; ++ir) {
        MPI_Request_free(i2o_recv_rqst_ + ir);
        MPI_Request_free(o2i_send_rqst_ + ir);
    }

    // free the request arrays
    m_free(i2o_send_rqst_);
    m_free(i2o_recv_rqst_);
    m_free(o2i_send_rqst_);
    m_free(o2i_recv_rqst_);
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Send and receive the non-blocking calls, overlaping with the shuffle execution
 *
 * @param v
 * @param sign
 */
void SwitchTopoX_nb::execute(opt_double_ptr v, const int sign) const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    m_profStarti(prof_, "Switchtopo%d", idswitchtopo_);
    if (sign == FLUPS_FORWARD) {
        SendRecv(i2o_nchunks_, i2o_send_rqst_, i2o_chunks_,
                 o2i_nchunks_, i2o_recv_rqst_, o2i_chunks_,
                 topo_in_, topo_out_, v, prof_);
    } else {
        SendRecv(o2i_nchunks_, o2i_send_rqst_, o2i_chunks_,
                 i2o_nchunks_, o2i_recv_rqst_, i2o_chunks_,
                 topo_out_, topo_in_, v, prof_);
    }
    m_profStopi(prof_, "Switchtopo%d", idswitchtopo_);

    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_nb::disp() const {
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
void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    //..........................................................................
    auto send_my_rqst = [=]( int count, MemChunk *chunk, MPI_Request *request) {
        FLUPS_INFO("sending request to rank %d of size %d %d %d",chunk->dest_rank,chunk->isize[0],chunk->isize[1],chunk->isize[2]);
        // copy here the chunk from the input topo to the chunk
        for(int ic = 0; ic < count; ++ic){
            CopyData2Chunk(nmem_in, mem, chunk + ic);
        }
        // start the request
        MPI_Startall(count, request);
    };
    auto recv_my_rqst = [=](MPI_Request *request, MemChunk *chunk) {
        FLUPS_INFO("recving request from rank %d of size %d %d %d",chunk->dest_rank,chunk->isize[0],chunk->isize[1],chunk->isize[2]);
        // shuffle the data
        DoShuffleChunk(chunk);
        // CopyChunk2Data(chunk, nmem_out, mem);
    };
    //..........................................................................
    const int send_batch = MPI_BATCH_SEND;
    int send_cntr = 0;
    int recv_cntr = 0;
    int *completed_id = reinterpret_cast<int *>(m_calloc(n_recv_rqst * sizeof(int)));

    FLUPS_INFO("starting %d recv request",n_recv_rqst);
    // start all the recv request
    MPI_Startall(n_recv_rqst, recv_rqst);

    // while we still have to send or recv something, we continue
    m_profStarti(prof, "send/recv");
    while ((send_cntr < n_send_rqst) || (recv_cntr < n_recv_rqst)) {
        // perform a batch of  sends
        int count_send = m_min(n_send_rqst - send_cntr, send_batch);
        FLUPS_INFO("sending request from %d to %d /%d", send_cntr, send_cntr + count_send, n_send_rqst);
        m_profStarti(prof, "copy data 2 chunk");
        send_my_rqst(count_send, send_chunks + send_cntr, send_rqst + send_cntr);
        m_profStopi(prof, "copy data 2 chunk");
        // increment the counter
        send_cntr += count_send;
        
        // if we have some requests to recv, test it
        int n_completed = 0;
        if (recv_cntr < n_recv_rqst) {
            MPI_Testsome(n_recv_rqst, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);
        }

        // Maintain the difference between the pending send and receive request to send_bacth
        count_send = m_min(n_send_rqst - send_cntr, n_completed);
        FLUPS_INFO("sending request from %d to %d /%d", send_cntr, send_cntr + count_send, n_send_rqst);
        m_profStarti(prof, "copy data 2 chunk");
        send_my_rqst(count_send, send_chunks + send_cntr, send_rqst + send_cntr);
        m_profStopi(prof, "copy data 2 chunk");   
        // increment the counter
        send_cntr += count_send;
        
        // for each of the completed request, treat it
        for (int id = 0; id < n_completed; ++id) {
            m_profStarti(prof, "shuffle");
            FLUPS_INFO("recving request %d/%d", recv_cntr + id, n_recv_rqst);
            const int rqst_id = completed_id[id];
            recv_my_rqst(recv_rqst + rqst_id, recv_chunks + rqst_id);
            m_profStopi(prof, "shuffle");
        }
    
        // increment the counter of recv request
        recv_cntr += n_completed;
    

    }
    m_profStopi(prof, "send/recv");

    // we need to officially close the send requests
    m_profStarti(prof, "Wait all");
    MPI_Waitall(n_send_rqst, send_rqst, MPI_STATUSES_IGNORE);
    m_profStopi(prof, "Wait all");

    // reset the memory to 0.0 as we do inplace computations
    const size_t reset_size = topo_out->memsize();
    std::memset(mem, 0, reset_size * sizeof(double));

    // once all the send has been done we can overwrite the received info
    m_profStarti(prof, "copy chunk 2 data");
    for (int ic = 0; ic < n_recv_rqst; ++ic) {
        CopyChunk2Data(recv_chunks + ic, nmem_out, mem);
    }
    m_profStopi(prof, "copy chunk 2 data");

    // m_free(completed_status);
    m_free(completed_id);
    //--------------------------------------------------------------------------
    END_FUNC;
}
