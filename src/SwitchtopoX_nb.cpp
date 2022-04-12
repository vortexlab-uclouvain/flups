#include "SwitchTopoX_nb.hpp"

void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem);

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
    if (sign == FLUPS_FORWARD) {
        SendRecv(i2o_nchunks_, i2o_send_rqst_, i2o_chunks_,
                 o2i_nchunks_, i2o_recv_rqst_, o2i_chunks_,
                 topo_in_, topo_out_, v);
    } else {
        SendRecv(o2i_nchunks_, o2i_send_rqst_, o2i_chunks_,
                 o2i_nchunks_, o2i_recv_rqst_, i2o_chunks_,
                 topo_out_, topo_in_, v);
    }

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
    // FLUPS_INFO("  - istart = %d %d %d", istart_[0], istart_[1], istart_[2]);
    // FLUPS_INFO("  - iend = %d %d %d", iend_[0], iend_[1], iend_[2]);
    FLUPS_INFO("--- OUTPUT");
    FLUPS_INFO("  - output axis = %d", topo_out_->axis());
    FLUPS_INFO("  - output local = %d %d %d", topo_out_->nloc(0), topo_out_->nloc(1), topo_out_->nloc(2));
    FLUPS_INFO("  - output global = %d %d %d", topo_out_->nglob(0), topo_out_->nglob(1), topo_out_->nglob(2));
    // FLUPS_INFO("  - ostart = %d %d %d", ostart_[0], ostart_[1], ostart_[2]);
    // FLUPS_INFO("  - oend = %d %d %d", oend_[0], oend_[1], oend_[2]);
    FLUPS_INFO("--- Chunks");
    // FLUPS_INFO("  - nByBlock  = %d %d %d", nByBlock_[0], nByBlock_[1], nByBlock_[2]);
    FLUPS_INFO("  - i2o n chunks = %d", i2o_nchunks_);
    FLUPS_INFO("  - o2i n chunks = %d", o2i_nchunks_);
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
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    //..........................................................................
    auto send_my_rqst = [=](MemChunk *chunk, MPI_Request *request) {
        // copy here the chunk from the input topo to the chunk
        CopyData2Chunk(nmem_in, mem, chunk);
        // start the request
        MPI_Start(request);
    };
    auto recv_my_rqst = [=](MPI_Request *request, MemChunk *chunk) {
        // shuffle the data
        DoShuffleChunk(chunk);
        CopyChunk2Data(chunk, nmem_out, mem);
    };
    //..........................................................................
    int         send_cntr        = 0;
    int         recv_cntr        = 0;
    int        *completed_id     = reinterpret_cast<int *>(m_calloc(n_recv_rqst * sizeof(int)));
    // MPI_Status *completed_status = reinterpret_cast<MPI_Status *>(m_calloc(n_recv_rqst * sizeof(MPI_Status)));

    FLUPS_INFO("starting %d recv request",n_recv_rqst);
    // start all the recv request
    MPI_Startall(n_recv_rqst, recv_rqst);

    // while we still have to send or recv something, we continue
    while (send_cntr < n_send_rqst && recv_cntr < n_recv_rqst) {
        // if we have some requests to recv, test it
        if (recv_cntr < n_recv_rqst) {
            int n_completed;
            MPI_Testsome(n_recv_rqst, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);

            // for each of the completed request, treat it
            for (int id = 0; id < n_completed; ++id) {
                const int rqst_id = completed_id[id];
                recv_my_rqst(recv_rqst + rqst_id, recv_chunks + rqst_id);
            }
            // increment the counter
            recv_cntr += n_completed;
        }

        // perform one of the send
        if (send_cntr < n_send_rqst) {
            send_my_rqst(send_chunks + send_cntr, send_rqst + send_cntr);
            // increment the counter
            send_cntr += 1;
        }
    }

    // m_free(completed_status);
    m_free(completed_id);
    //--------------------------------------------------------------------------
    END_FUNC;
}
