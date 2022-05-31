#include "SwitchTopoX_nb.hpp"

void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks, const int self_send_idx,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks, const int self_recv_idx,
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
    FLUPS_CHECK(i2o_selfcomm_ * o2i_selfcomm_ >= 0, "The selfcomm i2o and o2i must have the same sign ");
    //--------------------------------------------------------------------------
    // first setup the basic stuffs
    this->SwitchTopoX::setup_buffers(sendData, recvData);

    int sub_rank;
    MPI_Comm_rank(subcomm_, &sub_rank);

    //..........................................................................
    // once we have the MemChunks we can allocate the requests
    // i2o transfert goes from i2o_chunks to o2i_chunks
    // We do not use MPI to process the self request so we have to remove it from the number of request
    i2o_send_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Request)));
    i2o_recv_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Request)));
    // o2i transfert goes from o2i_chunks to i2o_chunks
    o2i_send_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Request)));
    o2i_recv_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Request)));

    //..........................................................................
    // this is the loop over the input topo and the associated chunks
    for (int ir = 0; ir < i2o_nchunks_ - (i2o_selfcomm_ >= 0 ); ++ir) {
        MemChunk      *cchunk = i2o_chunks_ + ir;
        opt_double_ptr buf    = cchunk->data;
        size_t         count  = cchunk->size_padded * cchunk->nda;

        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        MPI_Send_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, sub_rank, subcomm_, i2o_send_rqst_ + ir);
        MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, sub_rank, subcomm_, o2i_recv_rqst_ + ir);
    }

    // here we go for the output topo and the associated chunks
    for (int ir = 0; ir < o2i_nchunks_ - (o2i_selfcomm_ >= 0 ); ++ir) {
        MemChunk      *cchunk = o2i_chunks_ + ir;
        opt_double_ptr buf    = cchunk->data;
        size_t         count  = cchunk->size_padded * cchunk->nda;
        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, cchunk->dest_rank, subcomm_, i2o_recv_rqst_ + ir);
        MPI_Send_init(buf, (int)(count), MPI_DOUBLE, cchunk->dest_rank, cchunk->dest_rank, subcomm_, o2i_send_rqst_ + ir);
    }

    //..........................................................................
    // Create the self requests
    if(i2o_selfcomm_ >= 0){
        // Setting everything that use the i2o_chunks 
        {
            MemChunk      *cchunk = i2o_chunks_ + i2o_selfcomm_;
            opt_double_ptr buf    = cchunk->data;
            size_t         count  = cchunk->size_padded * cchunk->nda;
            FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
            FLUPS_CHECK(cchunk->dest_rank == sub_rank,"the dest rank should be equal to the subrank when performing self request");
            MPI_Send_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, i2o_send_rqst_ + i2o_selfcomm_);
            MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, o2i_recv_rqst_ + i2o_selfcomm_);
        }

        // Setting everything that use the o2i_chunks 
        {
            MemChunk *cchunk      = o2i_chunks_ + o2i_selfcomm_;
            opt_double_ptr buf    = cchunk->data;
            size_t         count  = cchunk->size_padded * cchunk->nda;
            FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
            FLUPS_CHECK(cchunk->dest_rank == sub_rank,"the dest rank should be equal to the subrank when performing self request");
            MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, i2o_recv_rqst_ + o2i_selfcomm_);
            MPI_Send_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, o2i_send_rqst_ + o2i_selfcomm_);
        }

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
    m_profStarti(prof_, "Switchtopo%d_%s", idswitchtopo_, (FLUPS_FORWARD == sign) ? "forward" : "backward");
    if (sign == FLUPS_FORWARD) {
        SendRecv(i2o_nchunks_, i2o_send_rqst_, i2o_chunks_, i2o_selfcomm_,
                 o2i_nchunks_, i2o_recv_rqst_, o2i_chunks_, o2i_selfcomm_,
                 topo_in_, topo_out_, v, prof_);
    } else {
        SendRecv(o2i_nchunks_, o2i_send_rqst_, o2i_chunks_, o2i_selfcomm_,
                 i2o_nchunks_, o2i_recv_rqst_, i2o_chunks_, i2o_selfcomm_,
                 topo_out_, topo_in_, v, prof_);
    }
    m_profStopi(prof_, "Switchtopo%d_%s", idswitchtopo_, (FLUPS_FORWARD == sign) ? "forward" : "backward");

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
void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks, const int self_send_idx,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks, const int self_recv_idx,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof) {
    BEGIN_FUNC;
    FLUPS_CHECK(self_send_idx*self_recv_idx >= 0, "the self idx must be the same for the send and receive" );
    //--------------------------------------------------------------------------
    // Get the memory arrangement
    const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    // Get the number of request we have to perform with other ranks than myself
    const int n_other_send = n_send_rqst - (self_send_idx >= 0 );
    const int n_other_recv = n_recv_rqst - (self_recv_idx >= 0 );

    //..........................................................................
    // Define the tasks to perform while proceeding to the MPI send and recv
    auto send_my_rqst = [=](int count, MemChunk *chunk, MPI_Request *request) {
        FLUPS_INFO("sending request to rank %d of size %d %d %d", chunk->dest_rank, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
        // copy here the chunk from the input topo to the chunk
        m_profStart(prof, "copy");
        for (int ic = 0; ic < count; ++ic) {
            CopyData2Chunk(nmem_in, mem, chunk + ic);
        }
        m_profStop(prof, "copy");
        // start the request
        MPI_Startall(count, request);
    };
    auto recv_my_rqst = [prof](MPI_Request *request, MemChunk *chunk) {
        FLUPS_INFO("recving request from rank %d of size %d %d %d", chunk->dest_rank, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
        // shuffle the data
        DoShuffleChunk(chunk);
    };

    // Define the send of the batch
    auto send_my_batch = [&send_my_rqst, prof](int n_ttl_to_send, int *n_already_send, int n_to_send, MemChunk *chunks, MPI_Request *request) {
        int count_send = m_min(n_ttl_to_send - *n_already_send, n_to_send);
        FLUPS_INFO("sending %d/%d request -- already send = %d", count_send, n_ttl_to_send, *n_already_send);
        if (count_send > 0) {
            // m_profStarti(prof, "copy data 2 chunk");
            send_my_rqst(count_send, chunks + *n_already_send, request + *n_already_send);
            // m_profStopi(prof, "copy data 2 chunk");
        }
        *n_already_send += count_send;
    };

    //..........................................................................
    // Define the counter needed to perform the send and receive
    const int send_batch   = FLUPS_MPI_BATCH_SEND;
    int       send_cntr    = 0;
    int       recv_cntr    = 0;
    int      *completed_id = reinterpret_cast<int *>(m_calloc(n_other_recv * sizeof(int)));

    //..........................................................................
    m_profStart(prof, "send/recv");

    m_profStart(prof, "pre-send");
    m_profInitLeave(prof, "copy");
    // m_profInitLeave(prof, "shuffle");
    {
        // According to the standard 3.1 pg 77, MPI should start all the requests in the array
        // so we start all the other request and the self request using the same start.
        FLUPS_INFO("starting %d recv request", n_recv_rqst);
        MPI_Startall(n_recv_rqst, recv_rqst);

        // Start a first batch of send request
        send_my_batch(n_other_send, &send_cntr, send_batch, send_chunks, send_rqst);
    }
    m_profStop(prof, "pre-send");

    //..........................................................................
    // make sure everybody enter all the profiler counters
    m_profStart(prof, "self");
    m_profInitLeave(prof, "copy");
    m_profInitLeave(prof, "shuffle");
    // process the self request
    if (0 <= self_send_idx) {
        FLUPS_INFO("sending and processing the self request");
        // send my request
        send_my_rqst(1, send_chunks + self_send_idx, send_rqst + self_send_idx);
        // wait for the completion of the recv
        MPI_Wait(recv_rqst + self_recv_idx, MPI_STATUSES_IGNORE);

        // process the recv
        m_profStart(prof, "shuffle");
        recv_my_rqst(recv_rqst + self_recv_idx, recv_chunks + self_recv_idx);
        m_profStop(prof, "shuffle");

        // finish the send
        MPI_Wait(send_rqst + self_send_idx, MPI_STATUSES_IGNORE);
    }
    m_profStop(prof, "self");
    
    // while we still have to send or recv something, we continue
    m_profStart(prof, "while loop");
    m_profInitLeave(prof, "copy");
    m_profInitLeave(prof, "shuffle");
    while ((send_cntr < n_other_send) || (recv_cntr < n_other_recv)) {
        // if we have some requests to recv, test it
        int n_completed = 0;
        if (recv_cntr < n_recv_rqst) {
            MPI_Testsome(n_recv_rqst, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);
        }

        // Maintain the difference between the pending send and receive request to send_batch
        // if we completed n_completed recvs it means that the send are completed on the other rank.
        // to maintain the total balance accross the network we start another chunk of n_completed
        send_my_batch(n_other_send, &send_cntr, n_completed, send_chunks, send_rqst);

        // for each of the completed request, treat it
        for (int id = 0; id < n_completed; ++id) {
            FLUPS_INFO("recving request %d/%d", recv_cntr + id, n_recv_rqst);
            const int rqst_id = completed_id[id];
            m_profStart(prof, "shuffle");
            recv_my_rqst(recv_rqst + rqst_id, recv_chunks + rqst_id);
            m_profStop(prof, "shuffle");
        }

        // increment the counter of recv request
        recv_cntr += n_completed;
    }
    m_profStop(prof, "while loop");
    FLUPS_INFO("Send recv completed");

    //..........................................................................
    // we need to officially close the send requests
    // m_profStart(prof, "Wait all");
    MPI_Waitall(n_other_send, send_rqst, MPI_STATUSES_IGNORE);
    // m_profStop(prof, "Wait all");
    FLUPS_INFO("Wait all completed");

    //..........................................................................
    // reset the memory to 0.0 as we do inplace computations
    const size_t reset_size = topo_out->memsize();
    std::memset(mem, 0, reset_size * sizeof(double));
    FLUPS_INFO("reset mem done ");

    // once all the send has been done we can overwrite the received info
    FLUPS_INFO("Copying data");
    m_profStart(prof, "final copy");
    m_profStart(prof, "copy");
    for (int ic = 0; ic < n_recv_rqst; ++ic) {
        CopyChunk2Data(recv_chunks + ic, nmem_out, mem);
    }
    m_profStop(prof, "copy");
    m_profStop(prof, "final copy");
    m_profStop(prof, "send/recv");

    FLUPS_INFO("Copying data completed");

    // m_free(completed_status);
    m_free(completed_id);
    //--------------------------------------------------------------------------
    END_FUNC;
}
