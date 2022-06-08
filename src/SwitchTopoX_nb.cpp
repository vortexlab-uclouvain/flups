#include "SwitchTopoX_nb.hpp"

void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks, const int *send_order_list,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks, int *completed_id,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler *prof) ;

SwitchTopoX_nb::SwitchTopoX_nb(const Topology *topo_in, const Topology *topo_out, const int shift[3], H3LPR::Profiler *prof)
    : SwitchTopoX(topo_in, topo_out, shift, prof) {
    BEGIN_FUNC;
    FLUPS_CHECK(sendData != nullptr, "The send data must be != to nullptr");
    FLUPS_CHECK(recvData != nullptr, "The recv data must be != to nullptr");
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

    //..........................................................................
    // get information on shared rank
    int sub_rank;
    MPI_Comm_rank(subcomm_, &sub_rank);
    MPI_Comm_split_type(subcomm_, MPI_COMM_TYPE_SHARED, sub_rank, MPI_INFO_NULL, &shared_comm_);

#if (0 == FLUPS_OLD_MPI)
    // apply some fancy parameters to allow faster MPI calls if we have a MPI-4.0 compliant version
    // the info is NOT transfered from one comm to another
    MPI_Info info;
    MPI_Info_create(&info);
    // MPI_Info_set(info, "mpi_assert_exact_length", "true");
    MPI_Info_set(info, "mpi_assert_allow_overtaking", "true");
    MPI_Info_set(info, "mpi_assert_no_any_tag", "true");
    MPI_Info_set(info, "mpi_assert_no_any_source", "true");
    MPI_Comm_set_info(shared_comm_, info);
    MPI_Info_free(&info);
#endif

    int shared_rank;
    MPI_Comm_rank(shared_comm_, &shared_rank);

    //..........................................................................
    // once we have the MemChunks we can allocate the requests
    // i2o transfert goes from i2o_chunks to o2i_chunks
    // We do not use MPI to process the self request so we have to remove it from the number of request
    i2o_send_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Request)));
    i2o_recv_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Request)));
    // o2i transfert goes from o2i_chunks to i2o_chunks
    o2i_send_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Request)));
    o2i_recv_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Request)));

    i2o_send_order_ = reinterpret_cast<int *>(m_calloc(i2o_nchunks_ * sizeof(int)));
    o2i_send_order_ = reinterpret_cast<int *>(m_calloc(o2i_nchunks_ * sizeof(int)));

    // allocate the completed_id array
    const int n_rqst = m_max(i2o_nchunks_, o2i_nchunks_);
    completed_id_    = reinterpret_cast<int *>(m_calloc(n_rqst * sizeof(int)));

    //..........................................................................
    // try to narrow down the comm group, if possible
    MPI_Group sub_group;
    MPI_Group shared_group;
    MPI_Comm_group(subcomm_, &sub_group);
    MPI_Comm_group(shared_comm_, &shared_group);

    // we store the priority requests in front of the order list
    int i2o_prior_idx = 0;
    int o2i_prior_idx = 0;
    // we store the non-priority chunks at the end of the order list
    int i2o_noprior_idx = i2o_nchunks_ - 1;
    int o2i_noprior_idx = o2i_nchunks_ - 1;

    auto opinit = [=](const int nchunks, const int selfcomm,
                      MemChunk *chunks, MPI_Request *send_rqst, MPI_Request *recv_rqst,
                      int *send_order, int *prior_idx, int *noprior_idx) {
        //......................................................................
        const int n_requests = nchunks - (selfcomm >= 0);
        for (int ir = 0; ir < n_requests; ++ir) {
            // we offset the starting index to avoid congestion
            const int ichunk = (ir + sub_rank) % n_requests;
            // get the chunk informations
            MemChunk      *cchunk = chunks + ichunk;
            opt_double_ptr buf    = cchunk->data;
            size_t         count  = cchunk->size_padded * cchunk->nda;

            // test if the destination rank is inside the shared group
            int dest_shared_rank;
            MPI_Group_translate_ranks(sub_group, 1, &(cchunk->dest_rank), shared_group, &dest_shared_rank);
            bool is_in_shared = (MPI_UNDEFINED != dest_shared_rank);
            // if the dest_rank is inside the shared group, use the shared comm
            int      dest_rank = is_in_shared ? dest_shared_rank : cchunk->dest_rank;
            int      send_tag  = is_in_shared ? shared_rank : sub_rank;
            int      recv_tag  = is_in_shared ? dest_shared_rank : cchunk->dest_rank;
            MPI_Comm dest_comm = is_in_shared ? shared_comm_ : subcomm_;

            FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
            MPI_Send_init(buf, (int)(count), MPI_DOUBLE, dest_rank, send_tag, dest_comm, send_rqst + ichunk);
            MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, dest_rank, recv_tag, dest_comm, recv_rqst + ichunk);

            // store the id in the send order list
            if (!is_in_shared) {
                send_order[prior_idx[0]] = ichunk;
                (prior_idx[0])++;
            } else {
                send_order[noprior_idx[0]] = ichunk;
                (noprior_idx[0])--;
            }
        }
        //......................................................................
        // Create the self requests first
        if (selfcomm >= 0) {
            const int      ichunk = selfcomm;
            MemChunk      *cchunk = chunks + ichunk;
            opt_double_ptr buf    = cchunk->data;
            size_t         count  = cchunk->size_padded * cchunk->nda;
            FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
            FLUPS_CHECK(cchunk->dest_rank == sub_rank, "the dest rank should be equal to the subrank when performing self request");
            MPI_Send_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, send_rqst + ichunk);
            MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, recv_rqst + ichunk);

            // store the id in the send order list
            // send_order[prior_idx[0]] = ichunk;
            // (prior_idx[0])++;
            send_order[noprior_idx[0]] = ichunk;
            (noprior_idx[0])--;
        }
        FLUPS_CHECK(noprior_idx[0] == (prior_idx[0] - 1), "the prior index = %d should be = %d + 1 = %d", prior_idx[0], noprior_idx[0], noprior_idx[0] + 1);
    };

    opinit(i2o_nchunks_, i2o_selfcomm_, i2o_chunks_, i2o_send_rqst_, o2i_recv_rqst_, i2o_send_order_, &i2o_prior_idx, &i2o_noprior_idx);
    opinit(o2i_nchunks_, o2i_selfcomm_, o2i_chunks_, o2i_send_rqst_, i2o_recv_rqst_, o2i_send_order_, &o2i_prior_idx, &o2i_noprior_idx);

    // free the groups
    MPI_Group_free(&sub_group);
    MPI_Group_free(&shared_group);

    // //..........................................................................
    // // this is the loop over the input topo and the associated chunks
    // const int i2o_n_requests = i2o_nchunks_ - (i2o_selfcomm_ >= 0);
    // for (int ir = 0; ir < i2o_n_requests; ++ir) {
    //     // we offset the starting index to avoid congestion
    //     const int      ichunk = (ir + sub_rank) % i2o_n_requests;
    //     MemChunk      *cchunk = i2o_chunks_ + ichunk;
    //     opt_double_ptr buf    = cchunk->data;
    //     size_t         count  = cchunk->size_padded * cchunk->nda;

    //     // test if the destination rank is inside the shared group
    //     int dest_shared_rank;
    //     MPI_Group_translate_ranks(sub_group, 1, &(cchunk->dest_rank), shared_group, &dest_shared_rank);
    //     bool is_in_shared = (MPI_UNDEFINED != dest_shared_rank);
    //     // if the dest_rank is inside the shared group, use the shared comm
    //     int      dest_rank = is_in_shared ? dest_shared_rank : cchunk->dest_rank;
    //     int      tag       = is_in_shared ? shared_rank : sub_rank;
    //     MPI_Comm dest_comm = is_in_shared ? shared_comm_ : subcomm_;

    //     FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
    //     MPI_Send_init(buf, (int)(count), MPI_DOUBLE, dest_rank, tag, dest_comm, i2o_send_rqst_ + ichunk);
    //     MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, dest_rank, tag, dest_comm, o2i_recv_rqst_ + ichunk);

    //     // store the id in the send order list
    //     if (!is_in_shared) {
    //         i2o_send_order_[i2o_prior_idx] = ichunk;
    //         i2o_prior_idx++;
    //     } else {
    //         i2o_send_order_[i2o_noprior_idx] = ichunk;
    //         i2o_noprior_idx--;
    //     }
    // }

    // // here we go for the output topo and the associated chunks
    // const int o2i_n_requests = o2i_nchunks_ - (o2i_selfcomm_ >= 0);
    // for (int ir = 0; ir < o2i_n_requests; ++ir) {
    //     // we offset the starting index to avoid congestion
    //     const int      ichunk =  (ir + sub_rank) % o2i_n_requests;
    //     MemChunk      *cchunk = o2i_chunks_ + ichunk;
    //     opt_double_ptr buf    = cchunk->data;
    //     size_t         count  = cchunk->size_padded * cchunk->nda;

    //     // test if the destination rank is inside the shared group
    //     int dest_shared_rank;
    //     MPI_Group_translate_ranks(sub_group, 1, &(cchunk->dest_rank), shared_group, &dest_shared_rank);
    //     bool is_in_shared = (MPI_UNDEFINED != dest_shared_rank);
    //     // if the dest_rank is inside the shared group, use the shared comm
    //     int      dest_rank = is_in_shared ? dest_shared_rank : cchunk->dest_rank;
    //     int      tag       = is_in_shared ? dest_shared_rank : cchunk->dest_rank;
    //     MPI_Comm dest_comm = is_in_shared ? shared_comm_ : subcomm_;

    //     FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
    //     MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, dest_rank, tag, dest_comm, i2o_recv_rqst_ + ichunk);
    //     MPI_Send_init(buf, (int)(count), MPI_DOUBLE, dest_rank, tag, dest_comm, o2i_send_rqst_ + ichunk);

    //     // store the id in the send order list
    //     if (!is_in_shared) {
    //         o2i_send_order_[o2i_prior_idx] = ichunk;
    //         o2i_prior_idx++;
    //     } else {
    //         o2i_send_order_[o2i_noprior_idx] = ichunk;
    //         o2i_noprior_idx--;
    //     }
    // }
    // //..........................................................................
    // // Create the self requests first
    // if (i2o_selfcomm_ >= 0) {
    //     // Setting everything that use the i2o_chunks
    //     {
    //         const int      ichunk = i2o_selfcomm_;
    //         MemChunk      *cchunk = i2o_chunks_ + ichunk;
    //         opt_double_ptr buf    = cchunk->data;
    //         size_t         count  = cchunk->size_padded * cchunk->nda;
    //         FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
    //         FLUPS_CHECK(cchunk->dest_rank == sub_rank, "the dest rank should be equal to the subrank when performing self request");
    //         MPI_Send_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, i2o_send_rqst_ + ichunk);
    //         MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, o2i_recv_rqst_ + ichunk);

    //         // store the id in the send order list
    //         // i2o_send_order_[i2o_prior_idx] = ichunk;
    //         // i2o_prior_idx++;
    //         i2o_send_order_[i2o_noprior_idx] = ichunk;
    //         i2o_noprior_idx--;
    //     }

    //     // Setting everything that use the o2i_chunks
    //     {
    //         const int      ichunk = o2i_selfcomm_;
    //         MemChunk      *cchunk = o2i_chunks_ + o2i_selfcomm_;
    //         opt_double_ptr buf    = cchunk->data;
    //         size_t         count  = cchunk->size_padded * cchunk->nda;
    //         FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
    //         FLUPS_CHECK(cchunk->dest_rank == sub_rank, "the dest rank should be equal to the subrank when performing self request");
    //         MPI_Recv_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, i2o_recv_rqst_ + ichunk);
    //         MPI_Send_init(buf, (int)(count), MPI_DOUBLE, 0, 0, MPI_COMM_SELF, o2i_send_rqst_ + ichunk);

    //         // store the id in the send order list
    //         // o2i_send_order_[o2i_prior_idx] = ichunk;
    //         // o2i_prior_idx++;
    //         o2i_send_order_[o2i_noprior_idx] = ichunk;
    //         o2i_noprior_idx--;
    //     }
    // }
    // FLUPS_CHECK(i2o_noprior_idx == (i2o_prior_idx - 1), "the prior index = %d should be = %d + 1 = %d", i2o_prior_idx, i2o_noprior_idx, i2o_noprior_idx + 1);
    // FLUPS_CHECK(o2i_noprior_idx == (o2i_prior_idx - 1), "the prior index = %d should be = %d + 1 = %d", o2i_prior_idx, o2i_noprior_idx, o2i_noprior_idx + 1);

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

    m_free(i2o_send_order_);
    m_free(o2i_send_order_);
    m_free(completed_id_);

    MPI_Comm_free(&shared_comm_);
    // int comp;
    // MPI_Comm_compare(shared_comm_,inComm_,&comp);
    // if(comp!=MPI_IDENT){
    //     MPI_Comm_free(&subcomm_);
    // }
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
        SendRecv(i2o_nchunks_, i2o_send_rqst_, i2o_chunks_, i2o_send_order_,
                 o2i_nchunks_, i2o_recv_rqst_, o2i_chunks_, completed_id_,
                 topo_in_, topo_out_, v, prof_);
    } else {
        SendRecv(o2i_nchunks_, o2i_send_rqst_, o2i_chunks_, o2i_send_order_,
                 i2o_nchunks_, o2i_recv_rqst_, i2o_chunks_, completed_id_,
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
 * @param send_order_list
 * @param n_recv_rqst
 * @param recv_rqst
 * @param recv_chunks
 * @param topo_in
 * @param topo_out
 * @param mem
 */
void SendRecv(const int n_send_rqst, MPI_Request *send_rqst, MemChunk *send_chunks, const int *send_order_list,
              const int n_recv_rqst, MPI_Request *recv_rqst, MemChunk *recv_chunks, int *completed_id,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler *prof) {
    BEGIN_FUNC;
    // FLUPS_CHECK(self_send_idx*self_recv_idx >= 0, "the self idx must be the same for the send and receive" );
    //--------------------------------------------------------------------------
    // Get the memory arrangement
    const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    // Get the number of request we have to perform with other ranks than myself
    // const int n_other_send = n_send_rqst - (self_send_idx >= 0);
    // const int n_other_recv = n_recv_rqst - (self_recv_idx >= 0);

    //..........................................................................
    // Define the tasks to perform while proceeding to the MPI send and recv
    // auto send_my_rqst = [=](const int count, MemChunk *chunk, MPI_Request *request) {
    //     FLUPS_INFO("sending request to rank %d of size %d %d %d", chunk->dest_rank, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
    //     // copy here the chunk from the input topo to the chunk
    //     for (int ic = 0; ic < count; ++ic) {
    //         m_profStart(prof, "copy");
    //         CopyData2Chunk(nmem_in, mem, chunk + ic);
    //         m_profStop(prof, "copy");
    //     }
    //     // start the request
    //     m_profStart(prof, "start");
    //     MPI_Startall(count, request);
    //     m_profStop(prof, "start");
    // };
    auto recv_my_rqst = [prof](MPI_Request *request, MemChunk *chunk) {
        FLUPS_INFO("recving request from rank %d of size %d %d %d", chunk->dest_rank, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
        // shuffle the data
        DoShuffleChunk(chunk);
    };

    // Define the send of a batch of requests
    auto send_my_batch = [=](const int n_ttl_to_send, int *n_already_send, const int n_batch,
                             const int *my_send_order, MemChunk *chunks, MPI_Request *request) {
        // determine how many requests are left to send
        int count_send = m_min(n_ttl_to_send - n_already_send[0], n_batch);
        FLUPS_CHECK(count_send >= 0, "count send = %d cannot be negative", count_send);
        FLUPS_INFO("sending %d/%d request -- already send = %d", count_send, n_ttl_to_send, *n_already_send);

        // for all of the requests, loop one to one and send them
        // get the starting index of the request
        for (int ir = 0; ir < count_send; ++ir) {
            // get the request id
            const int   *id_to_send = my_send_order + n_already_send[0] + ir;
            MPI_Request *c_rqst     = request + id_to_send[0];
            MemChunk    *c_chunk    = chunks + id_to_send[0];
            
            // copy the memory
            m_profStart(prof, "copy");
            CopyData2Chunk(nmem_in, mem, c_chunk);
            m_profStop(prof, "copy");

            // start the send
            m_profStart(prof, "start");
            MPI_Start(c_rqst);
            m_profStop(prof, "start");
        }
        // increment the send counter
        n_already_send[0] += count_send;
    };

    //..........................................................................
    // Define the counter needed to perform the send and receive
    const int send_batch   = FLUPS_MPI_BATCH_SEND;
    int       send_cntr    = 0;
    int       recv_cntr    = 0;

    //..........................................................................
    m_profStart(prof, "send/recv");

    m_profStart(prof, "pre-send");
    m_profInitLeave(prof, "copy");
    m_profInitLeave(prof, "start");
    {
        // According to the standard 3.1 pg 77, MPI should start all the requests in the array, independently from the communicator used
        // so we start all the other request and the self request using the same start.
        FLUPS_INFO("starting %d recv request", n_recv_rqst);
        m_profStart(prof, "start");
        MPI_Startall(n_recv_rqst, recv_rqst);
        m_profStop(prof, "start");

        // Start a first batch of send request
        // here we use n_other_send as the self will be processed!
        // send_my_batch(n_other_send, &send_cntr, send_batch, send_chunks, send_rqst);
        send_my_batch(n_send_rqst, &send_cntr, send_batch, send_order_list, send_chunks, send_rqst);
    }
    m_profStop(prof, "pre-send");

    //..........................................................................
    // make sure everybody enter all the profiler counters
    // m_profStart(prof, "self");
    // m_profInitLeave(prof, "copy");
    // m_profInitLeave(prof, "start");
    // // process the self send request
    // if (0 <= self_send_idx) {
    //     FLUPS_INFO("sending the self request");
    //     // we cannot count the self request into the send-counter :-) (yes, I know, it's cryptic)
    //     send_my_rqst(1, send_chunks + self_send_idx, send_rqst + self_send_idx);
    // }
    // m_profStop(prof, "self");

    // while we still have to send or recv something, we continue
    m_profStart(prof, "while loop");
    m_profInitLeave(prof, "copy");
    m_profInitLeave(prof, "start");
    m_profInitLeave(prof, "shuffle");
    // while we have to send msgs to others or recv msg, we keep going
    while ((send_cntr < n_send_rqst) || (recv_cntr < n_recv_rqst)) {
        // if we have some requests to recv, test it
        int n_completed = 0;
        if (recv_cntr < n_recv_rqst) {
            // MPI_Waitsome(n_recv_rqst, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);
            MPI_Testsome(n_recv_rqst, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);
        }

        // here we use the n_completed information as an estimation of the speed at which I receive requests
        // if I haven't received any, just re-send a batch.
        // if I have received many, then the throughput is great and I should resend many
        const int n_to_resend = m_max(n_completed,send_batch);
        send_my_batch(n_send_rqst, &send_cntr, n_to_resend, send_order_list, send_chunks, send_rqst);

        // for each of the completed request, treat it
        for (int id = 0; id < n_completed; ++id) {
            FLUPS_INFO("recving request %d/%d with id = %d", recv_cntr + id, n_recv_rqst, completed_id[id]);
            const int rqst_id = completed_id[id];
            m_profStart(prof, "shuffle");
            recv_my_rqst(recv_rqst + rqst_id, recv_chunks + rqst_id);
            m_profStop(prof, "shuffle");
        }

        // increment the counter of recv request
        recv_cntr += n_completed;
    }
    m_profStop(prof, "while loop");

    //..........................................................................
    // we need to officially close ALL the send requests
    MPI_Waitall(n_send_rqst, send_rqst, MPI_STATUSES_IGNORE);

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
    //--------------------------------------------------------------------------
    END_FUNC;
}
