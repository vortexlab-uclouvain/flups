#include "SwitchTopoX_isr.hpp"

void SendRecv(const int n_send_chunk, MPI_Request *send_rqst, MemChunk *send_chunks, MPI_Datatype *send_dtype, size_t *send_offset,
              const int n_recv_chunk, MPI_Request *recv_rqst, MemChunk *recv_chunks, MPI_Datatype *recv_dtype,
              const int *send_order_list, int *completed_id, MPI_Comm comm,
              const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler *prof);

SwitchTopoX_isr::SwitchTopoX_isr(const Topology *topo_in, const Topology *topo_out, const int shift[3], H3LPR::Profiler *prof)
    : SwitchTopoX(topo_in, topo_out, shift, prof) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // nothing special to do here
    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_isr::setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) {
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
    // once we have the MemChunks we can allocate the communication information
    const int n_rqst = m_max(i2o_nchunks_, o2i_nchunks_);
    completed_id_    = reinterpret_cast<int *>(m_calloc(n_rqst * sizeof(int)));
    send_rqst_       = reinterpret_cast<MPI_Request *>(m_calloc(n_rqst * sizeof(MPI_Request)));
    recv_rqst_       = reinterpret_cast<MPI_Request *>(m_calloc(n_rqst * sizeof(MPI_Request)));

    i2o_send_order_ = reinterpret_cast<int *>(m_calloc(i2o_nchunks_ * sizeof(int)));
    o2i_send_order_ = reinterpret_cast<int *>(m_calloc(o2i_nchunks_ * sizeof(int)));

    i2o_offset_     = reinterpret_cast<size_t *>(m_calloc(i2o_nchunks_ * sizeof(size_t)));
    o2i_offset_     = reinterpret_cast<size_t *>(m_calloc(o2i_nchunks_ * sizeof(size_t)));

    i2o_dtype_      = reinterpret_cast<MPI_Datatype *>(m_calloc(i2o_nchunks_ * sizeof(MPI_Datatype)));
    o2i_dtype_      = reinterpret_cast<MPI_Datatype *>(m_calloc(o2i_nchunks_ * sizeof(MPI_Datatype)));

    //..........................................................................
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

    auto opinit = [=](const int nmem[3],                     /* topology info */
                      size_t *offsets, MPI_Datatype *dtypes, /* types and offsets */
                      const int nchunks, MemChunk *chunks,   /* chunk info */
                      int *send_order, int *prior_idx, int *noprior_idx) {
        for (int ir = 0; ir < nchunks; ++ir) {
            // we offset the starting index to avoid congestion
            const int ichunk = (ir + sub_rank) % nchunks;
            // get the chunk informations
            MemChunk *cchunk = chunks + ichunk;
            // opt_double_ptr buf    = cchunk->data;
            // size_t         count  = cchunk->size_padded * cchunk->nda;

            // test if the destination rank is inside the shared group
            int dest_shared_rank;
            MPI_Group_translate_ranks(sub_group, 1, &(cchunk->dest_rank), shared_group, &dest_shared_rank);
            const bool is_in_shared = (MPI_UNDEFINED != dest_shared_rank);

            // store the datatype
            ChunkToMPIDataType(nmem, cchunk, offsets + ichunk, dtypes + ichunk);

            // store the id in the send order list
            if (!is_in_shared) {
                send_order[prior_idx[0]] = ichunk;
                (prior_idx[0])++;
            } else {
                send_order[noprior_idx[0]] = ichunk;
                (noprior_idx[0])--;
            }
        }
    };

    const int in_nmem[3]  = {topo_in_->nmem(0), topo_in_->nmem(1), topo_in_->nmem(2)};
    const int out_nmem[3] = {topo_out_->nmem(0), topo_out_->nmem(1), topo_out_->nmem(2)};
    opinit(in_nmem, i2o_offset_, i2o_dtype_, i2o_nchunks_, i2o_chunks_, i2o_send_order_, &i2o_prior_idx, &i2o_noprior_idx);
    opinit(out_nmem, o2i_offset_, o2i_dtype_, o2i_nchunks_, o2i_chunks_, o2i_send_order_, &o2i_prior_idx, &o2i_noprior_idx);
    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX_isr::~SwitchTopoX_isr(){
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    for (int ir = 0; ir < i2o_nchunks_; ++ir) {
        MPI_Type_free(i2o_dtype_ + ir);
    }
    // here we go for the output topo and the associated chunks
    for (int ir = 0; ir < o2i_nchunks_; ++ir) {
        MPI_Type_free(o2i_dtype_ + ir);
    }

    // free the request arrays
    m_free(send_rqst_);
    m_free(recv_rqst_);
    m_free(i2o_dtype_);
    m_free(i2o_offset_);
    m_free(o2i_dtype_);
    m_free(o2i_offset_);

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
void SwitchTopoX_isr::execute(opt_double_ptr v, const int sign) const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    m_profStarti(prof_, "Switchtopo%d_%s", idswitchtopo_, (FLUPS_FORWARD == sign) ? "forward" : "backward");
    if (sign == FLUPS_FORWARD) {
        SendRecv(i2o_nchunks_, send_rqst_, i2o_chunks_, i2o_dtype_, i2o_offset_,
                 o2i_nchunks_, recv_rqst_, o2i_chunks_, o2i_dtype_,
                 i2o_send_order_, completed_id_, subcomm_,
                 topo_out_, v, prof_);

    } else {
        SendRecv(o2i_nchunks_, recv_rqst_, o2i_chunks_, o2i_dtype_, o2i_offset_,
                 i2o_nchunks_, send_rqst_, i2o_chunks_, i2o_dtype_,
                 o2i_send_order_, completed_id_, subcomm_,
                 topo_in_, v, prof_);
    }
    m_profStopi(prof_, "Switchtopo%d_%s", idswitchtopo_, (FLUPS_FORWARD == sign) ? "forward" : "backward");
    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_isr::disp() const {
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

void SendRecv(const int n_send_chunk, MPI_Request *send_rqst, MemChunk *send_chunks, MPI_Datatype *send_dtype, size_t *send_offset,
              const int n_recv_chunk, MPI_Request *recv_rqst, MemChunk *recv_chunks, MPI_Datatype *recv_dtype,
              const int *send_order_list, int *completed_id, MPI_Comm comm,
              const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler *prof) {
    BEGIN_FUNC;
    // FLUPS_CHECK(self_send_idx*self_recv_idx >= 0, "the self idx must be the same for the send and receive" );
    //--------------------------------------------------------------------------
    // Get the memory arrangement
    // const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    // const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    // Get the number of request we have to perform with other ranks than myself
    // const int n_other_send = n_send_rqst - (self_send_idx >= 0);
    // const int n_other_recv = n_recv_rqst - (self_recv_idx >= 0);

    int rank;
    MPI_Comm_rank(comm,&rank);

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
    auto send_my_batch = [=](const int n_ttl_to_send, int *n_already_send, const int n_batch) {
        // determine how many requests are left to send
        int count_send = m_min(n_ttl_to_send - n_already_send[0], n_batch);
        FLUPS_CHECK(count_send >= 0, "count send = %d cannot be negative", count_send);
        FLUPS_INFO("sending %d/%d request -- already send = %d", count_send, n_ttl_to_send, *n_already_send);

        // for all of the requests, loop one to one and send them
        // get the starting index of the request
        for (int ir = 0; ir < count_send; ++ir) {
            const int     ridx      = n_already_send[0] + ir;
            const int     chunk_idx = send_order_list[ridx];
            size_t        offset    = send_offset[chunk_idx];
            MPI_Datatype *dtype     = send_dtype + chunk_idx;
            MemChunk     *c_chunk   = send_chunks + chunk_idx;

            // send is done directly from the memory to MPI
            m_profStart(prof, "start");
            MPI_Isend(mem + offset, 1, dtype[0], c_chunk->dest_rank, rank, comm, send_rqst + ridx);
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
        FLUPS_INFO("starting %d recv request", n_recv_chunk);
        m_profStart(prof, "start");
        for (int ir = 0; ir < n_recv_chunk; ++ir) {
            MPI_Datatype *dtype   = recv_dtype + ir;
            MemChunk     *c_chunk = recv_chunks + ir;
            MPI_Irecv(c_chunk->data, 1, dtype[0], c_chunk->dest_rank, c_chunk->dest_rank, comm, recv_rqst + ir);
        }
        m_profStop(prof, "start");

        // Start a first batch of send request
        // here we use n_other_send as the self will be processed!
        // send_my_batch(n_other_send, &send_cntr, send_batch, send_chunks, send_rqst);
        send_my_batch(n_send_chunk, &send_cntr, send_batch);
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
    while ((send_cntr < n_send_chunk) || (recv_cntr < n_recv_chunk)) {
        // if we have some requests to recv, test it
        int n_completed = 0;
        if (recv_cntr < n_recv_chunk) {
            // MPI_Waitsome(n_recv_rqst, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);
            MPI_Testsome(n_recv_chunk, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);
        }

        // here we use the n_completed information as an estimation of the speed at which I receive requests
        // if I haven't received any, just re-send a batch.
        // if I have received many, then the throughput is great and I should resend many
        const int n_to_resend = m_max(n_completed,send_batch);
        send_my_batch(n_send_chunk, &send_cntr, n_to_resend);

        // for each of the completed request, treat it
        for (int id = 0; id < n_completed; ++id) {
            FLUPS_INFO("recving request %d/%d with id = %d", recv_cntr + id, n_recv_chunk, completed_id[id]);
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
    MPI_Waitall(n_send_chunk, send_rqst, MPI_STATUSES_IGNORE);

    //..........................................................................
    // reset the memory to 0.0 as we do inplace computations
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};
    const size_t reset_size = topo_out->memsize();
    std::memset(mem, 0, reset_size * sizeof(double));
    FLUPS_INFO("reset mem done ");

    // once all the send has been done we can overwrite the received info
    FLUPS_INFO("Copying data");
    m_profStart(prof, "final copy");
    m_profStart(prof, "copy");
    for (int ic = 0; ic < n_recv_chunk; ++ic) {
        CopyChunk2Data(recv_chunks + ic, nmem_out, mem);
    }
    m_profStop(prof, "copy");
    m_profStop(prof, "final copy");
    m_profStop(prof, "send/recv");

    FLUPS_INFO("Copying data completed");
    //--------------------------------------------------------------------------
    END_FUNC;
}
