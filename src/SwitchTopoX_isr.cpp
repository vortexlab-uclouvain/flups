/**
 * @file SwitchTopoX_isr.cpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/
#include "SwitchTopoX_isr.hpp"

void SendRecv(const int n_send_chunk, MPI_Request *send_rqst, MemChunk *send_chunks,
              const int n_recv_chunk, MPI_Request *recv_rqst, MemChunk *recv_chunks,
              const int *send_order_list, int *completed_id, int* recv_order_list,
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
    FLUPS_CHECK(sendData == nullptr, "The send data must be = to nullptr");
    FLUPS_CHECK(recvData != nullptr, "The recv data must be != to nullptr");
    //--------------------------------------------------------------------------
    // first setup the basic stuffs
    this->SwitchTopoX::setup_buffers(sendData, recvData);

    //..........................................................................
    // once we have the MemChunks we can allocate the communication information
    const int n_rqst = m_max(i2o_nchunks_, o2i_nchunks_);
    completed_id_    = reinterpret_cast<int *>(m_calloc(n_rqst * sizeof(int)));
    recv_order_      = reinterpret_cast<int *>(m_calloc(n_rqst * sizeof(int)));
    send_rqst_       = reinterpret_cast<MPI_Request *>(m_calloc(n_rqst * sizeof(MPI_Request)));
    recv_rqst_       = reinterpret_cast<MPI_Request *>(m_calloc(n_rqst * sizeof(MPI_Request)));
    i2o_send_order_  = reinterpret_cast<int *>(m_calloc(i2o_nchunks_ * sizeof(int)));
    o2i_send_order_  = reinterpret_cast<int *>(m_calloc(o2i_nchunks_ * sizeof(int)));

    //..........................................................................
    // get information on shared rank and setup the priority list
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

    MPI_Group shared_group;
    MPI_Comm_group(shared_comm_, &shared_group);

    auto setup_priority = [=](const int nchunks, MemChunk *chunks, int *send_order, int *prior_idx, int *noprior_idx) {
        for (int ir = 0; ir < nchunks; ++ir) {
            // we offset the starting index to avoid congestion
#if (FLUPS_ROLLING_RANK)
            const int ichunk = (ir + sub_rank) % nchunks;
#else
            const int ichunk = ir;
#endif
            // get the chunk informations
            MemChunk *cchunk = chunks + ichunk;

            bool is_in_shared;
            ChunkToNewComm(shared_comm_, shared_group, cchunk, &is_in_shared);

            // store the id in the send order list and update the comm + dest_rank if needed
#if (FLUPS_PRIORITYLIST)
            if (is_in_shared) {
                // store the priority list
                send_order[noprior_idx[0]] = ichunk;
                (noprior_idx[0])--;
            } else {
                // store the priority
                send_order[prior_idx[0]] = ichunk;
                (prior_idx[0])++;
            }
#else
            send_order[prior_idx[0]] = ichunk;
            (prior_idx[0])++;
#endif 
        }
    };
    // we store the priority requests in front of the order list and the non-priority ones at the end of the list
    int i2o_prior_idx   = 0;
    int o2i_prior_idx   = 0;
    int i2o_noprior_idx = i2o_nchunks_ - 1;
    int o2i_noprior_idx = o2i_nchunks_ - 1;
    setup_priority(i2o_nchunks_, i2o_chunks_, i2o_send_order_, &i2o_prior_idx, &i2o_noprior_idx);
    setup_priority(o2i_nchunks_, o2i_chunks_, o2i_send_order_, &o2i_prior_idx, &o2i_noprior_idx);

    // free the groups
    MPI_Group_free(&shared_group);
    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX_isr::~SwitchTopoX_isr() {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // free the request arrays
    m_free(send_rqst_);
    m_free(recv_rqst_);
    m_free(i2o_send_order_);
    m_free(o2i_send_order_);
    m_free(completed_id_);
    m_free(recv_order_);

    MPI_Comm_free(&shared_comm_);
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
        SendRecv(i2o_nchunks_, send_rqst_, i2o_chunks_,
                 o2i_nchunks_, recv_rqst_, o2i_chunks_,
                 i2o_send_order_, completed_id_, recv_order_,
                 topo_out_, v, prof_);
    } else {
        SendRecv(o2i_nchunks_, send_rqst_, o2i_chunks_,
                 i2o_nchunks_, recv_rqst_, i2o_chunks_,
                 o2i_send_order_, completed_id_, recv_order_,
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

void SendRecv(const int n_send_chunk, MPI_Request *send_rqst, MemChunk *send_chunks,
              const int n_recv_chunk, MPI_Request *recv_rqst, MemChunk *recv_chunks,
              const int *send_order_list, int *completed_id, int* recv_order_list,
              const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler *prof) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // Define the send of a batch of requests
    auto send_my_batch = [=](const int n_ttl_to_send, int *n_already_send, const int n_batch) {
        // determine how many requests are left to send
        int count_send = m_min(n_ttl_to_send - n_already_send[0], n_batch);
        FLUPS_CHECK(count_send >= 0, "count send = %d cannot be negative", count_send);
        FLUPS_INFO("sending %d/%d request -- already send = %d", count_send, n_ttl_to_send, *n_already_send);

        // for all of the requests, loop one to one and send them
        // get the starting index of the request
        for (int ir = 0; ir < count_send; ++ir) {
            const int ridx      = n_already_send[0] + ir;
            const int chunk_idx = send_order_list[ridx];
            FLUPS_CHECK(chunk_idx < n_send_chunk, "the chunk id = %d must be < n_send_chunk = %d", chunk_idx, n_send_chunk);

            // send is done directly from the memory to MPI
            MemChunk *c_chunk = send_chunks + chunk_idx;
            int       rank_in_chunk;
            MPI_Comm_rank(c_chunk->comm, &rank_in_chunk);

            // start the Isend and store it using ridx to make sure we can test it later
            m_profStart(prof, "start");
            MPI_Isend(mem + c_chunk->offset, 1, c_chunk->dtype, c_chunk->dest_rank, rank_in_chunk, c_chunk->comm, send_rqst + ridx);
            m_profStop(prof, "start");
        }
        // increment the send counter
        n_already_send[0] += count_send;
        FLUPS_INFO(" I am all done here, moving on");
    };

    //..........................................................................
    // Define the counter needed to perform the send and receive
    const int send_batch    = FLUPS_MPI_BATCH_SEND;  // number of sends done at the same time
    int       send_cntr     = 0;                     // counter the number of send done
    int       recv_cntr     = 0;                     // count the number of recv completed
    int       copy_cntr     = 0;                     // count the number of processed received
    int       finished_send = 0;                     // count the number of completed send
    bool      is_mem_reset  = false;                 // track if the mem has been reset

    //..........................................................................
    m_profStart(prof, "send/recv");
    m_profStart(prof, "pre-send");
    m_profInitLeave(prof, "start");
    {
        // According to the standard 3.1 pg 77, MPI should start all the requests in the array, independently from the communicator used
        // so we start all the other request and the self request using the same start.
        FLUPS_INFO("starting %d recv request", n_recv_chunk);
        m_profStart(prof, "start");
        for (int ir = 0; ir < n_recv_chunk; ++ir) {
            MemChunk *c_chunk   = recv_chunks + ir;
            MPI_Irecv(c_chunk->data, 1, c_chunk->dest_dtype, c_chunk->dest_rank, c_chunk->dest_rank, c_chunk->comm, recv_rqst + ir);
        }
        m_profStop(prof, "start");

        
        send_my_batch(n_send_chunk, &send_cntr, send_batch);
    }
    m_profStop(prof, "pre-send");

    //..........................................................................
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};
    // while we still have to send or recv something, we continue
    m_profStart(prof, "while loop");
    m_profInitLeave(prof, "copy");
    m_profInitLeave(prof, "start");
    m_profInitLeave(prof, "shuffle");
    // while we have to send msgs to others or recv msg, we keep going
    while ((send_cntr < n_send_chunk) || (recv_cntr < n_recv_chunk) || (copy_cntr < n_recv_chunk)) {
        FLUPS_INFO("sent %d/%d - recvd %d/%d - copied %d/%d - reset done? %d", send_cntr, n_send_chunk, recv_cntr, n_recv_chunk, copy_cntr, n_recv_chunk, is_mem_reset);
        //......................................................................
        // [1] test if we have finished some send requests and start new ones
        //......................................................................
        if (finished_send < n_send_chunk) {
            // completed id can be reused here as it has been allocated on the max of send and recv
            int n_send_completed;
            MPI_Testsome(send_cntr, send_rqst, &n_send_completed, completed_id, MPI_STATUSES_IGNORE);

            // this is the total number of send that have completed
            finished_send += n_send_completed;
            const int still_ongoing_send = send_cntr - finished_send;
            const int n_to_resend        = m_min(FLUPS_MPI_MAX_NBSEND - still_ongoing_send, send_batch);
            FLUPS_CHECK(n_to_resend >= 0, " You need to send a positive number of request");
            send_my_batch(n_send_chunk, &send_cntr, n_to_resend);

            // if all the send have completed I can reset the memory to 0
            is_mem_reset = (finished_send == n_send_chunk);
            if (is_mem_reset) {
                const size_t reset_size = topo_out->memsize();
                std::memset(mem, 0, reset_size * sizeof(double));
                FLUPS_INFO("reset mem done ");
            }
        }
        //......................................................................
        // [2] test if we have finished some recv requests and shuffle them
        //......................................................................
        // if we have some requests to recv, test it
        if (recv_cntr < n_recv_chunk) {
            int n_completed = 0;
            MPI_Testsome(n_recv_chunk, recv_rqst, &n_completed, completed_id, MPI_STATUSES_IGNORE);

            // for each of the completed request save its id for processing later
            for (int id = 0; id < n_completed; ++id) {
                const int rqst_id = completed_id[id];
                FLUPS_INFO("shuffling request %d/%d with id = %d", recv_cntr + id, n_recv_chunk, rqst_id);
                // shuffle the data
                m_profStart(prof, "shuffle");
                DoShuffleChunk(recv_chunks + rqst_id);
                m_profStop(prof, "shuffle");
                // save the id for copy'ing it later
                recv_order_list[recv_cntr] = completed_id[id];
                recv_cntr++;
            }
        }

        //......................................................................
        // [3] test if we have finished some recv requests and shuffle them
        //......................................................................
        // for each of the ready  and not copied yet request, copy them
        if (is_mem_reset && (copy_cntr < recv_cntr)) {
            for(int icpy =  copy_cntr; icpy < recv_cntr; ++icpy){
                const int rqst_id = recv_order_list[icpy];
                FLUPS_INFO("treating recv request %d/%d with id = %d", icpy, n_recv_chunk, rqst_id);
                // copy the data
                m_profStart(prof, "copy");
                CopyChunk2Data(recv_chunks + rqst_id, nmem_out, mem);
                m_profStop(prof, "copy");
                ++copy_cntr;
            }
        }
    }
    m_profStop(prof, "while loop");
    m_profStop(prof, "send/recv");
    //--------------------------------------------------------------------------
    END_FUNC;
}
