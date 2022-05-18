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
    // Chunks are organised by rank so we loop over them and compute the counts
    // there is only one chunk per cpu so the displacement is obvious
    for (int ic = 0; ic < i2o_nchunks_; ++ic) {
        MemChunk *cchunk = i2o_chunks_ + ic;
        size_t    count  = cchunk->size_padded * cchunk->nda;
        int       drank  = cchunk->dest_rank;

        // add the a number of data to the destination rank
        i2o_count_[drank] = count;
        i2o_disp_[drank]  = cchunk->data - sendData;

        // do some sanity checks to make sure everything is coherent and that the chunk will write at the correct memory location
        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        FLUPS_CHECK(drank < sub_size, "Destination rank of the chunk should be inside the subcomm %d vs %d", drank, sub_size);
    }
    for (int ic = 0; ic < o2i_nchunks_; ++ic) {
        MemChunk *cchunk = o2i_chunks_ + ic;
        size_t    count  = cchunk->size_padded * cchunk->nda;
        int       drank  = cchunk->dest_rank;

        o2i_count_[drank] = count;
        o2i_disp_[drank]  = (cchunk->data - recvData);

        // again here we check that the chunk will write at the correct memory location
        FLUPS_CHECK(count < std::numeric_limits<int>::max(), "message is too big: %ld vs %d", count, std::numeric_limits<int>::max());
        FLUPS_CHECK(drank < sub_size, "Destination rank of the chunk should be inside the subcomm %d vs %d", drank, sub_size);
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
        FLUPS_INFO("sending request to rank %d of size %d %d %d", chunk->dest_rank, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
        // copy here the chunk from the input topo to the chunk
        CopyData2Chunk(nmem_in, mem, chunk);
    };

    auto complete_recv = [=](MemChunk *chunk) {
        FLUPS_INFO("recving request from rank %d of size %d %d %d", chunk->dest_rank, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
        // shuffle the data
        if (chunk->axis != chunk->dest_axis) {
            // // Print chunk data
            // printf("Chunk   = [");
            // for (int i = 0; i < chunk->isize[0]*chunk->isize[1]*chunk->isize[2]; i++) {
            //     if (chunk->nf == 2) printf("(% 02.2f, % 02.2f) ", chunk->data[2*i], chunk->data[2*i+1]);
            //     else                printf("(% 02.2f,   0.00) ", chunk->data[i]);
            // }
            // printf("]\n");
            // Print chunk data - post shuffle
            // printf("Shuffle = [");
            // for (int i = 0; i < chunk->isize[0]*chunk->isize[1]*chunk->isize[2]; i++) {
            //     if (chunk->nf == 2) printf("(% 02.2f, % 02.2f) ", chunk->data[2*i], chunk->data[2*i+1]);
            //     else                printf("(% 02.2f,   0.00) ", chunk->data[i]);
            // }
            // printf("]\n");
            // Print mem data
            // printf("Mem     = [");
            // for (int i = 0; i < nmem_out[0]*nmem_out[1]*nmem_out[2]; i++) {
            //     if (chunk->nf == 2) printf("(% 02.2f, % 02.2f) ", mem[2*i], mem[2*i+1]);
            //     else                printf("(% 02.2f,   0.00) ", mem[i]);
            // }
            // // printf("]\n

            // The chunk arrives in its dest topology and must be shuffled into its correct topology
            int iaxis[3] = {chunk->dest_axis, (chunk->dest_axis + 1) % 3, (chunk->dest_axis + 2) % 3};
            int oaxis[3] = {chunk->axis, (chunk->axis + 1) % 3, (chunk->axis + 2) % 3};

            int n0_in = 1, n1_in = 1;
            int n0_out = 1, n1_out = 1;
            for (int i = 0; i < 3; i++) {
                if (i == chunk->dest_axis)  n0_in *= chunk->isize[iaxis[i]];
                else                        n1_in *= chunk->isize[iaxis[i]];
                if (i == chunk->axis)       n0_out *= chunk->isize[oaxis[i]];
                else                        n1_out *= chunk->isize[oaxis[i]];
            }

            double alpha[2] = {1., 0.};

            // Call to tmp THEN copy...
            cblas_zomatcopy(CblasColMajor, CblasTrans, chunk->n1, chunk->n0, alpha, chunk->data, chunk->n1, chunk->shuffle_buff, chunk->n0);
            memcpy(chunk->data, chunk->shuffle_buff, chunk->size_padded*sizeof(double));

            // Print topo & chunk info
            // chunk->isize is supposed to be the mem size in the "input topology", which is not the case ??
            // chunk->axis is the axis in the input topology ( = the topo_out here)
            printf("Shuffling from axis %d to axis %d (chunk info)\n", chunk->dest_axis, chunk->axis);
            printf("                    %d         %d (topo  info)\n", topo_in->axis(), topo_out->axis());
            printf("chunk xyz = (%d, %d, %d)\n", chunk->isize[0], chunk->isize[1], chunk->isize[2]);
            printf("chunk IN  = (%d, %d, %d)\n", chunk->isize[iaxis[0]], chunk->isize[iaxis[1]], chunk->isize[iaxis[2]]);
            printf("chunk OUT = (%d, %d, %d)\n", chunk->isize[oaxis[0]], chunk->isize[oaxis[1]], chunk->isize[oaxis[2]]);
            printf("mem       = (%d, %d, %d)\n", nmem_out[oaxis[0]], nmem_out[oaxis[1]], nmem_out[oaxis[2]]);
            printf("PLAN size: (%d, %d)\n", chunk->n0, chunk->n1);
            printf("IN   size: (%d, %d)\n", n0_in, n1_in);
            printf("OUT  size: (%d, %d)\n", n0_out, n1_out);
        }
        printf("Copyin' data\n");
        CopyChunk2Data(chunk, nmem_out, mem);
        printf("\n");

    };
    //..........................................................................
    // Prepare the send buffer
    {
        m_profStarti(prof, "copy data 2 chunk");
        int count = 0;
        for (int ir = 0; ir < sub_size; ++ir) {
            // we might have nothing to send to that rank in the subcomm
            if (count_send[ir] > 0) {
                MemChunk *cchunk = send_chunks + count;
                FLUPS_CHECK(cchunk->dest_rank == ir, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ir);
                set_sendbuf(cchunk);
                count += 1;
            }
        }
        m_profStopi(prof, "copy data 2 chunk");
    }

    m_profStarti(prof, "all2all");
    MPI_Alltoallv(send_buf, count_send, disp_send, MPI_DOUBLE, recv_buf, count_recv, disp_recv, MPI_DOUBLE, subcomm);
    m_profStopi(prof, "all2all");

    // reset the memory to 0.0 as we do inplace computations
    const size_t reset_size = topo_out->memsize();
    std::memset(mem, 0, reset_size * sizeof(double));

    // Copy back the recveived data
    {
        m_profStarti(prof, "shuffle and copy chunk 2 data");
        int count = 0;
        for (int ir = 0; ir < sub_size; ++ir) {
            if (count_recv[ir] > 0) {
                MemChunk *cchunk = recv_chunks + count;
                FLUPS_CHECK(cchunk->dest_rank == ir, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ir);
                complete_recv(cchunk);
                count += 1;
            }
        }
        m_profStopi(prof, "shuffle and copy chunk 2 data");
    }

    //--------------------------------------------------------------------------
    END_FUNC;
}
