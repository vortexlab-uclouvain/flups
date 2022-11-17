/**
 * @file SwitchTopoX_a2a.cpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/
#include "SwitchTopoX_a2a.hpp"

#if FLUPS_OLD_MPI
void All2Allv(MemChunk *send_chunks, const int *count_send, const int *disp_send,
              MemChunk *recv_chunks, const int *count_recv, const int *disp_recv, 
              opt_double_ptr send_buf, opt_double_ptr recv_buf, MPI_Request* all2all_rqst, MPI_Comm subcomm,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof);
#else 

#include <mpi-ext.h>

void All2Allv(MemChunk *send_chunks, const int *count_send, 
             MemChunk *recv_chunks, const int *count_recv, 
             MPI_Request * all2all_rqst, MPI_Comm subcomm,
             const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof);
#endif

void PrintCountArr(const std::string filename, const int* count_arr, int array_size, MPI_Comm incomm);

void PrintCountArr(const std::string filename, const int* count_arr, int array_size, MPI_Comm incomm);

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
    i2o_rqst_  = reinterpret_cast<MPI_Request *>(m_calloc(1 * sizeof(MPI_Request)));    
    i2o_rqst_[0] = MPI_REQUEST_NULL;
    
    o2i_rqst_ = reinterpret_cast<MPI_Request *>(m_calloc(1 * sizeof(MPI_Request)));
    o2i_rqst_[0] = MPI_REQUEST_NULL;

#if FLUPS_OLD_MPI == 0     
    MPI_Info i2o_info;
    MPI_Info_create(&i2o_info);
    MPIX_Alltoallv_init(send_buf_, i2o_count_, i2o_disp_, MPI_DOUBLE, recv_buf_, o2i_count_, o2i_disp_, MPI_DOUBLE, subcomm_, i2o_info, i2o_rqst_);
    MPI_Info_free(&i2o_info);
    MPI_Info o2i_info;
    MPI_Info_create(&o2i_info);
    MPIX_Alltoallv_init(recv_buf_, o2i_count_, o2i_disp_, MPI_DOUBLE, send_buf_, i2o_count_, i2o_disp_, MPI_DOUBLE, subcomm_, o2i_info, o2i_rqst_);
    MPI_Info_free(&o2i_info);
#endif 

    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX_a2a::~SwitchTopoX_a2a(){
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // free the request arrays
    if(i2o_rqst_ != NULL){
        if (i2o_rqst_[0] != MPI_REQUEST_NULL) {
            FLUPS_INFO("Freeing i2o rqst");
            MPI_Request_free(i2o_rqst_);
        }
        m_free(i2o_rqst_);
    }

    if(o2i_rqst_ != NULL){
        if (o2i_rqst_[0] != MPI_REQUEST_NULL) {
            FLUPS_INFO("Freeing o2i rqst");
            MPI_Request_free(o2i_rqst_);
        }
        m_free(o2i_rqst_);
    }

    // Free the alltoall array
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
    m_profStarti(prof_, "Switchtopo%d_%s", idswitchtopo_, (FLUPS_FORWARD == sign) ? "forward" : "backward");

#if FLUPS_OLD_MPI    
    if (sign == FLUPS_FORWARD) { 
        All2Allv(i2o_chunks_, i2o_count_, i2o_disp_,
                 o2i_chunks_, o2i_count_, o2i_disp_,
                 send_buf_, recv_buf_, i2o_rqst_, subcomm_,
                 topo_in_, topo_out_, v, prof_);
    } else {
        All2Allv(o2i_chunks_, o2i_count_, o2i_disp_,
                 i2o_chunks_, i2o_count_, i2o_disp_,
                 recv_buf_, send_buf_, o2i_rqst_, subcomm_,
                 topo_out_, topo_in_, v, prof_);
    }
#else

    if (sign == FLUPS_FORWARD) { 
        All2Allv(i2o_chunks_, i2o_count_, 
                 o2i_chunks_, o2i_count_,
                 i2o_rqst_, subcomm_,
                 topo_in_, topo_out_, v, prof_);
    } else {
        All2Allv(o2i_chunks_, o2i_count_,
                 i2o_chunks_, i2o_count_, 
                 o2i_rqst_, subcomm_,
                 topo_out_, topo_in_, v, prof_);
    } 
#endif

    m_profStopi(prof_, "Switchtopo%d_%s", idswitchtopo_, (FLUPS_FORWARD == sign) ? "forward" : "backward");
    //--------------------------------------------------------------------------
    END_FUNC;
}

void SwitchTopoX_a2a::disp() const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    END_FUNC;
}


#if FLUPS_OLD_MPI
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
              opt_double_ptr send_buf, opt_double_ptr recv_buf, MPI_Request* all2all_rqst, MPI_Comm subcomm,
              const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof) {
#else 
void All2Allv(MemChunk *send_chunks, const int *count_send, 
             MemChunk *recv_chunks, const int *count_recv, 
             MPI_Request * all2all_rqst, MPI_Comm subcomm,
             const Topology *topo_in, const Topology *topo_out, opt_double_ptr mem, H3LPR::Profiler* prof) {
#endif

    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int nmem_in[3]  = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
    const int nmem_out[3] = {topo_out->nmem(0), topo_out->nmem(1), topo_out->nmem(2)};

    int sub_size;
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
        DoShuffleChunk(chunk);
        CopyChunk2Data(chunk, nmem_out, mem);
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
                // FLUPS_CHECK(cchunk->dest_rank == ir, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ir);
                set_sendbuf(cchunk);
                count += 1;
            }
        }
        m_profStopi(prof, "copy data 2 chunk");
    }


    m_profStarti(prof, "all2all - start");
#if FLUPS_OLD_MPI

    MPI_Ialltoallv(send_buf, count_send, disp_send, MPI_DOUBLE, recv_buf, count_recv, disp_recv, MPI_DOUBLE, subcomm, all2all_rqst);

#else 
    MPI_Start(all2all_rqst);
#endif
    m_profStopi(prof, "all2all - start");

    // reset the memory to 0.0 as we do inplace computations
    const size_t reset_size = topo_out->memsize();
    std::memset(mem, 0, reset_size * sizeof(double));

    m_profStarti(prof, "all2all - wait");
    MPI_Wait(all2all_rqst, MPI_STATUS_IGNORE);    
    m_profStopi(prof, "all2all - wait");

    // Copy back the recveived data
    {
        m_profStarti(prof, "shuffle and copy chunk 2 data");
        int count = 0;
        for (int ir = 0; ir < sub_size; ++ir) {
            if (count_recv[ir] > 0) {
                MemChunk *cchunk = recv_chunks + count;
                // FLUPS_CHECK(cchunk->dest_rank == ir, "Destination rank of the chunk should correspond to the chunk indexing %d vs %d", cchunk->dest_rank, ir);
                complete_recv(cchunk);
                count += 1;
            }
        }
        m_profStopi(prof, "shuffle and copy chunk 2 data");
    }

    //--------------------------------------------------------------------------
    END_FUNC;
}


void SwitchTopoX_a2a::print_info() const {
    BEGIN_FUNC;
    FLUPS_CHECK(i2o_count_!=NULL, "The setup must be initialised before printring their information");
    //--------------------------------------------------------------------------
    // Call the base class
    SwitchTopoX::print_info();

    // Get MPI information 
    int world_size;
    MPI_Comm_size(inComm_, &world_size);
    int sub_size;
    MPI_Comm_size(subcomm_, &sub_size);

    // Get filename  
    std::string filename_forward = "./prof/Nrank_" + std::to_string(world_size)+ "_SwitchTopo_" + std::to_string(idswitchtopo_) + "_forward_messages.txt";
    PrintCountArr(filename_forward, i2o_count_, sub_size, inComm_);
    
    std::string filename_backward = "./prof/Nrank_" + std::to_string(world_size)+ "_SwitchTopo_" + std::to_string(idswitchtopo_) + "_backward_messages.txt";
    PrintCountArr(filename_backward, o2i_count_, sub_size, inComm_);
    //--------------------------------------------------------------------------
    END_FUNC;

}


void PrintCountArr(const std::string filename, const int* count_arr, int array_size, MPI_Comm incomm){
    BEGIN_FUNC;
    FLUPS_CHECK(count_arr!=NULL, "The setup must be initialised before printring their information");
    //--------------------------------------------------------------------------
    // Get MPI information 
    int world_rank;
    MPI_Comm_rank(incomm, &world_rank);

    // Create the message for this rank
    std::string msg = std::string((world_rank == 0) ? "" : "\n") + "rank " + std::to_string(world_rank);
    for(int i = 0; i < array_size; ++i){
        msg += " " +  std::to_string((count_arr[i]/1000.0)*sizeof(double));
    }
    size_t msg_size = msg.length();


    // Open the file using MPI IO functions 
    FLUPS_INFO("opening file <%s>", filename.c_str());

    MPI_File   mpi_file;     
    int err = MPI_File_open(incomm, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL, MPI_INFO_NULL, &mpi_file);
    FLUPS_CHECK(err == MPI_SUCCESS, "ERROR while opening  <%s>, MPI_File_open failed (error = %d)", filename.c_str(), err);

    // Get the current offset
    size_t offset_ttl = 0; 
    FLUPS_CHECK(sizeof(unsigned long) == sizeof(size_t), "The mpi type should correspond to the standard size");
    MPI_Scan(&msg_size, &offset_ttl, 1, MPI_UNSIGNED_LONG, MPI_SUM, incomm);

    // the current position of current proc
    //MPI_Offset offset = offset_ttl - msg_size;
    //MPI_File_seek(mpi_file, offset, MPI_SEEK_SET);

    MPI_Status status;
    MPI_Offset offset = offset_ttl - msg_size;
    err = MPI_File_write_at(mpi_file, offset, msg.c_str(), msg_size, MPI_CHAR, &status);
    FLUPS_CHECK(err == MPI_SUCCESS, "ERROR while writing in  <%s> (error = %d)", filename.c_str(), err);

    MPI_File_close(&mpi_file);
    //--------------------------------------------------------------------------
    END_FUNC;
}
