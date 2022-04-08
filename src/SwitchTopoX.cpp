#include "SwitchTopoX.hpp"

#include "h3lpr/macros.hpp"

using namespace std;

SwitchTopoX::SwitchTopoX(const int shift[3], Topology* topo_in, Topology* topo_out, H3LPR::Profiler* prof)
    : i2o_shift_{shift[0], shift[1], shift[2]}, o2i_shift_{-shift[0], -shift[1], -shift[2]}, topo_in_(topo_in), topo_out_(topo_out) {
    BEGIN_FUNC;
    FLUPS_CHECK(topo_in->nf() == topo_out->nf(), "The two topos must both be complex or both real");
    FLUPS_CHECK(topo_in->lda() == topo_out->lda(), "The two topos must both be complex or both real");
    //--------------------------------------------------------------------------
#ifndef NDEBUG
    // TG: I have removed the different inComm and outComm support as its not clear to me how to properly do that
    // also, I cannot think about any situation where it would be useful.
    // supporting it would be cool though!
    // in order to do so, we need to improse a convention on which communicator is used for the communication
    int comp;
    MPI_Comm_compare(inComm_, outComm_, &comp);
    FLUPS_CHECK(comp == MPI_IDENT, "we do NOT support different communicators in and out for the moment");
#endif
    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX::~SwitchTopoX() {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // deallocate the plans
    for (int ic = 0; ic < i2o_nchunks_; ic++) {
        // the shuffle happens in the "out" topology
        fftw_destroy_plan(i2o_chunks_[ic].shuffle);
    }
    for (int ic = 0; ic < o2i_nchunks_; ic++) {
        // the shuffle happens in the "in" topology
        fftw_destroy_plan(o2i_chunks_[ic].shuffle);
    }

    // free the buffers
    m_free(send_buf_);
    m_free(recv_buf_);

    // free the MemChunks
    m_free(i2o_chunks_);
    m_free(o2i_chunks_);

    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Setup the chunks and the subcommunicator for the communications
 *
 */
void SwitchTopoX::setup() {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    PopulateChunk(i2o_shift_, topo_in_, topo_out_, &i2o_nchunks_, i2o_chunks_);
    PopulateChunk(o2i_shift_, topo_out_, topo_in_, &o2i_nchunks_, o2i_chunks_);

    //..........................................................................
    // get the buffer sizes: send lives in the input topo, recv in the output one
    size_t send_buff_size = get_bufMemSize(topo_in_->lda(), i2o_nchunks_, i2o_chunks_);
    size_t recv_buff_size = get_bufMemSize(topo_out_->lda(), o2i_nchunks_, o2i_chunks_);

    send_buf_ = reinterpret_cast<opt_double_ptr>(m_calloc(send_buff_size * sizeof(double)));
    recv_buf_ = reinterpret_cast<opt_double_ptr>(m_calloc(send_buff_size * sizeof(double)));

    // assign the chunks with the relevant memory address
    size_t size_counter = 0;
    for (int ic = 0; ic < i2o_nchunks_; ic++) {
        FLUPS_CHECK(i2o_chunks_[ic].size_padded > 0 && i2o_chunks_[ic].nda > 0, "the size of the chunk cannot be null here");
        i2o_chunks_[ic].data = send_buf_ + size_counter;
        // update the counter
        size_counter += i2o_chunks_[ic].size_padded * i2o_chunks_[ic].nda;
    }
    size_counter = 0;
    for (int ic = 0; ic < o2i_nchunks_; ic++) {
        FLUPS_CHECK(o2i_chunks_[ic].size_padded > 0 && o2i_chunks_[ic].nda > 0, "the size of the chunk cannot be null here");
        o2i_chunks_[ic].data = recv_buf_ + size_counter;
        size_counter += o2i_chunks_[ic].size_padded * o2i_chunks_[ic].nda;
    }

    //..........................................................................
    // compute the subcomm and start the gather of each rank
    SubCom_SplitComm();

    int inrank, subrank, worldsize;
    MPI_Comm_size(inComm_, &worldsize);
    MPI_Comm_rank(subcomm_, &subrank);
    MPI_Comm_rank(inComm_, &inrank);

    // get the ranks of everybody in all communicators
    MPI_Request subrank_rqst;
    int*        subRanks = reinterpret_cast<int*>(m_calloc(worldsize * sizeof(int)));
    MPI_Iallgather(&subrank, 1, MPI_INT, subRanks, 1, MPI_INT, inComm_, &subrank_rqst);
    //..........................................................................
    // init the shuffle, this might take some time
    for (int ic = 0; ic < i2o_nchunks_; ic++) {
        // the shuffle happens in the "out" topology
        PlanShuffleChunk(topo_out_->nf() == 2, i2o_chunks_ + ic);
    }
    for (int ic = 0; ic < o2i_nchunks_; ic++) {
        // the shuffle happens in the "in" topology
        PlanShuffleChunk(topo_out_->nf() == 2, o2i_chunks_ + ic);
    }

    //..........................................................................
    // finish the rank assignement
    MPI_Status subrank_status;
    MPI_Wait(&subrank_rqst,&subrank_status);
    // replace the old ranks by the newest ones
    for (int ic = 0; ic < i2o_nchunks_; ++ic) {
        // get the destination rank in the old inComm_
        const int in_dest_rank = i2o_chunks_[ic].dest_rank;
        // replace by the one in the subcommunicator
        i2o_chunks_[ic].dest_rank = subRanks[in_dest_rank];
    }
    for (int ic = 0; ic < o2i_nchunks_; ++ic) {
        // get the destination rank in the old inComm_
        const int in_dest_rank = o2i_chunks_[ic].dest_rank;
        // replace by the one in the subcommunicator
        o2i_chunks_[ic].dest_rank = subRanks[in_dest_rank];
    }
    // free the allocated array
    m_free(subRanks);
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief returns the communication buffer as the max between i2o and o2i
 *
 * @warning The size ensures that the start of every chunk is aligned (if the correctly offset)
 *
 * @return size_t
 */
size_t SwitchTopoX::get_bufMemSize(const size_t lda, const int nchunks, const MemChunk* chunks) const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // the nf is the maximum between in and out
    const int nf = std::max(topo_in_->nf(), topo_out_->nf());
    // nultiply by the number of blocks
    size_t total = 0;
    for (int ib = 0; ib < nchunks; ib++) {
        total += chunks[ib].size_padded *  chunks[ib].nda;
    }
    return total * lda;
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief split the inComm_ communicator into subcomms
 *
 * We here find the colors of the comm, i.e. ranks communicating together have the same color.
 * Once the color are known, we divide the current communicator into subcomms.
 *
 *
 */
void SwitchTopoX::SubCom_SplitComm() {
    //--------------------------------------------------------------------------
    BEGIN_FUNC;
    // get my rank and use-it as the initial color
    int comm_size, rank;
    MPI_Comm_rank(inComm_, &rank);
    MPI_Comm_size(inComm_, &comm_size);

    //-------------------------------------------------------------------------
    /** - Set the starting color and determine who I wish to get in my group */
    //-------------------------------------------------------------------------
    // allocate colors and inMyGroup array
    int*  colors    = (int*)m_calloc(comm_size * sizeof(int));
    bool* inMyGroup = (bool*)m_calloc(comm_size * sizeof(bool));

    for (int ir = 0; ir < comm_size; ir++) {
        inMyGroup[ir] = false;
    }
    inMyGroup[rank] = true;

    // do a first pass and give a color + who is in my group
    int mycolor = rank;
    for (int ib = 0; ib < i2o_nchunks_; ib++) {
        const int chunk_dest_rank  = i2o_chunks_[ib].dest_rank;
        mycolor                    = m_min(mycolor, chunk_dest_rank);
        inMyGroup[chunk_dest_rank] = true;
    }
    for (int ib = 0; ib < o2i_nchunks_; ib++) {
        const int chunk_dest_rank  = o2i_chunks_[ib].dest_rank;
        mycolor                    = m_min(mycolor, chunk_dest_rank);
        inMyGroup[chunk_dest_rank] = true;
    }

    //-------------------------------------------------------------------------
    /** - count how much ranks are in my group and assumes they don't have the same color as I do */
    //-------------------------------------------------------------------------
    int n_left       = 0;  // the global counter of wrong colors
    int n_wrongColor = 0;  // the local counter of wrong colors

    for (int ir = 0; ir < comm_size; ir++) {
        if (inMyGroup[ir]) {
            n_wrongColor += 1;
        }
    }
    // compute among everybody, if we need to continue
    MPI_Allreduce(&n_wrongColor, &n_left, 1, MPI_INT, MPI_SUM, inComm_);

    //-------------------------------------------------------------------------
    /** - Among everybody in group, get the minimum color.
     * The loop is used to force information to travel:
     * If 0 is in the group of 1 and 1 in the group of 2 and 2 in the group of 3, etc.
     * we need to spread the color 0 among the group
     */
    //-------------------------------------------------------------------------
    int iter = 0;
    while (n_left > 0 && iter < comm_size) {
        // gather the color info from everyone
        MPI_Allgather(&mycolor, 1, MPI_INT, colors, 1, MPI_INT, inComm_);
        // iterate on the proc
        n_wrongColor = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            // if the rank is in my group and it has not the same color
            if (inMyGroup[ir] && (colors[ir] != mycolor)) {
                // we first increment the counter flagging that one is missing
                n_wrongColor += 1;
                // then we change the color if we are able to do so....
                // remove 1 if we are able to solve the color issue <=> my color > colors[ir]
                n_wrongColor = n_wrongColor - (colors[ir] < mycolor);
                // changing the color if possible
                mycolor = m_min(mycolor, colors[ir]);
            }
        }
        // compute among everybody, if we need to continue
        MPI_Allreduce(&n_wrongColor, &n_left, 1, MPI_INT, MPI_SUM, inComm_);
        iter++;
    }
    // if we failed to create the subcom, uses the default one
    if (n_left > 0) {
        subcomm_ = inComm_;
        FLUPS_WARNING("I failed to create the subcomm: max iter reached, every group is not complete");
    } else {
        // If there is only 1 color left on all procs, it is 0, and I can still use COMM_WORLD
        int sumColor = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            sumColor += colors[ir];
        }
        // if nleft = 0 -> everybody is inside the same color = the rank = 0
        // we do not create a new comm if it is not necessary
        if (sumColor == 0) {
            // avoids the creation of a communicator
            subcomm_ = inComm_;
            FLUPS_INFO("I did not create a new comm since I did not find a way to subdivise the initial comm");
        } else {
            // create the communicator and give a name
            MPI_Comm_split(inComm_, mycolor, rank, &subcomm_);
            std::string commname = "comm-" + std::to_string(mycolor);
            MPI_Comm_set_name(subcomm_, commname.c_str());
        }
    }
    // free the vectors
    m_free(colors);
    m_free(inMyGroup);

    END_FUNC;
}

// /**
//  * @brief start the rank update 
//  * 
//  */
// void SwitchTopoX::SubCom_UpdateRanks_Start() {
//     BEGIN_FUNC;
//     //--------------------------------------------------------------------------
//     int inrank, subrank, worldsize;
//     MPI_Comm_size(inComm_, &worldsize);
//     MPI_Comm_rank(subcomm_, &subrank);
//     MPI_Comm_rank(inComm_, &inrank);

//     // get the ranks of everybody in all communicators
//     int* subRanks = (int*)m_calloc(worldsize * sizeof(int));
//     MPI_Allgather(&subrank, 1, MPI_INT, subRanks, 1, MPI_INT, inComm_);

//     // replace the old ranks by the newest ones
//     for (int ic = 0; ic < i2o_nchunks_; ++ic) {
//         // get the destination rank in the old inComm_
//         const int in_dest_rank = i2o_chunks_[ic].dest_rank;
//         // replace by the one in the subcommunicator
//         i2o_chunks_[ic].dest_rank = subRanks[in_dest_rank];
//     }
//     for (int ic = 0; ic < o2i_nchunks_; ++ic) {
//         // get the destination rank in the old inComm_
//         const int in_dest_rank = o2i_chunks_[ic].dest_rank;
//         // replace by the one in the subcommunicator
//         o2i_chunks_[ic].dest_rank = subRanks[in_dest_rank];
//     }
//     // free the allocated array
//     m_free(subRanks);
//     //--------------------------------------------------------------------------
//     END_FUNC;
// }

// /**
//  * @brief compute the start and count arrays needed for the all to all communication
//  *
//  * @param comm the communicator to use
//  * @param nBlock the number of block
//  * @param lda leading dimension of array
//  * @param blockSize the block sizes
//  * @param destRank the destination rank of each block
//  * @param count the count array
//  * @param start the start array
//  */
// void SwitchTopo::cmpt_start_and_count_(MPI_Comm comm, const int nBlock, const int lda, int* blockSize[3], int* destRank, int** count, int** start) {
//     BEGIN_FUNC;
//     const int nf = std::max(topo_in_->nf(),topo_out_->nf());
//     int size;
//     MPI_Comm_size(comm, &size);
//     // count the number of blocks to each process
//     (*count) = (int*)m_calloc(size * sizeof(int));
//     (*start) = (int*)m_calloc(size * sizeof(int));
//     std::memset((*count), 0, size * sizeof(int));
//     std::memset((*start), 0, size * sizeof(int));

//     // count the number of blocks per rank
//     for (int ib = 0; ib < nBlock; ib++) {
//         // get the size per block
//         const int blockMem = get_blockMemSize(ib,nf,blockSize);
//         (*count)[destRank[ib]] += blockMem * lda;
//     }
//     // compute the start indexes
//     if (start != NULL) {
//         (*start)[0] = 0;
//         for (int ir = 1; ir < size; ir++) {
//             (*start)[ir] = (*start)[ir - 1] + (*count)[ir - 1];
//         }
//     }

//     END_FUNC;
// }

/**
 * @brief setup the suffle plan to do the reordering of the indexes inside a block array form the axis of topo_in to topo_out
 *
 * The 3D array is split into a rectangular 2D array:
 * - the dimension of the current FRI
 * - the dimension of the targeted FRI
 *
 * The last dimension is aggregated with eather the current FRI or the targeted, whichever comes on its left
 * e.g.:
 * - if the shuffle is between 2 and 1, the default order is then 2 0 1, hence the dimensions will be (2 * 0) x (1)
 * - if the suffle is between 1 and 2, the  default order is then 1 2 0, hence the dimensions will be (1) x (2 * 0)
 * - if the suffle is between 0 and 1, the  default order is then 1 2 0, hence the dimensions will be (0) x (1 * 2)
 *
 * @param bSize the block size
 * @param topo_in the topo_in with the current axis
 * @param topo_out the topo_out with the desired axis
 * @param data the data on which to apply the transformation
 * @param shuffle the suffle plan
 */
// void SwitchTopoX::setup_shuffle_(const int bSize[3], const Topology* topo_in, const Topology* topo_out, double* data, fftw_plan* shuffle) {
//     BEGIN_FUNC;

//     // the nf will always be the max of both topologies !!
//     const int nf = std::max(topo_in->nf(), topo_out->nf());

//     // enable the multithreading for this plan
//     fftw_plan_with_nthreads(omp_get_max_threads());

//     fftw_iodim dims[2];
//     // dim[0] = dimension of the targeted FRI (FFTW-convention)
//     dims[0].n  = 1;
//     dims[0].is = 1;
//     dims[0].os = 1;
//     // dim[1] = dimension of the current FRI (FFTW-convention)
//     dims[1].n  = 1;
//     dims[1].is = 1;
//     dims[1].os = 1;

//     int iaxis[3] = {topo_in->axis(), (topo_in->axis() + 1) % 3, (topo_in->axis() + 2) % 3};
//     int oaxis[3] = {topo_out->axis(), (topo_out->axis() + 1) % 3, (topo_out->axis() + 2) % 3};

//     // compute the size and the stride of the array
//     for (int id = 0; id < 3; id++) {
//         if (iaxis[id] != topo_out->axis()) {
//             dims[1].n  = dims[1].n * bSize[iaxis[id]];
//             dims[0].is = dims[0].is * bSize[iaxis[id]];
//         } else {
//             break;
//         }
//     }
//     for (int id = 0; id < 3; id++) {
//         if (oaxis[id] != topo_in->axis()) {
//             dims[0].n  = dims[0].n * bSize[oaxis[id]];
//             dims[1].os = dims[1].os * bSize[oaxis[id]];
//         } else {
//             break;
//         }
//     }
//     // display some info
//     FLUPS_INFO_3("shuffle: setting up the shuffle form %d to %d", topo_in->axis(), topo_out->axis());
//     FLUPS_INFO_3("shuffle: nf = %d, blocksize = %d %d %d", nf, bSize[0], bSize[1], bSize[2]);
//     FLUPS_INFO_3("shuffle: DIM 0: n = %d, is=%d, os=%d", dims[0].n, dims[0].is, dims[0].os);
//     FLUPS_INFO_3("shuffle: DIM 1: n = %d, is=%d, os=%d", dims[1].n, dims[1].is, dims[1].os);

//     // plan the real or complex plan
//     // the nf is driven by the OUT topology ALWAYS
//     if (nf == 1) {
//         *shuffle = fftw_plan_guru_r2r(0, NULL, 2, dims, data, data, NULL, FFTW_FLAG);
//         FLUPS_CHECK(*shuffle != NULL, "Plan has not been setup");
//     } else if (nf == 2) {
//         *shuffle = fftw_plan_guru_dft(0, NULL, 2, dims, (fftw_complex*)data, (fftw_complex*)data, FLUPS_FORWARD, FFTW_FLAG);
//         FLUPS_CHECK(*shuffle != NULL, "Plan has not been setup");
//     }

//     END_FUNC;
// }
