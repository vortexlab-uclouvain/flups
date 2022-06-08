#include "SwitchTopoX.hpp"

#include "h3lpr/macros.hpp"

using namespace std;

SwitchTopoX::SwitchTopoX(const Topology* topo_in, const Topology* topo_out, const int shift[3], H3LPR::Profiler* prof)
    : i2o_shift_{shift[0], shift[1], shift[2]}, o2i_shift_{-shift[0], -shift[1], -shift[2]}, topo_in_(topo_in), topo_out_(topo_out), inComm_(topo_in->get_comm()), outComm_(topo_out->get_comm()), prof_(prof) {
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
    MPI_Comm_compare(topo_in->get_comm(), topo_out->get_comm(), &comp);
    FLUPS_CHECK(comp == MPI_IDENT, "we do NOT support different communicators in and out for the moment");
#endif
    idswitchtopo_ = topo_out->axproc(topo_out->axis());
    //--------------------------------------------------------------------------
    END_FUNC;
}

SwitchTopoX::~SwitchTopoX() {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // deallocate the plans
    for (int ic = 0; ic < i2o_nchunks_; ic++) {
        // the shuffle happens in the "out" topology
        MPI_Type_free(&i2o_chunks_[ic].dtype);
        fftw_destroy_plan(i2o_chunks_[ic].shuffle);
    }
    for (int ic = 0; ic < o2i_nchunks_; ic++) {
        // the shuffle happens in the "in" topology
        MPI_Type_free(&o2i_chunks_[ic].dtype);
        fftw_destroy_plan(o2i_chunks_[ic].shuffle);
    }

    // free the MemChunks
    m_free(i2o_chunks_);
    m_free(o2i_chunks_);

    // deallocate the subcom, only if it's NOT the inComm
    // MPI_Comm_free(&subcomm_);
    int comp;
    MPI_Comm_compare(subcomm_,inComm_,&comp);
    if(comp!=MPI_IDENT){
        MPI_Comm_free(&subcomm_);
    }

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
    // The input topo may have been reset to real, even if this switchtopo is a complex2complex.
    // We create a tmp input topo which is complex if needed, for the computation of start and end.
    int tmp_nglob[3], tmp_nproc[3], tmp_axproc[3];
    for (int i = 0; i < 3; i++){
        tmp_nglob[i] = topo_in_->nglob(i);
        tmp_nproc[i] = topo_in_->nproc(i);
        tmp_axproc[i] = topo_in_->axproc(i);
    }
    
    Topology * topo_in_tmp = new Topology(topo_in_->axis(), topo_in_->lda(), tmp_nglob, tmp_nproc, topo_in_->isComplex(), tmp_axproc, FLUPS_ALIGNMENT, topo_in_->get_comm());
    // If the output topo is complex while the input topo is real, switch the input topo as a complex one 
    if(topo_out_->isComplex() && !topo_in_->isComplex()){
        topo_in_tmp->switch2complex();
    }
    
    // Populate the arrays of memory chunks
    FLUPS_INFO("I2O chunks");
    PopulateChunk(i2o_shift_, topo_in_tmp, topo_out_, &i2o_nchunks_, &i2o_chunks_, &i2o_selfcomm_);
    FLUPS_INFO("O2I chunks");
    PopulateChunk(o2i_shift_, topo_out_, topo_in_tmp, &o2i_nchunks_, &o2i_chunks_, &o2i_selfcomm_);

    delete(topo_in_tmp);

    // Split the communication according to the destination of each chunk in the MPI_COMM_WORLD
    SubCom_SplitComm();
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Set the up buffers object and plan the fftw shuffle 
 * 
 * @param sendData 
 * @param recvData 
 */
void SwitchTopoX::setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData){
    BEGIN_FUNC;
    //..........................................................................
    int sub_rank; 
    MPI_Comm_rank(subcomm_, &sub_rank);

    send_buf_ = sendData; //reinterpret_cast<opt_double_ptr>(m_calloc(send_buff_size * sizeof(double)));
    recv_buf_ = recvData; //reinterpret_cast<opt_double_ptr>(m_calloc(recv_buff_size * sizeof(double)));

    // assign the chunks with the relevant memory address
    size_t size_counter = 0;
    for (int ic = 0; ic < i2o_nchunks_; ic++) {
        // The self communication is the last one in the chunk array but correspond to its 
        // rank number in the physical memory space. Therefore we need to check if the chunk index is 
        // the same as the current rank number in the subcommunicator. If its the case, we have to register the 
        // pointer to the buffer in the last chunk 
        // Once we have perform that for the self communication, we must continue as usual. Yet, there is 
        // a difference of 1 between ic and the real chunk index
        bool is_self_idx = ((i2o_selfcomm_>=0) && (ic == sub_rank));
        bool is_greater_self_idx = ((i2o_selfcomm_>=0) && (ic > sub_rank));
        int i_chunk = is_self_idx ? i2o_selfcomm_ : ic - (int) is_greater_self_idx; 
        
        FLUPS_CHECK(i2o_chunks_[i_chunk].size_padded > 0 && i2o_chunks_[i_chunk].nda > 0, "the size of the chunk cannot be null here");
        i2o_chunks_[i_chunk].data = send_buf_ + size_counter;
        // update the counter
        size_counter += i2o_chunks_[i_chunk].size_padded * i2o_chunks_[i_chunk].nda;
    }
    size_counter = 0;
    for (int ic = 0; ic < o2i_nchunks_; ic++) {
        bool is_self_idx = ((o2i_selfcomm_>=0) && (ic == sub_rank));
        bool is_greater_self_idx = ((o2i_selfcomm_>=0) && (ic > sub_rank));
        int i_chunk = is_self_idx ? o2i_selfcomm_ : ic - (int) is_greater_self_idx; 

        FLUPS_CHECK(o2i_chunks_[i_chunk].size_padded > 0 && o2i_chunks_[i_chunk].nda > 0, "the size of the chunk cannot be null here");
        o2i_chunks_[i_chunk].data = recv_buf_ + size_counter;
        size_counter += o2i_chunks_[i_chunk].size_padded * o2i_chunks_[i_chunk].nda;
    }

    FLUPS_INFO("memory addresses have been assigned");

    //..........................................................................
    // compute the subcomm and start the gather of each rank
    // SubCom_SplitComm();
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
size_t SwitchTopoX::get_ChunkArraysMemSize(const int lda, const int nchunks, const MemChunk* chunks) const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // nultiply by the number of blocks
    size_t total = 0;
    for (int ib = 0; ib < nchunks; ib++) {
        FLUPS_CHECK(chunks[ib].nda == lda, "The number of component given must be equal to the one of the chunks --> %d vs %d", chunks[ib].nda, lda );
        total += chunks[ib].size_padded;
    }
    return total * lda;
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
size_t SwitchTopoX::get_bufMemSize() const {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // get the buffer sizes: send lives in the input topo, recv in the output one
    size_t send_buff_size = get_ChunkArraysMemSize(topo_in_->lda(), i2o_nchunks_, i2o_chunks_);
    size_t recv_buff_size = get_ChunkArraysMemSize(topo_out_->lda(), o2i_nchunks_, o2i_chunks_);
    //--------------------------------------------------------------------------
    END_FUNC;
    
    return std::max(send_buff_size, recv_buff_size);
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
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // get my rank and use-it as the initial color
    FLUPS_INFO("Get my rank");
    int comm_size, rank;
    MPI_Comm_rank(inComm_, &rank);
    MPI_Comm_size(inComm_, &comm_size);

    // // the easier way to get the color of subcomm is to use the id of the non-transposed direction
    // const int in_axis    = topo_in_->axis();
    // const int out_axis   = topo_out_->axis();
    // const int other_axis = (3) - (in_axis + out_axis);
    // // const int color      = topo_in_->rankd(other_axis);
    // // const int key        = topo_in_->rankd(in_axis) + topo_in_->rankd(out_axis);
    // // MPI_Comm_split(inComm_, color, key, &subcomm_);

    // const int period[3] = {0, 0, 0};
    // const int nproc[3]  = {topo_in_->nproc(0), topo_in_->nproc(1), topo_in_->nproc(2)};
    // MPI_Comm  comm_temp;
    // MPI_Cart_create(inComm_, 3, nproc, period, true, &comm_temp);
    // // false if we split along that direction
    // const int remain_dim[3] = {other_axis != 0, other_axis != 1, other_axis != 2};
    // MPI_Cart_sub(comm_temp, remain_dim, &subcomm_);

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
    // after this I have a clear view on who is communicating with me
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
        n_wrongColor += inMyGroup[ir];
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
        // subcomm_ = inComm_;
        MPI_Comm_dup(inComm_,&subcomm_);
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
            // // avoids the creation of a communicator
            // subcomm_ = inComm_;
            MPI_Comm_dup(inComm_,&subcomm_);
            FLUPS_INFO("I did not create a new comm since I did not find a way to subdivise the initial comm");
        } else {
            // create the communicator and give a name
            MPI_Comm_split(inComm_, mycolor, rank, &subcomm_);
            std::string commname = "comm-" + std::to_string(mycolor);
            MPI_Comm_set_name(subcomm_, commname.c_str());

#if (0 == FLUPS_OLD_MPI)
            // apply some fancy parameters to allow faster MPI calls if we have a MPI-4.0 compliant version
            MPI_Info info;
            MPI_Info_create(&info);
            // MPI_Info_set(info, "mpi_assert_exact_length", "true");
            MPI_Info_set(info, "mpi_assert_allow_overtaking", "true");
            MPI_Info_set(info, "mpi_assert_no_any_tag", "true");
            MPI_Info_set(info, "mpi_assert_no_any_source", "true");
            MPI_Comm_set_info(subcomm_, info);
            MPI_Info_free(&info);
#endif
        }
    }
    // free the vectors
    m_free(colors);
    m_free(inMyGroup);
    //--------------------------------------------------------------------------
    END_FUNC;
}


void SwitchTopoX::print_info() const {
    BEGIN_FUNC;
    FLUPS_CHECK(topo_in_ != NULL, "You must initialise the switchtopo before printing the information");
    //--------------------------------------------------------------------------
    // Get info about the in communicator 
    int rank_world, size_world;
    MPI_Comm_rank(inComm_, &rank_world);
    MPI_Comm_size(inComm_, &size_world);

    // Get info about the subcommunicator 
    int rank_subcomm, size_subcomm, num_subcomm;
    MPI_Comm_rank(subcomm_, &rank_subcomm);
    MPI_Comm_size(subcomm_, &size_subcomm);

    // Initiate the global information
    int max_size_subcomm, min_size_subcomm;
    int ttl_ni2ochunks, max_ni2ochunks, min_ni2ochunks;
    int ttl_no2ichunks, max_no2ichunks, min_no2ichunks;
    size_t ttl_size_i2ochunks;
    size_t ttl_size_o2ichunks;

    // get the local information 
    size_t loc_ttl_size_i2ochunk = 0;
    size_t loc_ttl_size_o2ichunk = 0;

    // define the MPI structure used to get the location of the max and the min
    struct{
        size_t val = 0;
        int   rank = 0;
    } max_size_i2ochunks, max_size_o2ichunks,\
      loc_max_size_i2ochunk, loc_max_size_o2ichunk;

    struct{
        size_t val = std::numeric_limits<size_t>::max();
        int   rank = 0;
    }  min_size_i2ochunks,min_size_o2ichunks,\
      loc_min_size_i2ochunk, loc_min_size_o2ichunk;

    for (int ir = 0; ir < i2o_nchunks_; ++ir) {
        MemChunk      *cchunk = i2o_chunks_ + ir;
        loc_ttl_size_i2ochunk += cchunk->size_padded * cchunk->nda;
        loc_max_size_i2ochunk.val  =  std::max((cchunk->size_padded * cchunk->nda), loc_max_size_i2ochunk.val);
        loc_min_size_i2ochunk.val  =  std::min((cchunk->size_padded * cchunk->nda), loc_min_size_i2ochunk.val);
    }
    loc_max_size_i2ochunk.rank = rank_world;
    loc_min_size_i2ochunk.rank = rank_world;

    
    for (int ir = 0; ir < o2i_nchunks_; ++ir) {
        MemChunk      *cchunk = o2i_chunks_ + ir;
        loc_ttl_size_o2ichunk += cchunk->size_padded * cchunk->nda;
        loc_max_size_o2ichunk.val  = std::max((cchunk->size_padded * cchunk->nda), loc_max_size_o2ichunk.val);
        loc_min_size_o2ichunk.val  = std::min((cchunk->size_padded * cchunk->nda), loc_min_size_o2ichunk.val);
    }
    loc_max_size_o2ichunk.rank = rank_world;
    loc_min_size_o2ichunk.rank = rank_world;

    // Get all the global information
    {   
        // Get the number of rank communicator inside the sub_communicator
        // To do that, we compute the number of rank who are identified as root (rank 0)
        // in their communicator. The sum of number of root rank is equal to the number of communicator
        // /!\ Only the root rank of inComm have its num_subcomm up-to-date
        int is_root = (int)(rank_subcomm == 0);
        MPI_Reduce(&is_root, &num_subcomm, 1, MPI_INT, MPI_SUM, 0, inComm_);
        MPI_Reduce(&size_subcomm, &max_size_subcomm, 1, MPI_INT, MPI_MAX, 0, inComm_);
        MPI_Reduce(&size_subcomm, &min_size_subcomm, 1, MPI_INT, MPI_MIN, 0, inComm_);

        // Get information about nchunks per rank 
        MPI_Reduce(&i2o_nchunks_, &ttl_ni2ochunks, 1, MPI_INT, MPI_SUM, 0, inComm_);
        MPI_Reduce(&i2o_nchunks_, &max_ni2ochunks, 1, MPI_INT, MPI_MAX, 0, inComm_);
        MPI_Reduce(&i2o_nchunks_, &min_ni2ochunks, 1, MPI_INT, MPI_MIN, 0, inComm_);
        
        MPI_Reduce(&o2i_nchunks_, &ttl_no2ichunks, 1, MPI_INT, MPI_SUM, 0, inComm_);
        MPI_Reduce(&o2i_nchunks_, &max_no2ichunks, 1, MPI_INT, MPI_MAX, 0, inComm_);
        MPI_Reduce(&o2i_nchunks_, &min_no2ichunks, 1, MPI_INT, MPI_MIN, 0, inComm_);

        MPI_Reduce(&loc_ttl_size_i2ochunk, &ttl_size_i2ochunks, 1, MPI_LONG, MPI_SUM, 0, inComm_);
        MPI_Reduce(&loc_max_size_i2ochunk, &max_size_i2ochunks, 1, MPI_LONG_INT, MPI_MAXLOC, 0, inComm_);
        MPI_Reduce(&loc_min_size_i2ochunk, &min_size_i2ochunks, 1, MPI_LONG_INT, MPI_MINLOC, 0, inComm_);

        MPI_Reduce(&loc_ttl_size_o2ichunk, &ttl_size_o2ichunks, 1, MPI_LONG, MPI_SUM, 0, inComm_);
        MPI_Reduce(&loc_max_size_o2ichunk, &max_size_o2ichunks, 1, MPI_LONG_INT, MPI_MAXLOC, 0, inComm_);
        MPI_Reduce(&loc_min_size_o2ichunk, &min_size_o2ichunks, 1, MPI_LONG_INT, MPI_MINLOC, 0, inComm_);
    }
    
    if(rank_world == 0){
        // Create the prof directory if it does not exist
        {
            struct stat st = {0};
            if (stat("./prof", &st) == -1) {
                mkdir("./prof", 0700);
            }
        }
        // Get the file name
        std::string name = "./prof/SwitchTopo_" + std::to_string(idswitchtopo_) + "_info.txt";

        // Open the file and print general information about the subcommunicators
        FILE *file = fopen(name.c_str(), "a+");

        // Print the information
        fprintf(file,"Switchtopo %d - from axis in = %d to axis out = %d -- total number of ranks = %d\n", idswitchtopo_, topo_in_->axis(), topo_out_->axis(), size_world);
        fprintf(file, "-----------------------------------------------------------------------------------------------------\n");
        fprintf(file,"Topo in: lda   = %d -- iscomplex = %d\n", topo_in_->lda(), topo_in_->isComplex());
        fprintf(file,"         nglob = %d %d %d \n", topo_in_->nglob(0), topo_in_->nglob(1), topo_in_->nglob(2));
        fprintf(file,"         nproc = %d %d %d \n", topo_in_->nproc(0), topo_in_->nproc(1), topo_in_->nproc(2));
        fprintf(file,"Topo out: lda   = %d -- iscomplex = %d\n", topo_out_->lda(), topo_out_->isComplex());
        fprintf(file,"         nglob = %d %d %d \n", topo_out_->nglob(0), topo_out_->nglob(1), topo_out_->nglob(2));
        fprintf(file,"         nproc = %d %d %d \n", topo_out_->nproc(0), topo_out_->nproc(1), topo_out_->nproc(2));
        fprintf(file, "-----------------------------------------------------------------------------------------------------\n");
        fprintf(file, "number of subcom               = %d \n", num_subcomm);
        fprintf(file, "mean number of rank per subcom = %d \n", size_world/num_subcomm);
        fprintf(file, "max  number of rank per subcom = %d \n", max_size_subcomm);
        fprintf(file, "min  number of rank per subcom = %d \n", min_size_subcomm);
        fprintf(file, "-----------------------------------------------------------------------------------------------------\n");
        fprintf(file, "total  number of chunks in the input topo = %d \n", ttl_ni2ochunks);
        fprintf(file, "mean  number of chunks in the input topo = %d \n", ttl_ni2ochunks/size_world);
        fprintf(file, "max   number of chunks in the input topo = %d \n", max_ni2ochunks);
        fprintf(file, "min   number of chunks in the input topo = %d \n", min_ni2ochunks);
        fprintf(file, "\n");
        fprintf(file, "mean  chunk size in the input topo = %lu \n", ttl_size_i2ochunks/ttl_ni2ochunks);
        fprintf(file, "max   chunk size in the input topo = %lu belongs to rank %d in Comm_WORLD \n", max_size_i2ochunks.val, max_size_i2ochunks.rank);
        fprintf(file, "min   chunk size in the input topo = %lu belongs to rank %d in Comm_WORLD \n", min_size_i2ochunks.val, min_size_i2ochunks.rank);
        fprintf(file, "-----------------------------------------------------------------------------------------------------\n");
        fprintf(file, "total  number of chunks in the output topo = %d \n", ttl_no2ichunks);
        fprintf(file, "mean  number of chunks in the output topo = %d \n", ttl_no2ichunks/size_world);
        fprintf(file, "max   number of chunks in the output topo = %d \n", max_no2ichunks);
        fprintf(file, "min   number of chunks in the output topo = %d \n", min_no2ichunks);
        fprintf(file, "\n");
        fprintf(file, "mean  chunk size in the output topo = %lu \n", ttl_size_o2ichunks/ttl_no2ichunks);
        fprintf(file, "max   chunk size in the output topo = %lu belongs to rank %d in Comm_WORLD \n", max_size_o2ichunks.val, max_size_o2ichunks.rank);
        fprintf(file, "min   chunk size in the output topo = %lu belongs to rank %d in Comm_WORLD \n", min_size_o2ichunks.val, min_size_o2ichunks.rank);
        fclose(file);
    }
    
    MPI_Barrier(inComm_);
    //--------------------------------------------------------------------------
    END_FUNC;
}