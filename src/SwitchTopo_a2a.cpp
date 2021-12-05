/**
 * @file SwitchTopo_a2a.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * 
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include "SwitchTopo.hpp"
#include "SwitchTopo_a2a.hpp"

/**
 * @brief Construct a Switch Topo object
 * 
 * Let us consider the switch from the TOPO_IN to the TOPO_OUT.
 *
 * ```
 * +------------------------------------+
 * |  TOPO_OUT  |                       |
 * |            |                       |
 * |            |  n=5                  |
 * |            |                       |
 * |            v                       |
 * |  --------> +-------------+         |
 * |    n=3     | TOPO_IN     |         |
 * |            |             |         |
 * |            |             |         |
 * |            |             |         |
 * |            +-------------+         |
 * |                                    |
 * |                                    |
 * |                                    |
 * +------------------------------------+
 * ```
 * 
 * The shift argument will then be (3,5) since we need to add (3,5) points in the topo_output
 * to reach the (0,0,0) point in the topo_input.
 * 
 * The switch between topologies works using blocks. 
 * A block is defined as a memory block on one proc that goes on another proc.
 * The size of the block will always have the same size on every process.
 * The number of block changes from one process to another.
 * Therefore we have to initialize the block structure and then use it during the execute.
 * 
 * 
 * @param topo_input the input topology
 * @param topo_output the output topology 
 * @param shift the shift is the position of the (0,0,0) of topo_input in the topo_output indexing (in XYZ-indexing)
 * @param prof the profiler to use to profile the execution of the SwitchTopo
 */

SwitchTopo_a2a::SwitchTopo_a2a(const Topology* topo_input, const Topology* topo_output, const int shift[3], Profiler* prof) {
    BEGIN_FUNC;

    FLUPS_CHECK(topo_input->isComplex() == topo_output->isComplex(), "both topologies have to be the same kind", LOCATION);
    FLUPS_CHECK(topo_input->lda() == topo_output->lda(),"Both lda's topologies must match: in=%d vs out=%d",topo_input->lda(), topo_output->lda(),LOCATION);

    topo_in_  = topo_input;
    topo_out_ = topo_output;

    // we store the communicators at init
    inComm_  = topo_in_->get_comm();
    outComm_ = topo_out_->get_comm();

#ifdef PROF    
    prof_     = prof;
    iswitch_ = topo_out_->axproc(topo_out_->axis());
#endif
    //-------------------------------------------------------------------------
    /** - compute block information: sizes, start ids, ... */
    //-------------------------------------------------------------------------
    // get the blockshift
    shift_[0] = shift[0];
    shift_[1] = shift[1];
    shift_[2] = shift[2];

    init_blockInfo_(topo_in_, topo_out_);

    //-------------------------------------------------------------------------
    /** - initialize the profiler    */
    //-------------------------------------------------------------------------
#ifdef PROF    
    if (prof_ != NULL) {
        prof_->create("reorder", "solve");

        prof_->create("switch0", "reorder");
        prof_->create("mem2buf0", "switch0");
        prof_->create("buf2mem0", "switch0");
        prof_->create("all_2_all0", "switch0");
        prof_->create("all_2_all_v0", "switch0");

        prof_->create("switch1", "reorder");
        prof_->create("mem2buf1", "switch1");
        prof_->create("buf2mem1", "switch1");
        prof_->create("all_2_all1", "switch1");
        prof_->create("all_2_all_v1", "switch1");

        prof_->create("switch2", "reorder");
        prof_->create("mem2buf2", "switch2");
        prof_->create("buf2mem2", "switch2");
        prof_->create("all_2_all2", "switch2");
        prof_->create("all_2_all_v2", "switch2");
    }
#endif
    END_FUNC;
}

/**
 * @brief initialize the communication blocks
 * 
 * First, we compute nByBlock[3], the smallest size of unknowns that goes from one proc to another.
 * This small nByBlock is the same accross each rank.
 * 
 * Then, for each of this unit block (of size nByBlock[3]), we compute their destination rank.
 * 
 * Afterwards, using the rank of those unit blocks, we try to gather them by destination ranks.
 * all the kernels blocks that have the same destination will be packed together for the communication.
 * 
 */
void SwitchTopo_a2a::init_blockInfo_(const Topology* topo_in, const Topology* topo_out){
    BEGIN_FUNC;
    
    int comm_size,ocomm_size;
    MPI_Comm_size(inComm_, &comm_size);
    MPI_Comm_size(outComm_, &ocomm_size);

    FLUPS_CHECK(ocomm_size==comm_size,"In and out communicators must have the same size.",LOCATION);

    //-------------------------------------------------------------------------
    /** - allocate temporary arrays */
    //-------------------------------------------------------------------------
    int  inBlockv[3];
    int  onBlockv[3];
    int  istart[3];
    int  ostart[3];
    int  iend[3];
    int  oend[3];
    int  nByBlock[3];

    //-------------------------------------------------------------------------
    /** - Compute intersection ids */
    //-------------------------------------------------------------------------
    //recompute start_ and end_
    topo_in->cmpt_intersect_id(shift_, topo_out, istart, iend);
    int mshift[3] = {-shift_[0], -shift_[1], -shift_[2]};
    topo_out->cmpt_intersect_id(mshift, topo_in, ostart, oend);

    //-------------------------------------------------------------------------
    /** - get the block size as the GCD of the memory among every process between send and receive */
    //-------------------------------------------------------------------------
    cmpt_nByBlock_(istart,iend,ostart,oend,nByBlock);

    cmpt_blockIndexes_(istart, iend, nByBlock, topo_in, inBlockv);
    cmpt_blockIndexes_(ostart, oend, nByBlock, topo_out, onBlockv);

    // allocate the destination ranks
    i2o_destRank_ = (int*)flups_malloc(inBlockv[0] * inBlockv[1] * inBlockv[2] * sizeof(int));
    o2i_destRank_ = (int*)flups_malloc(onBlockv[0] * onBlockv[1] * onBlockv[2] * sizeof(int));

    // get the ranks
    // shift if the root position of the topo_in in the topo_out
    cmpt_blockDestRank_(inBlockv,nByBlock,shift_,istart,topo_in,topo_out,i2o_destRank_);
    cmpt_blockDestRank_(onBlockv,nByBlock,mshift,ostart,topo_out,topo_in,o2i_destRank_);

    // try to gather blocks together if possible, rewrittes the sizes, the blockistart, the number of blocks, the ranks and the tags
    gather_blocks_(topo_in, nByBlock, istart,iend, inBlockv, iBlockSize_, iBlockiStart_, &inBlock_, &i2o_destRank_);
    gather_blocks_(topo_out, nByBlock, ostart,oend, onBlockv, oBlockSize_, oBlockiStart_, &onBlock_, &o2i_destRank_);

    END_FUNC;
}

/**
 * @brief free the blocks information: their number, their size and their source/destination
 * 
 */
void SwitchTopo_a2a::free_blockInfo_(){
    if (i2o_destRank_ != NULL) flups_free(i2o_destRank_);
    if (o2i_destRank_ != NULL) flups_free(o2i_destRank_);

    i2o_destRank_ = NULL;
    o2i_destRank_ = NULL;

    for (int id = 0; id < 3; id++) {
        if (iBlockSize_[id] != NULL) flups_free(iBlockSize_[id]);
        if (oBlockSize_[id] != NULL) flups_free(oBlockSize_[id]);
        if (iBlockiStart_[id] != NULL) flups_free(iBlockiStart_[id]);
        if (oBlockiStart_[id] != NULL) flups_free(oBlockiStart_[id]);

        iBlockSize_[id] = NULL;
        oBlockSize_[id] = NULL;
        iBlockiStart_[id] = NULL;
        oBlockiStart_[id] = NULL;
    }
}

/**
 * @brief setup the switchtopo
 * 
 */
void SwitchTopo_a2a::setup() {
    BEGIN_FUNC;

    int rank, comm_size;
    MPI_Comm inComm = topo_in_->get_comm();
    MPI_Comm outComm = topo_out_->get_comm();

    MPI_Comm_rank(inComm, &rank);
    MPI_Comm_size(inComm, &comm_size);

    //-------------------------------------------------------------------------
    /** - compare the current communicators from topo_in and topo_out.
     * If they have change, we need to recompute all the communication init
     */
    //-------------------------------------------------------------------------
    //Ensure that comms have not changed since init. Otherwise recompute the source/destination of blocks.
    int compIn, compOut;
    MPI_Comm_compare(inComm, inComm_, &compIn);
    MPI_Comm_compare(outComm, outComm_, &compOut);
    //if the graph communicator has the same numbering as the old commn we will skip the following
    if( compIn != MPI_CONGRUENT || compOut != MPI_CONGRUENT){
        if (rank == 0){
            FLUPS_WARNING("The inComm and/or outComm have changed since this switchtopo was created. I will recompute the communication scheme.",LOCATION);
        }

        inComm_ = inComm;
        outComm_ = outComm;

        //reinit the block information
        free_blockInfo_();

        //The input topo may have been reset to real, even if this switchtopo is a complex2complex. 
        //We create a tmp input topo which is complex if needed, for the computation of start and end.
        bool isC2C = topo_out_->isComplex();
        int tmp_nglob[3], tmp_nproc[3], tmp_axproc[3];
        for(int i = 0; i<3;i++){
            tmp_nglob[i] = topo_in_->nglob(i);
            tmp_nproc[i] = topo_in_->nproc(i);
            tmp_axproc[i] = topo_in_->axproc(i);
        }
        const Topology* topo_in_tmp = new Topology(topo_in_->axis(),topo_in_->lda(), tmp_nglob,tmp_nproc,isC2C,tmp_axproc,FLUPS_ALIGNMENT,topo_in_->get_comm());
        
        //recompute block info
        init_blockInfo_(topo_in_tmp,topo_out_);

        // free the temp topo
        delete(topo_in_tmp);
    }

    //-------------------------------------------------------------------------
    /** - Setup subcomm (if possible) */
    //-------------------------------------------------------------------------
    cmpt_commSplit_();
    
    // setup the dest rank, counts and starts
    setup_subComm_(inBlock_,topo_in_->lda() ,iBlockSize_, i2o_destRank_, &i2o_count_, &i2o_start_);
    setup_subComm_(onBlock_,topo_out_->lda(),oBlockSize_, o2i_destRank_, &o2i_count_, &o2i_start_);
    
    //-------------------------------------------------------------------------
    /** - determine if we are all to all */
    //-------------------------------------------------------------------------
    int subsize;
    MPI_Comm_size(subcomm_, &subsize);

    int tmp_size = i2o_count_[0];
    is_all2all_  = (tmp_size != 0);
    is_all2all_  = is_all2all_ && (tmp_size == o2i_count_[0]);
    for (int ir = 1; ir < subsize; ir++) {
        // if the count from and to the rank is the same, we can do an A2A
        is_all2all_ = is_all2all_ && (tmp_size == i2o_count_[ir]);
        is_all2all_ = is_all2all_ && (tmp_size == o2i_count_[ir]);
    }

    //-------------------------------------------------------------------------
    /** - Check that everybody is in the same communication mode*/
    //-------------------------------------------------------------------------
    // determine if every proc is in the all_to_all mode
    bool global_is_alltoall;
    MPI_Allreduce(&is_all2all_, &global_is_alltoall, 1, MPI_CXX_BOOL, MPI_LAND, subcomm_);
    // determine if at least one proc is in the all to all mode
    bool any_is_alltoall;
    MPI_Allreduce(&is_all2all_,&any_is_alltoall,1,MPI_CXX_BOOL,MPI_LOR,subcomm_);
    // generate an error if it is not compatible
    if (is_all2all_ && (!global_is_alltoall)){
        int rlen;
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(subcomm_, myname, &rlen);
        FLUPS_ERROR("communicator %s: at least one process is NOT in the all to all communication scheme",myname,LOCATION);
    }
    if((!is_all2all_) && any_is_alltoall){
        int rlen;
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(subcomm_, myname, &rlen);
        FLUPS_ERROR("communicator %s: at least one process is in the all to all communication scheme",myname,LOCATION);
    }
    
    // if we are all to all, clean the start array
    if (is_all2all_) {
        if (i2o_start_ != NULL) {
            flups_free(i2o_start_);
            i2o_start_ = NULL;
        }
        if (o2i_start_ != NULL) {
            flups_free(o2i_start_);
            o2i_start_ = NULL;
        }
    }

    //-------------------------------------------------------------------------
    /** - Display performance information if asked */
    //-------------------------------------------------------------------------    
#ifdef PERF_VERBOSE
    int rankworld;
    MPI_Comm_rank(inComm, &rankworld);
    // we display important information for the performance
    string name = "./prof/SwitchTopo_" + std::to_string(topo_in_->axis()) + "to" + std::to_string(topo_out_->axis()) + "rank_" + std::to_string(rankworld) + ".txt";
    FILE* file = fopen(name.c_str(),"a+");
    if(file != NULL){
        fprintf(file,"============================================================\n");
+       fprintf(file,"NX = %d - rank = %d - threads = %d\n",topo_in_->nglob(0),comm_size,omp_get_max_threads());
        if(is_all2all_) fprintf(file,"- is all to all\n");
        if(!is_all2all_) fprintf(file,"- is all to all VECTOR\n");
        
        int newrank;
        MPI_Comm_rank(subcomm_,&newrank);
        int rlen;
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(subcomm_, myname, &rlen);
        fprintf(file,"- in subcom %s with rank %d/%d\n",myname,newrank,subsize);
        fprintf(file,"- nglob = %d %d %d to %d %d %d\n",topo_in_->nglob(0),topo_in_->nglob(1),topo_in_->nglob(2),topo_out_->nglob(0),topo_out_->nglob(1),topo_out_->nglob(2));
        fprintf(file,"- nproc = %d %d %d to %d %d %d\n",topo_in_->nproc(0),topo_in_->nproc(1),topo_in_->nproc(2),topo_out_->nproc(0),topo_out_->nproc(1),topo_out_->nproc(2));
        // fprintf(file,"- start = %d %d %d to %d %d %d\n",istart_[0],istart_[1],istart_[2],ostart_[0],ostart_[1],ostart_[2]);
        // fprintf(file,"- end = %d %d %d to %d %d %d\n",iend_[0],iend_[1],iend_[2],oend_[0],oend_[1],oend_[2]);
        // int totalsize = (nByBlock_[0]+exSize_[0]%2)*(nByBlock_[1]+exSize_[1]%2)*(nByBlock_[2]+exSize_[2]%2);
        // fprintf(file,"- nByBlock = %d %d %d, real size = %d %d %d, alignement padding? %d vs %d\n",nByBlock_[0],nByBlock_[1],nByBlock_[2],(nByBlock_[0] == 1)?1:nByBlock_[0]+exSize_[0]%2,(nByBlock_[1] == 1)?1:nByBlock_[1]+exSize_[1]%2,(nByBlock_[2] == 1)?1:nByBlock_[2]+exSize_[2]%2,totalsize,get_blockMemSize());

        fprintf(file,"--------------------------\n");
        fprintf(file,"%d inblock: %d\n",newrank,inBlock_);
        fprintf(file,"%d I2O:\n",newrank);
        for(int ib=0; ib<inBlock_; ib++){
            fprintf(file," block[%d]: size= %d %d %d - istart = %d %d %d - destrank = %d \n",ib,iBlockSize_[0][ib],iBlockSize_[1][ib],iBlockSize_[2][ib],iBlockiStart_[0][ib],iBlockiStart_[1][ib],iBlockiStart_[2][ib],i2o_destRank_[ib]);
        }
        fprintf(file,"\n");
        fprintf(file,"--------------------------\n");
        fprintf(file,"%d onblock: %d \n",newrank,onBlock_);
        fprintf(file,"%d O2I:\n",newrank);
        for(int ib=0; ib<onBlock_; ib++){
            fprintf(file," block[%d]: size= %d %d %d - istart = %d %d %d - destrank = %d \n",ib,oBlockSize_[0][ib],oBlockSize_[1][ib],oBlockSize_[2][ib],oBlockiStart_[0][ib],oBlockiStart_[1][ib],oBlockiStart_[2][ib],o2i_destRank_[ib]);
        }
        fprintf(file,"\n");
        fclose(file);
    }
#endif
    END_FUNC;
}

/**
 * @brief Destroy the Switch Topo
 * 
 */
SwitchTopo_a2a::~SwitchTopo_a2a() {
    BEGIN_FUNC;

    int  rlen, comp;
    // if the subcom is not equal to the income we need to free it
    MPI_Comm_compare(subcomm_,inComm_,&comp);
    if(comp!=MPI_IDENT){
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(subcomm_, myname, &rlen);
        FLUPS_INFO("freeing the comm %s",myname);
        MPI_Comm_free(&subcomm_);
    }

    FLUPS_INFO("freeing the arrays");

    if (i2o_count_ != NULL) flups_free(i2o_count_);
    if (o2i_count_ != NULL) flups_free(o2i_count_);
    if (i2o_start_ != NULL) flups_free(i2o_start_);
    if (o2i_start_ != NULL) flups_free(o2i_start_);

    if (sendBuf_ != NULL) flups_free((double*)sendBuf_);
    if (recvBuf_ != NULL) flups_free((double*)recvBuf_);

    free_blockInfo_();

    if (i2o_shuffle_ != NULL) {
        for (int ib = 0; ib < onBlock_; ib++) {
            fftw_destroy_plan(i2o_shuffle_[ib]);
        }
        flups_free(i2o_shuffle_);
    }
    if (o2i_shuffle_ != NULL) {
        for (int ib = 0; ib < inBlock_; ib++) {
            fftw_destroy_plan(o2i_shuffle_[ib]);
        }
        flups_free(o2i_shuffle_);
    }

    END_FUNC;
}


/**
 * @brief Setup the communication buffer for every block
 * 
 * We associate #sendBuf_[i_block] with the correct adress inside the raw buffer.
 * This way, we can use only #sendBuf_ and #recvBuf_ for each block without any additional computation inside the execute.
 * Moreover, asking the user to allocate the data reduces the memory footprint as it can be shared among several SwitchTopo
 * 
 * @param sendData the "raw" communication buffer allocated at least at the size returned by get_bufMemSize 
 * @param recvData the "raw" communication buffer allocated at least at the size returned by get_bufMemSize 
 */
void SwitchTopo_a2a::setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) {
    BEGIN_FUNC;

    // determine the nf: since topo_in may have change, we take the max to have the correct one
    const int nf = std::max(topo_in_->nf(),topo_out_->nf());
    const int lda = std::max(topo_in_->lda(),topo_out_->lda());
    
    // store the buffers
    sendBufG_ = sendData;
    recvBufG_ = recvData;

    // conpute the subsizes
    int subsize;  //I can use subsize because, if there is no subcomm, this is mastercomm
    MPI_Comm_size(subcomm_, &subsize);

    // allocate the second layer of buffers
    sendBuf_ = (double**)flups_malloc(inBlock_ * sizeof(double*));
    recvBuf_ = (double**)flups_malloc(onBlock_ * sizeof(double*));

    // determine if we have to suffle or not
    const bool doShuffle=(topo_in_->axis() != topo_out_->axis());
    // allocate the plan if we have to suffle
    if (doShuffle) {
        i2o_shuffle_ = (fftw_plan*)flups_malloc(onBlock_ * sizeof(fftw_plan));
        o2i_shuffle_ = (fftw_plan*)flups_malloc(inBlock_ * sizeof(fftw_plan));
    } else {
        i2o_shuffle_ = NULL;
        o2i_shuffle_ = NULL;
    }
    
    // link the buff of every block to the data initialized
    int* countPerRank = (int*)flups_malloc(subsize * sizeof(int));
    std::memset(countPerRank, 0, subsize * sizeof(int));
    for (int ib = 0; ib < inBlock_; ib++) {
        // get the destination rank
        int destrank = i2o_destRank_[ib];
        // we count the number of unkowns in that will be in the buffer in front of me
        size_t memblocks = 0;
        for (int ir = 0; ir < destrank; ir++) {
            memblocks += (size_t)i2o_count_[ir]; //already accounts for lda
        }
        // place the block given the number of ranks bellow + the number of block already set to my rank
        sendBuf_[ib] = sendData + memblocks + countPerRank[destrank];
        // add the block size to the number of already added blocks
        countPerRank[destrank] += get_blockMemSize(ib, nf, iBlockSize_)*lda;

        // setup the suffle plan for the out 2 in transformation if needed
        if (doShuffle) {
            int tmp_size[3] = {iBlockSize_[0][ib], iBlockSize_[1][ib], iBlockSize_[2][ib]};
            setup_shuffle_(tmp_size, topo_out_, topo_in_, sendBuf_[ib], &o2i_shuffle_[ib]);
        }
    }

    std::memset(countPerRank, 0, subsize * sizeof(int));
    for (int ib = 0; ib < onBlock_; ib++) {
        // get the destination rank
        int destrank = o2i_destRank_[ib];
        // we count the number of unkowns in that will be in the buffer in front of me
        size_t memblocks = 0;
        for (int ir = 0; ir < destrank; ir++) {
            memblocks += (size_t)o2i_count_[ir];
        }
        // place the block given the number of ranks bellow + the number of block already set to my rank
        recvBuf_[ib] = recvData + memblocks + countPerRank[destrank];
        // add the block size to the number of already added blocks
        countPerRank[destrank] += get_blockMemSize(ib, nf, oBlockSize_)*lda;

        // setup the suffle plan for the in 2 out transformation
        if (doShuffle) {
            int tmp_size[3] = {oBlockSize_[0][ib], oBlockSize_[1][ib], oBlockSize_[2][ib]};
            setup_shuffle_(tmp_size,topo_in_,topo_out_, recvBuf_[ib], &i2o_shuffle_[ib]);
        }
    }

    flups_free(countPerRank);
    END_FUNC;
}

/**
 * @brief execute the switch from one topo to another
 * 
 * #### Buffer writting
 * The buffer memory writting is done according to the axis of the input topologies.
 * This allows to have a continuous memory access while filling the buffer.
 * 
 * We go through each block and we fill it using the local memory.
 * After a block has been filled it is send using the non-blocking send.
 * Since the writting of buffers is aligned with the topo_in axis, the loops are continuous in memory and fully vectorizable.
 * 
 * @warning
 * Let us note that the block is send with a tag which is its local index in the destination proc.
 * 
 * #### Buffer reading
 * We wait to receive one block of memory. Once one has been received, we do the copy.
 * The buffer reading has to follow the same order as in the buffer writting, so the axis of the topo_in in the inner loop.
 * 
 * The reading of the buffer is hence continuous but the writting inside the memory has an apriori unkown stride.
 * The stride may be computed using the difference of axis between the two topologies.
 * Hence the reading will be a bit slower since the writting due to memory discontinuities
 * 
 * @param v the memory to switch from one topo to another. It has to be large enough to contain both local data's
 * @param sign if the switch is forward (FLUPS_FORWARD) or backward (FLUPS_BACKWARD) w.r.t. the order defined at init.
 * 
 * -----------------------------------------------
 * We do the following:
 */
void SwitchTopo_a2a::execute(double* v, const int sign) const {
    BEGIN_FUNC;

    FLUPS_CHECK(topo_in_->isComplex() == topo_out_->isComplex(), "both topologies have to be complex or real", LOCATION);
    FLUPS_CHECK(topo_in_->lda() == topo_out_->lda(), "both topologies must have the same lda", LOCATION);
    FLUPS_CHECK(topo_in_->nf() <= 2, "the value of nf is not supported", LOCATION);
    FLUPS_CHECK(sendBuf_!=NULL && recvBuf_ != NULL, "both buffers have to be non NULL",LOCATION);

    int comm_size;
    // MPI_Comm_rank(subcomm_, &rank);
    MPI_Comm_size(subcomm_, &comm_size);

    PROF_START("reorder");

    //-------------------------------------------------------------------------
    /** - setup required memory arrays */
    //-------------------------------------------------------------------------

    const Topology* topo_in;
    const Topology* topo_out;

    int lda = topo_in_->lda();
    int send_nBlock;
    int recv_nBlock;

    // int istart[3];
    // int ostart[3];
    // int iend[3];
    // int oend[3];
    int inmem[3];
    int onmem[3];

    int* iBlockSize[3];
    int* oBlockSize[3];
    int* iBlockiStart[3];
    int* oBlockiStart[3];

    int* send_count;
    int* recv_count;
    int* send_start;
    int* recv_start;

    // const int nByBlock[3] = {nByBlock_[0], nByBlock_[1], nByBlock_[2]};

    fftw_plan* shuffle = NULL;

    opt_double_ptr* sendBuf;
    opt_double_ptr* recvBuf;
    opt_double_ptr sendBufG;
    opt_double_ptr recvBufG;

    if (sign == FLUPS_FORWARD) {
        topo_in  = topo_in_;
        topo_out = topo_out_;
        sendBuf  = sendBuf_;
        recvBuf  = recvBuf_;
        sendBufG  = sendBufG_;
        recvBufG  = recvBufG_;

        send_count = i2o_count_;
        recv_count = o2i_count_;
        send_start = i2o_start_;
        recv_start = o2i_start_;

        send_nBlock = inBlock_;
        recv_nBlock = onBlock_;

        shuffle = i2o_shuffle_;

        for (int id = 0; id < 3; id++) {
            // send_nBlock[id] = inBlock_[id];
            // recv_nBlock[id] = onBlock_[id];
            // istart[id]      = istart_[id];
            // iend[id]        = iend_[id];
            // ostart[id]      = ostart_[id];
            // oend[id]        = oend_[id];
            inmem[id]       = topo_in_->nmem(id);
            onmem[id]       = topo_out_->nmem(id);
            iBlockSize[id]  = iBlockSize_[id];
            oBlockSize[id]  = oBlockSize_[id];
            iBlockiStart[id]  = iBlockiStart_[id];
            oBlockiStart[id]  = oBlockiStart_[id];
        }
    } else if (sign == FLUPS_BACKWARD) {
        topo_in  = topo_out_;
        topo_out = topo_in_;
        sendBuf  = recvBuf_;
        recvBuf  = sendBuf_;
        sendBufG  = recvBufG_;
        recvBufG  = sendBufG_;

        send_count = o2i_count_;
        recv_count = i2o_count_;
        send_start = o2i_start_;
        recv_start = i2o_start_;

        shuffle     = o2i_shuffle_;
        send_nBlock = onBlock_;
        recv_nBlock = inBlock_;

        for (int id = 0; id < 3; id++) {
            // istart[id]      = ostart_[id];
            // iend[id]        = oend_[id];
            // ostart[id]      = istart_[id];
            // oend[id]        = iend_[id];
            inmem[id]         = topo_out_->nmem(id);
            onmem[id]         = topo_in_->nmem(id);
            iBlockSize[id]    = oBlockSize_[id];
            oBlockSize[id]    = iBlockSize_[id];
            iBlockiStart[id]  = oBlockiStart_[id];
            oBlockiStart[id]  = iBlockiStart_[id];
        }
    } else {
        FLUPS_CHECK(false, "the sign is not FLUPS_FORWARD nor FLUPS_BACKWARD", LOCATION);
    }

    FLUPS_INFO("switch: previous topo: %d,%d,%d axis=%d", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    FLUPS_INFO("switch: new topo: %d,%d,%d  axis=%d", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());
    FLUPS_INFO("switch: using %d blocks on send and %d on recv", send_nBlock, recv_nBlock);

    // check if we can return already, because the switchtopo would be useless
    {
        int rank;
        MPI_Comm_rank(subcomm_, &rank);

        bool cond = (topo_in_->axis() == topo_out_->axis()); //same axis
        cond &= (send_nBlock == 1); //only one block on this proc
        cond &= (recv_nBlock == 1);   
        cond &= (i2o_destRank_[0] == rank) ; //the only block will stay with me
        cond &= (o2i_destRank_[0] == rank) ;
        for (int i = 0; i < 3; i++) {
            cond &= (shift_[i] == 0); //no shift in memory
            cond &= (topo_in_->nloc(i) == topo_out_->nloc(i)); //same size of topology
        }
        cond &= (inmem[topo_in->axis()] == onmem[topo_out->axis()]); //same size in memory in the FRI (also for alignement)
        if(cond){
            FLUPS_INFO("I skip this switch because nothing needs to change.");
            PROF_STOP("reorder");
            return void();
        }
    };

    // define important constants
    const int iax0 = topo_in->axis();
    const int iax1 = (iax0 + 1) % 3;
    const int iax2 = (iax0 + 2) % 3;
    const int oax0 = topo_out->axis();
    const int oax1 = (oax0 + 1) % 3;
    const int oax2 = (oax0 + 2) % 3;
    const int nf   = topo_in->nf();

    PROF_STARTi("switch",iswitch_);
    PROF_STARTi("mem2buf",iswitch_);

    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    // const int nblocks_send = send_nBlock[0] * send_nBlock[1] * send_nBlock[2];

#pragma omp parallel proc_bind(close) default(none) firstprivate(send_nBlock, v, sendBuf, iBlockSize,iBlockiStart, nf, inmem, iax0, iax1, iax2, lda)
    for (int bid = 0; bid < send_nBlock; bid++) {
        for (int lia = 0; lia<lda; lia++){
            // // get the split index
            // int ibv[3];
            // localSplit(bid, send_nBlock, 0, ibv, 1);

            // // get the starting index in the global memory using !!nByBlock!!
            // // since only the last block may have a different size
            // const int loci0 = istart[iax0] + ibv[iax0] * nByBlock[iax0];
            // const int loci1 = istart[iax1] + ibv[iax1] * nByBlock[iax1];
            // const int loci2 = istart[iax2] + ibv[iax2] * nByBlock[iax2];

            // total size of a block, 1 component
            const size_t blockSize = iBlockSize[iax0][bid] * iBlockSize[iax1][bid] * iBlockSize[iax2][bid] * nf;

            // go inside the block
            const int id_max = iBlockSize[iax1][bid] * iBlockSize[iax2][bid];
            const size_t nmax = iBlockSize[iax0][bid] * nf;

            // the buffer is aligned if the starting id is aligned and if nmax is a multiple of the alignement
            const bool isBuffAligned = FLUPS_ISALIGNED(sendBuf[bid] + lia * blockSize) && nmax % FLUPS_ALIGNMENT == 0;
            // the data is aligned if the starting index is aligned AND if the gap between two entries, inmem[iax0] is a multiple of the alignment
            FLUPS_INFO_3("block %d: Moving the pointer by %d %d %d elements", bid, iBlockiStart[0][bid], iBlockiStart[1][bid], iBlockiStart[2][bid]);
            FLUPS_INFO_3("block %d: Tackling a block of size %d %d %d", bid, iBlockSize[0][bid], iBlockSize[1][bid], iBlockSize[2][bid]);
            double* my_v = v + localIndex(iax0, iBlockiStart[iax0][bid], iBlockiStart[iax1][bid], iBlockiStart[iax2][bid], iax0, inmem, nf, lia);

            const bool isVectorAligned = FLUPS_ISALIGNED(my_v) && inmem[iax0] % FLUPS_ALIGNMENT == 0;

            // we choose the best loop depending on the alignement
            if (isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / iBlockSize[iax1][bid];
                    const int i1 = id % iBlockSize[iax1][bid];

                    // get the local starting location for the buffer and the field.
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    const opt_double_ptr vloc = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf, 0);
                    opt_double_ptr dataloc    = sendBuf[bid] + id * nmax + lia * blockSize;
                    // set the alignment
                    FLUPS_ASSUME_ALIGNED(vloc, FLUPS_ALIGNMENT);
                    FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                    // do the copy -> vectorized
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        dataloc[i0] = vloc[i0];
                    }
                }
            } else if (isBuffAligned && !isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / iBlockSize[iax1][bid];
                    const int i1 = id % iBlockSize[iax1][bid];

                    // get the local starting location for the buffer and the field
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    const double* __restrict vloc = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf, 0);
                    opt_double_ptr dataloc        = sendBuf[bid] + lia * blockSize + id * nmax ;
                    // set the alignment
                    FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                    // do the copy -> vectorized
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        dataloc[i0] = vloc[i0];
                    }
                }
            } else if (!isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / iBlockSize[iax1][bid];
                    const int i1 = id % iBlockSize[iax1][bid];

                    // get the local starting location for the buffer and the field
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    const opt_double_ptr vloc  = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf, 0);
                    double* __restrict dataloc = sendBuf[bid] + lia * blockSize + id * nmax ;
                    // set the alignment
                    FLUPS_ASSUME_ALIGNED(vloc, FLUPS_ALIGNMENT);
                    // do the copy -> vectorized
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        dataloc[i0] = vloc[i0];
                    }
                }
            }else{
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / iBlockSize[iax1][bid];
                    const int i1 = id % iBlockSize[iax1][bid];

                    // get the local starting location for the buffer and the field
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    const double* __restrict vloc = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf, 0);
                    double* __restrict dataloc = sendBuf[bid] + lia * blockSize + id * nmax ;

                    // do the copy -> vectorized
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        dataloc[i0] = vloc[i0];
                    }
                }
            }
        }
    }
    
    PROF_STOPi("mem2buf",iswitch_);

    //-------------------------------------------------------------------------
    /** - Do the communication */
    //-------------------------------------------------------------------------
    if (is_all2all_) {
        PROF_STARTi("all_2_all",iswitch_);
        MPI_Alltoall(sendBufG, send_count[0], MPI_DOUBLE, recvBufG, recv_count[0], MPI_DOUBLE, subcomm_);
#ifdef PROF        
        if (prof_ != NULL) {
            string profName = "all_2_all"+to_string(iswitch_);
            prof_->stop(profName);
            int loc_mem = send_count[0] * comm_size;
            prof_->addMem(profName, loc_mem*sizeof(double));
        }
#endif

    } else {
        PROF_STARTi("all_2_all_v",iswitch_)
        MPI_Alltoallv(sendBufG, send_count, send_start, MPI_DOUBLE, recvBufG, recv_count, recv_start, MPI_DOUBLE, subcomm_);
#ifdef PROF        
        if (prof_ != NULL) {
            string profName = "all_2_all_v"+to_string(iswitch_);
            prof_->stop(profName);
            int loc_mem = 0;
            for (int ir = 0; ir < comm_size; ir++) {
                loc_mem += send_count[ir];
            }
            prof_->addMem(profName, loc_mem*sizeof(double));
        }
#endif        
    }

    //-------------------------------------------------------------------------
    /** - reset the memory to 0 */
    //-------------------------------------------------------------------------
    // reset the memory to 0
    const size_t nmax = topo_out->memsize();
    if (FLUPS_ISALIGNED(v)) {
        opt_double_ptr my_v = v;
        // tell the compiler about alignment
        FLUPS_ASSUME_ALIGNED(my_v, FLUPS_ALIGNMENT);
#pragma omp parallel for default(none) proc_bind(close) firstprivate(my_v, nmax)
        for (size_t id = 0; id < nmax; id++) {
            my_v[id] = 0.0;
        }
    } else {
        double* __restrict my_v = v;
#pragma omp parallel for default(none) proc_bind(close) firstprivate(my_v, nmax)
        for (size_t id = 0; id < nmax; id++) {
            my_v[id] = 0.0;
        }
    }

    //-------------------------------------------------------------------------
    /** - wait for a block and copy when it arrives */
    //-------------------------------------------------------------------------
    // get some counters
    // const int nblocks_recv = recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2];
    // for each block
    PROF_STARTi("buf2mem",iswitch_);

    
#pragma omp parallel default(none) proc_bind(close) firstprivate(shuffle, recv_nBlock, v, recvBuf, oBlockSize,oBlockiStart, nf, onmem, oax0, oax1, oax2, lda)
    for (int bid = 0; bid < recv_nBlock; bid++) {
        const size_t blockSize = oBlockSize[oax0][bid] * oBlockSize[oax1][bid] * oBlockSize[oax2][bid] * nf;

        // shuffle the block to get the correct index order
#pragma omp master
        {
            // only the master calls the fftw_execute which is executed in multithreading
            if (shuffle != NULL) {

//---> should shuffle all the buffer at once                
                for (int lia = 0; lia < lda; lia++){
                    // fftw_execute(shuffle[bid]);
                    if( nf == 1){
                        fftw_execute_r2r(shuffle[bid], recvBuf[bid] + lia * blockSize , recvBuf[bid] + lia * blockSize);
                    } else {
                        fftw_execute_dft(shuffle[bid], (opt_complex_ptr) (recvBuf[bid] + lia * blockSize),(opt_complex_ptr)  (recvBuf[bid] + lia * blockSize));
                    }
                }
            }
        }
#pragma omp barrier
        for (int lia = 0; lia < lda; lia++){
            // // get the indexing of the block in 012-indexing
            // int ibv[3];
            // localSplit(bid, recv_nBlock, 0, ibv, 1);

            // // get the starting index in the global memory using !!nByBlock!!
            // // since only the last block may have a different size
            // const int loci0 = ostart[oax0] + ibv[oax0] * nByBlock[oax0];
            // const int loci1 = ostart[oax1] + ibv[oax1] * nByBlock[oax1];
            // const int loci2 = ostart[oax2] + ibv[oax2] * nByBlock[oax2];

            // go inside the block
            const int id_max = oBlockSize[oax1][bid] * oBlockSize[oax2][bid];
            const size_t nmax = oBlockSize[oax0][bid] * nf;

            // the buffer is aligned if the starting id is aligned and if nmax is a multiple of the alignement
            const bool isBuffAligned = FLUPS_ISALIGNED(recvBuf[bid] + lia * blockSize) &&  nmax%FLUPS_ALIGNMENT == 0;
            // the data is aligned if the starting index is aligned AND if the gap between two entries, inmem[iax0] is a multiple of the alignment
            double*    my_v            = v + localIndex(oax0, oBlockiStart[oax0][bid], oBlockiStart[oax1][bid], oBlockiStart[oax2][bid], oax0, onmem, nf, lia);
            const bool isVectorAligned = FLUPS_ISALIGNED(my_v) && onmem[oax0] % FLUPS_ALIGNMENT == 0;

            //choose the correct loop to improve the efficiency
            if (isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / oBlockSize[oax1][bid];
                    const int i1 = id % oBlockSize[oax1][bid];
                    // get the local starting id for the buffer and the data
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    opt_double_ptr       vloc    = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf, 0);
                    const opt_double_ptr dataloc = recvBuf[bid] + lia * blockSize + id * nmax;
                    // tell the compiler about alignment
                    FLUPS_ASSUME_ALIGNED(vloc, FLUPS_ALIGNMENT);
                    FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                    // do the copy
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        vloc[i0] = dataloc[i0];
                    }
                }
            } else if (isBuffAligned && !isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / oBlockSize[oax1][bid];
                    const int i1 = id % oBlockSize[oax1][bid];
                    // get the local starting id for the buffer and the data
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    double* __restrict vloc      = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf, 0);
                    const opt_double_ptr dataloc = recvBuf[bid] + lia * blockSize + id * nmax;
                    // tell the compiler about alignment
                    FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                    // do the copy
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        vloc[i0] = dataloc[i0];
                    }
                }
            } else if (!isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / oBlockSize[oax1][bid];
                    const int i1 = id % oBlockSize[oax1][bid];
                    // get the local starting id for the buffer and the data
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    opt_double_ptr vloc              = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf, 0);
                    const double* __restrict dataloc = recvBuf[bid] + lia * blockSize + id * nmax;
                    // tell the compiler about alignment
                    FLUPS_ASSUME_ALIGNED(vloc, FLUPS_ALIGNMENT);
                    // do the copy
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        vloc[i0] = dataloc[i0];
                    }
                }
            } else {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / oBlockSize[oax1][bid];
                    const int i1 = id % oBlockSize[oax1][bid];
                    // get the local starting id for the buffer and the data
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    double* __restrict vloc          = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf, 0);
                    const double* __restrict dataloc = recvBuf[bid] + lia * blockSize + id * nmax;
                    // do the copy
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        vloc[i0] = dataloc[i0];
                    }
                }
            }
        }
    }

    PROF_STOPi("buf2mem",iswitch_);
    PROF_STOPi("switch",iswitch_);
    PROF_STOP("reorder");
    END_FUNC;
}

void SwitchTopo_a2a::disp() const {
    BEGIN_FUNC;
    FLUPS_INFO("------------------------------------------");
    if (is_all2all_) FLUPS_INFO("## Topo Swticher All to All !! MPI");
    if (!is_all2all_) FLUPS_INFO("## Topo Swticher All to All vector MPI");
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
    FLUPS_INFO("--- BLOCKS");
    // FLUPS_INFO("  - nByBlock  = %d %d %d", nByBlock_[0], nByBlock_[1], nByBlock_[2]);
    FLUPS_INFO("  - inBlock = %d", inBlock_);
    FLUPS_INFO("  - onBlock = %d", onBlock_);
    FLUPS_INFO("------------------------------------------");
}

/**
 * @brief 
 * 
 * Test of the switchtopo with fieldstart non zero, and thus exclusion of data after switch
 */
void SwitchTopo_a2a_test() {
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {18, 17, 18};
    // const int nglob[3] = {2, 2, 2};
    const int nproc[3] = {1, 1, 1};

    const int nglob_big[3] = {17, 17, 18};
    // const int nglob_big[3] = {2, 2, 2};
    const int nproc_big[3] = {1, 1, 1};
    {
        //===========================================================================
        // real numbers
        Topology* topo    = new Topology(0, 1, nglob, nproc, false, NULL, 1, MPI_COMM_WORLD);
        Topology* topobig = new Topology(1, 1, nglob_big, nproc_big, false, NULL, 1, MPI_COMM_WORLD);

        topo->disp();
        topobig->disp();

        double* data = (double*)flups_malloc(sizeof(double) * std::max(topo->memsize(), topobig->memsize()));

        const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
        for (int i2 = 0; i2 < topo->nloc(2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                    const size_t id = localIndex(0, i0, i1, i2, 0, nmem, 1, 0);

                    data[id] = (double)id;
                }
            }
        }
        // try the dump
        hdf5_dump(topo, "test_real", data);

        const int fieldstart[3] = {-1, 0, 0};
        // printf("\n=============================");
        SwitchTopo*    switchtopo = new SwitchTopo_a2a(topo, topobig, fieldstart, NULL);
        switchtopo->setup();
        size_t         max_mem    = switchtopo->get_bufMemSize();
        opt_double_ptr send_buff  = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
        opt_double_ptr recv_buff  = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
        std::memset(send_buff, 0, max_mem * sizeof(double));
        std::memset(recv_buff, 0, max_mem * sizeof(double));
        // associate the buffer
        switchtopo->setup_buffers(send_buff, recv_buff);
        switchtopo->disp();

        // printf("\n\n============ FORWARD =================");
        switchtopo->execute(data, FLUPS_FORWARD);
        hdf5_dump(topobig, "test_real_padd", data);

        // printf("\n\n============ BACKWARD =================");
        switchtopo->execute(data, FLUPS_BACKWARD);

        hdf5_dump(topo, "test_real_returned", data);
        MPI_Barrier(MPI_COMM_WORLD);

        flups_free(data);
        flups_free((double*)send_buff);
        flups_free((double*)recv_buff);
        delete (switchtopo);
        delete (topo);
        delete (topobig);
    }

    // //===========================================================================
    // complex numbers
    {
        Topology* topo    = new Topology(0, 1, nglob, nproc, true, NULL, 1, MPI_COMM_WORLD);
        Topology* topobig = new Topology(1, 1, nglob_big, nproc_big, true, NULL, 1, MPI_COMM_WORLD);

        double* data = (double*)flups_malloc(sizeof(double) * std::max(topo->memsize(), topobig->memsize()));

        const int nmem2[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
        for (int i2 = 0; i2 < topo->nloc(2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                    size_t id    = localIndex(0, i0, i1, i2, 0, nmem2, 2, 0);
                    data[id + 0] = 0;
                    data[id + 1] = id;
                }
            }
        }
        // try the dump
        hdf5_dump(topo, "test_complex", data);

        // topobig->switch2complex();
        // printf("as complex: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));
        // topobig->switch2real();
        // printf("as real: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));

        const int fieldstart2[3] = {-1, 0, 0};
        // printf("\n=============================");
        SwitchTopo* switchtopo               = new SwitchTopo_a2a(topo, topobig, fieldstart2, NULL);
        switchtopo->setup();
        switchtopo->disp();
        size_t         max_mem   = switchtopo->get_bufMemSize();
        opt_double_ptr send_buff = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
        opt_double_ptr recv_buff = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
        std::memset(send_buff, 0, max_mem * sizeof(double));
        std::memset(recv_buff, 0, max_mem * sizeof(double));
        // associate the buffer
        switchtopo->setup_buffers(send_buff, recv_buff);

        switchtopo->execute(data, FLUPS_FORWARD);

        hdf5_dump(topobig, "test_complex_padd", data);

        switchtopo->execute(data, FLUPS_BACKWARD);

        hdf5_dump(topo, "test_complex_returned", data);

        flups_free(data);
        delete (switchtopo);
        delete (topo);
        delete (topobig);
    }
    END_FUNC;
}

/**
 * @brief 
 * Test of different ranksplit, and switch to a topo with a graph_comm
 */
void SwitchTopo_a2a_test2() {
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // const int nglob[3] = {24, 24, 24};
    // const int nproc[3] = {1, 3, 2};

    // const int nglob_big[3] = {24, 24, 24};
    // const int nproc_big[3] = {1, 2, 3};
    // const int axproc[3] = {1,0,2};

    const int nglob[3] = {20, 20, 20};
    const int nproc[3] = {1, 3, 1};

    const int nglob_big[3] = {20, 20, 20};
    const int nproc_big[3] = {1, 1, 3};
    // const int nproc_big[3] = {1, 3, 2};
    const int axproc[3] = {0,1,2};

    {
        //===========================================================================
        // real numbers
        Topology* topo    = new Topology(0, 1, nglob, nproc, false, NULL, 1, MPI_COMM_WORLD);
        Topology* topobig = new Topology(1, 1, nglob_big, nproc_big, false, axproc, 1, MPI_COMM_WORLD);

        topo->disp();
        topobig->disp();

        const int fieldstart[3] = {0, 0, 0};

        //CREATE THE SWITCHTOPO BEFORE CHANGING THE TOPOS
        // printf("\n=============================");
        SwitchTopo*    switchtopo = new SwitchTopo_a2a(topo, topobig, fieldstart, NULL);
        
        MPI_Comm graph_comm = NULL;

#ifndef DEV_SIMULATE_GRAPHCOMM
        const int per[3] = {0,0,0};
        const int dims[3] = {topo->nproc(0),topo->nproc(1),topo->nproc(2)};
        MPI_Cart_create(MPI_COMM_WORLD, 3, dims,  per,  0,  &graph_comm);
#else
        int s;
        //simulate a new comm with reordered ranks:
        // int       outRanks[6] = {0, 3, 4, 1, 2, 5}; s=6;
        // int       outRanks[6] = {0, 1, 4, 2, 3, 5}; s=6;
            //CAUTION: rank i goes in posisition outRanks[i] in the new comm
            //the associated rank will be {0 1 3 4 2 5}

        int       outRanks[6] = {0, 2, 1}; s=3;
        MPI_Group group_in, group_out;
        MPI_Comm_group(MPI_COMM_WORLD, &group_in);                //get the group of the current comm
        MPI_Group_incl(group_in, s, outRanks, &group_out);        //manually reorder the ranks
        MPI_Comm_create(MPI_COMM_WORLD, group_out, &graph_comm);  // create the new comm
        
        MPI_Comm graph_comm2 = NULL;
        MPI_Comm_create(MPI_COMM_WORLD, group_out, &graph_comm2);  // create the new comm
#endif
        std::string commname = "graph_comm";
        MPI_Comm_set_name(graph_comm, commname.c_str());

        topobig->change_comm(graph_comm);
        // topo->change_comm(graph_comm2);
        topo->disp_rank();
        topobig->disp_rank();

        double* data = (double*)flups_malloc(sizeof(double) * std::max(topo->memsize(), topobig->memsize()));       

        //Filling data (AFTER having assigned topo to a new topo)
        const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
        int istart[3];
        topo->get_istart_glob(istart);
        for (int i2 = 0; i2 < topo->nloc(2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                    const size_t id = localIndex(0, i0, i1, i2, 0, nmem, 1, 0);

                    // data[id] = (double)(i0+istart[0] + i1+istart[1] + i2+istart[2]);
                    // data[id] = (double)((i0+istart[0])/4 + (i1+istart[1])/4 + (i2+istart[2])/4);
                    data[id] = (double)( (i1+istart[1])/4 + 6*((i2+istart[2])/4));
                }
            }
        }
        // try the dump
        hdf5_dump(topo, "test_real", data);


        // //CREATE THE SWITCHTOPO AFTER CHANGE IN TOPOS
        // // printf("\n=============================");
        // SwitchTopo*    switchtopo = new SwitchTopo_a2a(topo, topobig, fieldstart, NULL);
        

        // associate the buffer
        switchtopo->setup();
        size_t         max_mem    = switchtopo->get_bufMemSize();
        opt_double_ptr send_buff  = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
        opt_double_ptr recv_buff  = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
        std::memset(send_buff, 0, max_mem * sizeof(double));
        std::memset(recv_buff, 0, max_mem * sizeof(double));
        switchtopo->setup_buffers(send_buff, recv_buff);
        switchtopo->disp();

        MPI_Barrier(MPI_COMM_WORLD);

        // printf("\n\n============ FORWARD =================");
        switchtopo->execute(data, FLUPS_FORWARD);
        hdf5_dump(topobig, "test_FWD", data);

        // printf("\n\n============ BACKWARD =================");
        switchtopo->execute(data, FLUPS_BACKWARD);

        hdf5_dump(topo, "test_BCKWD", data);

        MPI_Barrier(MPI_COMM_WORLD);

        flups_free(data);
        flups_free((double*)send_buff);
        flups_free((double*)recv_buff);
        delete (switchtopo);
        delete (topo);
        delete (topobig);
    }
    END_FUNC;
}
