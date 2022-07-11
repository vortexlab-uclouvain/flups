/**
 * @file SwitchTopo_nb.cpp
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


#include "SwitchTopo_nb.hpp"

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

SwitchTopo_nb::SwitchTopo_nb(const Topology* topo_input, const Topology* topo_output, const int shift[3], H3LPR::Profiler* prof) {
    BEGIN_FUNC;

    FLUPS_CHECK(topo_input->isComplex() == topo_output->isComplex(), "both topologies have to be the same kind");
    FLUPS_CHECK(topo_input->lda() == topo_output->lda(),"Both lda's topologies must match: in=%d vs out=%d",topo_input->lda(), topo_output->lda());

    topo_in_  = topo_input;
    topo_out_ = topo_output;

    inComm_ = topo_in_->get_comm();
    outComm_ = topo_out_->get_comm();

    prof_     = prof;
    iswitch_ = topo_out_->axproc(topo_out_->axis());
    //-------------------------------------------------------------------------
    /** - compute block info */
    //-------------------------------------------------------------------------
    // get the blockshift
    shift_[0] = shift[0];
    shift_[1] = shift[1];
    shift_[2] = shift[2];

    init_blockInfo_(topo_in_, topo_out_);

    //-------------------------------------------------------------------------
    /** - initialize the profiler    */
    //-------------------------------------------------------------------------
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
 * Finally, we compute the destination tag of each block. It is defined as the local block id of the block in the received topology.
 * 
 */
void SwitchTopo_nb::init_blockInfo_(const Topology* topo_in, const Topology* topo_out){
    BEGIN_FUNC;
    
    int comm_size,ocomm_size;
    MPI_Comm_size(inComm_, &comm_size);
    MPI_Comm_size(outComm_, &ocomm_size);

    FLUPS_CHECK(ocomm_size==comm_size,"In and out communicators must have the same size.");

    //-------------------------------------------------------------------------
    /** - get the number of blocks and for each block get the size and the destination rank */
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
    //Compute start and end
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
    i2o_destRank_ = (int*)m_calloc(inBlockv[0] * inBlockv[1] * inBlockv[2] * sizeof(int));
    o2i_destRank_ = (int*)m_calloc(onBlockv[0] * onBlockv[1] * onBlockv[2] * sizeof(int));

    // get the ranks
    cmpt_blockDestRank_(inBlockv,nByBlock,shift_,istart,topo_in,topo_out,i2o_destRank_);
    cmpt_blockDestRank_(onBlockv,nByBlock,mshift,ostart,topo_out,topo_in,o2i_destRank_);

    // try to gather blocks together if possible, rewrittes the sizes, the blockistart, the number of blocks, the ranks and the tags
    gather_blocks_(topo_in, nByBlock, istart, iend, inBlockv, iBlockSize_, iBlockiStart_, &inBlock_, &i2o_destRank_);
    gather_blocks_(topo_out, nByBlock, ostart, oend, onBlockv, oBlockSize_, oBlockiStart_, &onBlock_, &o2i_destRank_);
    // get the tags for the gathered blocks
    gather_tags_(inComm_, inBlock_, onBlock_, i2o_destRank_, o2i_destRank_, &i2o_destTag_, &o2i_destTag_);

    // allocate the requests
    i2o_sendRequest_ = (MPI_Request*)m_calloc(inBlock_ * sizeof(MPI_Request));
    i2o_recvRequest_ = (MPI_Request*)m_calloc(onBlock_ * sizeof(MPI_Request));
    o2i_sendRequest_ = (MPI_Request*)m_calloc(onBlock_ * sizeof(MPI_Request));
    o2i_recvRequest_ = (MPI_Request*)m_calloc(inBlock_ * sizeof(MPI_Request));

    END_FUNC;
}

void SwitchTopo_nb::free_blockInfo_(){

    if (i2o_destRank_ != NULL) m_free(i2o_destRank_);
    if (o2i_destRank_ != NULL) m_free(o2i_destRank_);
    if (i2o_destTag_ != NULL) m_free(i2o_destTag_);
    if (o2i_destTag_ != NULL) m_free(o2i_destTag_);

    i2o_destRank_ = NULL;
    o2i_destRank_ = NULL;
    i2o_destTag_  = NULL;
    o2i_destTag_  = NULL;

    for (int id = 0; id < 3; id++) {
        if (iBlockSize_[id] != NULL) m_free(iBlockSize_[id]);
        if (oBlockSize_[id] != NULL) m_free(oBlockSize_[id]);
        if (iBlockiStart_[id] != NULL) m_free(iBlockiStart_[id]);
        if (oBlockiStart_[id] != NULL) m_free(oBlockiStart_[id]);
        iBlockSize_[id]   = NULL;
        oBlockSize_[id]   = NULL;
        iBlockiStart_[id] = NULL;
        oBlockiStart_[id] = NULL;
    }

    if (i2o_sendRequest_ != NULL) m_free(i2o_sendRequest_);
    if (i2o_recvRequest_ != NULL) m_free(i2o_recvRequest_);
    if (o2i_sendRequest_ != NULL) m_free(o2i_sendRequest_);
    if (o2i_recvRequest_ != NULL) m_free(o2i_recvRequest_);

    i2o_sendRequest_ = NULL;
    i2o_recvRequest_ = NULL;
    o2i_sendRequest_ = NULL;
    o2i_recvRequest_ = NULL;

}

/**
 * @brief 
 * 
 */
void SwitchTopo_nb::setup(){
    BEGIN_FUNC;

    int rank, comm_size;
    MPI_Comm inComm = topo_in_->get_comm();
    MPI_Comm outComm = topo_out_->get_comm();

    MPI_Comm_rank(inComm, &rank);
    MPI_Comm_size(inComm, &comm_size);

    //Ensure that comms have not changed since init. Otherwise recompute the source/destination of blocks.
    int compIn, compOut;
    MPI_Comm_compare(inComm, inComm_, &compIn);
    MPI_Comm_compare(outComm, outComm_, &compOut);
    //if the graph communicator has the same numbering as the old commn we will skip the following
    if( compIn != MPI_CONGRUENT || compOut != MPI_CONGRUENT){
        if (rank == 0){
            FLUPS_WARNING("The inComm and/or outComm have changed since this switchtopo was created. I will recompute the communication scheme.");
        }

        inComm_ = inComm;
        outComm_ = outComm;

        //reinit the block information
        free_blockInfo_();

        //The input topo may have been reset to real, even if this switchtopo is a complex2complex. 
        //We create a tmp input topo which is complex if needed, for the computation of start and end.
        int tmp_nglob[3], tmp_nproc[3], tmp_axproc[3];
        for(int i = 0; i<3;i++){
            tmp_nglob[i] = topo_in_->nglob(i);
            tmp_nproc[i] = topo_in_->nproc(i);
            tmp_axproc[i] = topo_in_->axproc(i);
        }
        Topology* topo_in_tmp = new Topology(topo_in_->axis(),topo_in_->lda(),tmp_nglob,tmp_nproc,topo_in_->isComplex(),tmp_axproc,FLUPS_ALIGNMENT,topo_in_->get_comm());
        if(topo_out_->isComplex() && !topo_in_->isComplex()){
            topo_in_tmp->switch2complex();
        }
        //recompute block info
        init_blockInfo_(topo_in_tmp, topo_out_);

        delete(topo_in_tmp);
    }

    //-------------------------------------------------------------------------
    /** - Setup subcomm */
    //-------------------------------------------------------------------------
    cmpt_commSplit_();
    // setup the dest rank, counts and starts
    setup_subComm_(inBlock_,topo_in_->lda() ,iBlockSize_, i2o_destRank_,NULL,NULL);
    setup_subComm_(onBlock_,topo_out_->lda(),oBlockSize_, o2i_destRank_,NULL,NULL);

    //-------------------------------------------------------------------------
    /** - Compute the self blocks in the new comms   */
    //-------------------------------------------------------------------------
    int newrank;
    MPI_Comm_rank(subcomm_,&newrank);
    selfBlockN_ = 0;
    for (int bid = 0; bid < inBlock_; bid++) {
        // for the send when doing input 2 output: send to rank i2o with tag i2o_destTag_[bid]
        if (i2o_destRank_[bid] == newrank) {
            selfBlockN_++;
        }
    }
    int temp = 0;
    for (int bid = 0; bid < onBlock_; bid++) {
        if (o2i_destRank_[bid] == newrank) {
            temp++;
        }
    }
    FLUPS_CHECK(temp == selfBlockN_, "the number of selfBlocks has to be the same in both TOPO!");
    iselfBlockID_ = (int*)m_calloc(selfBlockN_ * sizeof(int));
    oselfBlockID_ = (int*)m_calloc(selfBlockN_ * sizeof(int));
    //-------------------------------------------------------------------------
    /** - Display performance information if asked */
    //-------------------------------------------------------------------------
#ifdef PERF_VERBOSE
    int rankworld;
    MPI_Comm_rank(inComm, &rankworld);
    // we display important information for the performance
    std::string name = "./prof/SwitchTopo_" + std::to_string(topo_in_->axis()) + "to" + std::to_string(topo_out_->axis()) + "rank_" + std::to_string(rankworld) + ".txt";
    FILE* file = fopen(name.c_str(),"a+");
    if(file != NULL){
        fprintf(file,"============================================================\n");
+       fprintf(file,"NX = %d - rank = %d - threads = %d\n",topo_in_->nglob(0),comm_size,omp_get_max_threads());
        fprintf(file,"- non blocking and persistent communications\n");

        // int rlen;
        // char myname[MPI_MAX_OBJECT_NAME];
        // MPI_Comm_get_name(subcomm_, myname, &rlen);
        // fprintf(file,"- in subcom %s with rank %d/%d\n",myname,newrank,subsize);
        fprintf(file,"- nglob = %d %d %d to %d %d %d\n",topo_in_->nglob(0),topo_in_->nglob(1),topo_in_->nglob(2),topo_out_->nglob(0),topo_out_->nglob(1),topo_out_->nglob(2));
        fprintf(file,"- nproc = %d %d %d to %d %d %d\n",topo_in_->nproc(0),topo_in_->nproc(1),topo_in_->nproc(2),topo_out_->nproc(0),topo_out_->nproc(1),topo_out_->nproc(2));
        // fprintf(file,"- start = %d %d %d to %d %d %d\n",istart_[0],istart_[1],istart_[2],ostart_[0],ostart_[1],ostart_[2]);
        // fprintf(file,"- end = %d %d %d to %d %d %d\n",iend_[0],iend_[1],iend_[2],oend_[0],oend_[1],oend_[2]);
        // int totalsize = (nByBlock_[0]+exSize_[0]%2)*(nByBlock_[1]+exSize_[1]%2)*(nByBlock_[2]+exSize_[2]%2)*topo_out_->nf();
        // fprintf(file,"- nByBlock = %d %d %d, real size = %d %d %d, alignement padding? %d vs %d\n",nByBlock_[0],nByBlock_[1],nByBlock_[2],(nByBlock_[0] == 1)?1:nByBlock_[0]+exSize_[0]%2,(nByBlock_[1] == 1)?1:nByBlock_[1]+exSize_[1]%2,(nByBlock_[2] == 1)?1:nByBlock_[2]+exSize_[2]%2,totalsize,get_blockMemSize());

        fprintf(file,"--------------------------\n");
        fprintf(file,"%d inblock: %d\n",newrank,inBlock_);
        fprintf(file,"%d I2O:\n",newrank);
        for(int ib=0; ib<inBlock_; ib++){
            fprintf(file," block[%d]: size= %d %d %d - istart = %d %d %d - destrank = %d - desttag = %d\n",ib,iBlockSize_[0][ib],iBlockSize_[1][ib],iBlockSize_[2][ib],iBlockiStart_[0][ib],iBlockiStart_[1][ib],iBlockiStart_[2][ib],i2o_destRank_[ib],i2o_destTag_[ib]);
        }
        fprintf(file,"\n");
        fprintf(file,"--------------------------\n");
        fprintf(file,"%d onblock: %d\n",newrank,onBlock_);
        fprintf(file,"%d O2I:\n",newrank);
        for(int ib=0; ib<onBlock_; ib++){
            fprintf(file," block[%d]: size= %d %d %d - istart = %d %d %d - destrank = %d - desttag = %d \n",ib,oBlockSize_[0][ib],oBlockSize_[1][ib],oBlockSize_[2][ib],oBlockiStart_[0][ib],oBlockiStart_[1][ib],oBlockiStart_[2][ib],o2i_destRank_[ib],o2i_destTag_[ib]);
        }
        fprintf(file,"\n");
        fclose(file);
    }
#endif

    END_FUNC;
}

void SwitchTopo_nb::setup_buffers(opt_double_ptr sendData,opt_double_ptr recvData){
    BEGIN_FUNC;


    // determine the nf: since topo_in may have change, we take the max to have the correct one
    const int nf = std::max(topo_in_->nf(),topo_out_->nf());
    const int lda = std::max(topo_in_->lda(),topo_out_->lda());

    int newrank;
    MPI_Comm_rank(subcomm_, &newrank);
    // allocate the second layer of buffers
    sendBuf_ = (double**)m_calloc(inBlock_ * sizeof(double*));
    recvBuf_ = (double**)m_calloc(onBlock_ * sizeof(double*));

    const bool doShuffle=(topo_in_->axis() != topo_out_->axis());
    
    if (doShuffle) {
        i2o_shuffle_ = (fftw_plan*)m_calloc(onBlock_ * sizeof(fftw_plan));
        o2i_shuffle_ = (fftw_plan*)m_calloc(inBlock_ * sizeof(fftw_plan));
    } else {
        i2o_shuffle_ = NULL;
        o2i_shuffle_ = NULL;
    }
    
    //-------------------------------------------------------------------------
    /** - for each block we associate the the data buffer and the MPI requests or associate it to NULL */
    //-------------------------------------------------------------------------
    // reset the counter to 0
    int selfcount = 0;
    for (int bid = 0; bid < inBlock_; bid++) {
        //associate the pointer to the correct block
        sendBuf_[bid] = sendData;
        // compute the send size
        const size_t sendSize = get_blockMemSize(bid, nf, iBlockSize_) * topo_in_->lda();
        // update the pointer
        sendData = sendData + sendSize;
        
        // if I send to myself, store the info
        if (i2o_destRank_[bid] == newrank) {
            // save the bid to the self block list
            iselfBlockID_[selfcount] = bid;
            // associate the request to NULL
            i2o_sendRequest_[bid] = MPI_REQUEST_NULL;
            o2i_recvRequest_[bid] = MPI_REQUEST_NULL;
            // increment the counter
            selfcount++;
        } else {
            // get the send size without padding
            MPI_Send_init(sendBuf_[bid], sendSize, MPI_DOUBLE, i2o_destRank_[bid], i2o_destTag_[bid], subcomm_, &(i2o_sendRequest_[bid]));
            // for the send when doing output 2 input: send to rank o2i with tag o2i
            MPI_Recv_init(sendBuf_[bid], sendSize, MPI_DOUBLE, i2o_destRank_[bid], bid, subcomm_, &(o2i_recvRequest_[bid]));
        }

        // setup the suffle plan for the out 2 in transformation if needed
        if (doShuffle) {
            int tmp_size[3] = {iBlockSize_[0][bid], iBlockSize_[1][bid], iBlockSize_[2][bid]};
            setup_shuffle_(tmp_size, topo_out_, topo_in_, sendBuf_[bid], &o2i_shuffle_[bid]);
        }
    }
    FLUPS_CHECK(selfcount == selfBlockN_, "the number of counted block has to match the allocation number: %d vs %d", selfcount, selfBlockN_);

    // reset the self count
    selfcount = 0;
    for (int bid = 0; bid < onBlock_; bid++) {
        //associate the pointer with the correct block
        recvBuf_[bid] = recvData;
        //compute the recvsize
        const size_t recvSize = get_blockMemSize(bid,nf,oBlockSize_) * topo_out_->lda();
        // update the pointer
        recvData = recvData + recvSize;
        // create the request if needed
        if (o2i_destRank_[bid] == newrank) {
            // save the bid
            oselfBlockID_[selfcount] = bid;
            // associate the request to NULL
            i2o_recvRequest_[bid] = MPI_REQUEST_NULL;
            o2i_sendRequest_[bid] = MPI_REQUEST_NULL;
            // increment the counter
            selfcount++;
        } else {
            // for the reception when doing input 2 output: receive from the rank o2i with tag bid
            MPI_Recv_init(recvBuf_[bid], recvSize, MPI_DOUBLE, o2i_destRank_[bid], bid, subcomm_, &(i2o_recvRequest_[bid]));
            // for the send when doing output 2 input: send to rank o2i with tag o2i
            MPI_Send_init(recvBuf_[bid], recvSize, MPI_DOUBLE, o2i_destRank_[bid], o2i_destTag_[bid], subcomm_, &(o2i_sendRequest_[bid]));
        }

        // setup the suffle plan for the in 2 out transformation
        if (doShuffle) {
            int tmp_size[3] = {oBlockSize_[0][bid], oBlockSize_[1][bid], oBlockSize_[2][bid]};
            setup_shuffle_(tmp_size, topo_in_, topo_out_, recvBuf_[bid], &i2o_shuffle_[bid]);
        }
    }
    FLUPS_CHECK(selfcount == selfBlockN_, "the number of counted block has to match the allocation number: %d vs %d", selfcount, selfBlockN_);
    END_FUNC;
}

/**
 * @brief Destroy the Switch Topo
 * 
 */
SwitchTopo_nb::~SwitchTopo_nb() {
    BEGIN_FUNC;

    int  rlen, comp;
    MPI_Comm_compare(subcomm_,inComm_,&comp);
    if(comp!=MPI_IDENT){
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(subcomm_, myname, &rlen);
        FLUPS_INFO("freeing the comm %s",myname);
        MPI_Comm_free(&subcomm_);
    }

    for(int ib=0; ib< inBlock_; ib++){
        // if (sendBuf_[ib] != NULL) m_free(sendBuf_[ib]);
        if (i2o_sendRequest_[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(i2o_sendRequest_[ib]));
        if (o2i_recvRequest_[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(o2i_recvRequest_[ib]));
    }
    for(int ib=0; ib< onBlock_; ib++){
        // if (recvBuf_[ib] != NULL) m_free(recvBuf_[ib]);
        if (i2o_recvRequest_[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(i2o_recvRequest_[ib]));
        if (o2i_sendRequest_[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(o2i_sendRequest_[ib]));
    }

    free_blockInfo_();

    if (iselfBlockID_ != NULL) m_free(iselfBlockID_);
    if (oselfBlockID_ != NULL) m_free(oselfBlockID_);

    if (i2o_shuffle_ != NULL) {
        for (int ib = 0; ib < onBlock_; ib++) {
            fftw_destroy_plan(i2o_shuffle_[ib]);
        }
        m_free(i2o_shuffle_);
    }
    if (o2i_shuffle_ != NULL) {
        for (int ib = 0; ib < inBlock_; ib++) {
            fftw_destroy_plan(o2i_shuffle_[ib]);
        }
        m_free(o2i_shuffle_);
    }

    m_free((double*)sendBuf_);
    m_free((double*)recvBuf_);
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
void SwitchTopo_nb::execute(double* v, const int sign) const {
    BEGIN_FUNC;

    FLUPS_CHECK(topo_in_->isComplex() == topo_out_->isComplex(),"both topologies have to be complex or real");
    FLUPS_CHECK(topo_in_->lda() == topo_out_->lda(), "both topologies must have the same lda");
    FLUPS_CHECK(topo_in_->nf() <= 2, "the value of nf is not supported");
    FLUPS_CHECK(sendBuf_!=NULL && recvBuf_ != NULL, "both buffers have to be non NULL");

    m_profStarti(prof_, "reorder");
    int iswitch = iswitch_;

    //-------------------------------------------------------------------------
    /** - setup required memory arrays */
    //-------------------------------------------------------------------------

    const Topology* topo_in;
    const Topology* topo_out;

    MPI_Request* sendRequest;
    MPI_Request* recvRequest;

    int lda = topo_in_->lda();
    int send_nBlock;
    int recv_nBlock;

    int inmem[3];
    int onmem[3];

    int* iBlockSize[3];
    int* oBlockSize[3];
    int* iBlockiStart[3];
    int* oBlockiStart[3];
    
    int* oselfBlockID;
    int* destTag;

    // const int nByBlock[3] = {nByBlock_[0],nByBlock_[1],nByBlock_[2]};

    opt_double_ptr* sendBuf;
    opt_double_ptr* recvBuf;

    fftw_plan* shuffle = NULL;

    if (sign == FLUPS_FORWARD) {
        topo_in     = topo_in_;
        topo_out    = topo_out_;
        sendRequest = i2o_sendRequest_;
        recvRequest = i2o_recvRequest_;
        sendBuf     = sendBuf_;
        recvBuf     = recvBuf_;

        oselfBlockID = oselfBlockID_;
        destTag      = i2o_destTag_;

        shuffle = i2o_shuffle_;

        send_nBlock = inBlock_;
        recv_nBlock = onBlock_;

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
        topo_in     = topo_out_;
        topo_out    = topo_in_;
        sendRequest = o2i_sendRequest_;
        recvRequest = o2i_recvRequest_;
        sendBuf     = recvBuf_;
        recvBuf     = sendBuf_;

        oselfBlockID = iselfBlockID_;
        destTag      = o2i_destTag_;

        shuffle = o2i_shuffle_;

        send_nBlock = onBlock_;
        recv_nBlock = inBlock_;

        for (int id = 0; id < 3; id++) {
            // send_nBlock[id] = onBlock_[id];
            // recv_nBlock[id] = inBlock_[id];
            // istart[id]      = ostart_[id];
            // iend[id]        = oend_[id];
            // ostart[id]      = istart_[id];
            // oend[id]        = iend_[id];
            inmem[id]       = topo_out_->nmem(id);
            onmem[id]       = topo_in_->nmem(id);
            iBlockSize[id]  = oBlockSize_[id];
            oBlockSize[id]  = iBlockSize_[id];
            iBlockiStart[id]  = oBlockiStart_[id];
            oBlockiStart[id]  = iBlockiStart_[id];
            
        }
    } else {
        FLUPS_CHECK(false, "the sign is not FLUPS_FORWARD nor FLUPS_BACKWARD");
    }

    FLUPS_INFO("switch nb: previous topo: %d,%d,%d axis=%d", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    FLUPS_INFO("switch nb: new topo: %d,%d,%d  axis=%d", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());
    FLUPS_INFO("switch nb: using %d blocks on send and %d on recv",send_nBlock,recv_nBlock);

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
            m_profStopi(prof_, "reorder");
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

    //-------------------------------------------------------------------------
    /** - start the reception requests so we are ready to receive */
    //-------------------------------------------------------------------------
    for (int bid = 0; bid < recv_nBlock; bid++) {
        if (recvRequest[bid] != MPI_REQUEST_NULL){
            MPI_Start(&(recvRequest[bid]));
        }
    }

    m_profStarti(prof_, "switch%d",iswitch);
    m_profStarti(prof_, "mem2buf%d",iswitch);
    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    // const int nblocks_send = send_nBlock[0] * send_nBlock[1] * send_nBlock[2];

#if defined(__INTEL_COMPILER)
//possible need to add ```shared(ompi_request_null)``` depending on the compiler version
#pragma omp parallel proc_bind(close) default(none) firstprivate(send_nBlock, v, sendBuf, recvBuf, destTag, iBlockSize,iBlockiStart, nf, inmem, iax0, iax1,iax2,sendRequest, lda)
#elif defined(__GNUC__)
#pragma omp parallel proc_bind(close) default(none) shared(ompi_request_null) firstprivate(send_nBlock, v, sendBuf, recvBuf, destTag,iBlockSize,iBlockiStart, nf, inmem, iax0, iax1,iax2,sendRequest, lda)
#endif
    for (int bid = 0; bid < send_nBlock; bid++) {
        for (int lia = 0; lia < lda ; lia++){
            // total size of a block, 1 component
            const size_t blockSize = iBlockSize[iax0][bid] * iBlockSize[iax1][bid] * iBlockSize[iax2][bid] * nf;

            // get the buffer data for this block
            double* data;
            if(sendRequest[bid] == MPI_REQUEST_NULL){
                // if we are doing a self block the data is the recv buff
                // the new block ID is given by destTag[bid]
                data = recvBuf[destTag[bid]] + lia * blockSize;
            } else {
                // else we copy inside the sendbuffer
                data = sendBuf[bid] + lia * blockSize;
            }

            // go inside the block
            const int id_max = iBlockSize[iax1][bid] * iBlockSize[iax2][bid];
            const size_t nmax = iBlockSize[iax0][bid] * nf;

            // the buffer is aligned if the starting id is aligned and if nmax is a multiple of the alignement
            const bool isBuffAligned = FLUPS_ISALIGNED(data) &&  nmax%FLUPS_ALIGNMENT == 0;
            // the data is aligned if the starting index is aligned AND if the gap between two entries, inmem[iax0] is a multiple of the alignment
            double*    my_v            = v + localIndex(iax0, iBlockiStart[iax0][bid], iBlockiStart[iax1][bid], iBlockiStart[iax2][bid], iax0, inmem, nf, lia);
            const bool isVectorAligned = FLUPS_ISALIGNED(my_v) && inmem[iax0] % FLUPS_ALIGNMENT == 0;

            // we choose the best loop depending on the alignement
            if (isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
                for (int id = 0; id < id_max; id++) {
                    // get the id from a small modulo
                    const int i2 = id / iBlockSize[iax1][bid];
                    const int i1 = id % iBlockSize[iax1][bid];

                    // get the local starting location for the buffer and the field
                    //   my_v has already set the address in the right portion of lda, so now,
                    //   only running over the chunks as if lda=1
                    const opt_double_ptr vloc = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf, 0);
                    opt_double_ptr dataloc    = data + id * nmax;
                    FLUPS_ASSUME_ALIGNED(vloc,FLUPS_ALIGNMENT);
                    FLUPS_ASSUME_ALIGNED(dataloc,FLUPS_ALIGNMENT);
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
                    opt_double_ptr dataloc        = data + id * nmax;
                    FLUPS_ASSUME_ALIGNED(dataloc,FLUPS_ALIGNMENT);
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
                    double* __restrict dataloc = data + id * nmax;
                    FLUPS_ASSUME_ALIGNED(vloc,FLUPS_ALIGNMENT);
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
                    const double* __restrict vloc  = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf, 0);
                    double* __restrict dataloc = data + id * nmax;

                    // do the copy -> vectorized
                    for (size_t i0 = 0; i0 < nmax; i0++) {
                        dataloc[i0] = vloc[i0];
                    }
                }
            }
        }
        // the barrier after an OpenMP "for" block is implicit
        // start the send the block and continue
#pragma omp master
        {
            if (sendRequest[bid] != MPI_REQUEST_NULL) {
                MPI_Start(&(sendRequest[bid]));
            }
        }
    }

    m_profStopi(prof_, "mem2buf%d",iswitch);

    //-------------------------------------------------------------------------
    /** - reset the memory to 0 */
    //-------------------------------------------------------------------------
    // reset the memory to 0
    const size_t nmax = topo_out->memsize();
    if (FLUPS_ISALIGNED(v)) {
        opt_double_ptr my_v = v;
        FLUPS_ASSUME_ALIGNED(my_v,FLUPS_ALIGNMENT);
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
    // const int nblocks_recv  = recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2];

    // create the status as a shared variable
    MPI_Status status;

#pragma omp parallel default(none) proc_bind(close) shared(status) firstprivate(recv_nBlock, oselfBlockID, v, recvBuf, oBlockSize, oBlockiStart, nf, onmem, oax0, oax1, oax2, recvRequest, iswitch, shuffle, lda)
    for (int count = 0; count < recv_nBlock; count++) {
        // only the master receive the call
        int bid = -1;

        // if we are doing a self block
        if (count < selfBlockN_) {
            bid = oselfBlockID[count];
            
            // total size of a block, 1 component
            const size_t blockSize = oBlockSize[oax0][bid] * oBlockSize[oax1][bid] * oBlockSize[oax2][bid] * nf;
#pragma omp master
            {
                m_profStarti(prof_, "buf2mem%d",iswitch);
                // only the master call the fftw_execute which is executed in multithreading
                if (shuffle != NULL) {
                    for (int lia = 0; lia < lda; lia++){
                        // fftw_execute(shuffle[bid]);
                        if( nf == 1){
                            fftw_execute_r2r(shuffle[bid], recvBuf[bid] + lia * blockSize, recvBuf[bid] + lia * blockSize );
                        } else {
                            fftw_execute_dft(shuffle[bid], (opt_complex_ptr) (recvBuf[bid] + lia * blockSize), (opt_complex_ptr) (recvBuf[bid] + lia * blockSize) );
                        }
                    }
                }
            }
#pragma omp barrier
        } else {
#pragma omp master
            {
                m_profStarti(prof_, "waiting%d",iswitch);
                int request_index;
                MPI_Waitany(recv_nBlock, recvRequest, &request_index, &status);
                m_profStopi(prof_, "waiting%d",iswitch);
                m_profStarti(prof_, "buf2mem%d",iswitch);
                
                // bid is set for the master
                bid = status.MPI_TAG;
                // total size of a block, 1 component
                const size_t blockSize = oBlockSize[oax0][bid] * oBlockSize[oax1][bid] * oBlockSize[oax2][bid] * nf;

                // only the master call the fftw_execute which is executed in multithreading                
                if (shuffle != NULL) {
                    for (int lia = 0; lia < lda; lia++){
                        // fftw_execute(shuffle[bid]);
                        if( nf == 1){
                            fftw_execute_r2r(shuffle[bid], recvBuf[bid] + lia * blockSize, recvBuf[bid] + lia * blockSize );
                        } else {
                            fftw_execute_dft(shuffle[bid], (opt_complex_ptr) (recvBuf[bid] + lia * blockSize), (opt_complex_ptr) (recvBuf[bid] + lia * blockSize));
                        }
                    }
                }
            }
            // make sure that the master has received the status before going further
            // there is no implicit barrier after
#pragma omp barrier
            // get the block id = the tag
            bid = status.MPI_TAG;
        }
//         // add the bandwith info
//         #pragma omp master
//         {
// #ifdef PROF            
//             if (prof_ != NULL) {
//                 prof_->addMem("waiting"+to_string(iswitch), get_blockMemSize(bid,nf,oBlockSize)*lda*sizeof(double));
//             }
// #endif
//         }

        for (int lia = 0; lia < lda; lia++){
            // total size of a block, 1 component
            const size_t blockSize = oBlockSize[oax0][bid] * oBlockSize[oax1][bid] * oBlockSize[oax2][bid] * nf;

            // // get the indexing of the block in 012-indexing
            // int ibv[3];
            // localSplit(bid, recv_nBlock, 0, ibv, 1);

            // // get the starting index in the global memory using !!nByBlock!!
            // // since only the last block may have a different size
            // const int loci0 = ostart[oax0] + ibv[oax0] * nByBlock[oax0];
            // const int loci1 = ostart[oax1] + ibv[oax1] * nByBlock[oax1];
            // const int loci2 = ostart[oax2] + ibv[oax2] * nByBlock[oax2];

            // go inside the block
            const int    id_max = oBlockSize[oax1][bid] * oBlockSize[oax2][bid];
            const size_t nmax   = oBlockSize[oax0][bid] * nf;

            // the buffer is aligned if the starting id is aligned and if nmax is a multiple of the alignement
            const bool isBuffAligned = FLUPS_ISALIGNED(recvBuf[bid] + lia * blockSize) &&  nmax%FLUPS_ALIGNMENT == 0;
            // the data is aligned if the starting index is aligned AND if the gap between two entries, inmem[iax0] is a multiple of the alignment
            // double*    my_v            = v + localIndex(oax0, oBlockiStart[0][bid], oBlockiStart[1][bid], oBlockiStart[2][bid], oax0, onmem, nf);
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
                    FLUPS_ASSUME_ALIGNED(vloc,FLUPS_ALIGNMENT);
                    FLUPS_ASSUME_ALIGNED(dataloc,FLUPS_ALIGNMENT);
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
                    FLUPS_ASSUME_ALIGNED(dataloc,FLUPS_ALIGNMENT);
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
                    FLUPS_ASSUME_ALIGNED(vloc,FLUPS_ALIGNMENT);
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

#pragma omp master
        {
            m_profStopi(prof_, "buf2mem%d",iswitch);
        }
    }
    // now that we have received everything, close the send requests
    MPI_Waitall(send_nBlock, sendRequest,MPI_STATUSES_IGNORE);

    m_profStopi(prof_, "switch%d",iswitch);
    m_profStopi(prof_, "reorder");
    END_FUNC;
}

void SwitchTopo_nb::disp() const {
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
    FLUPS_INFO("--- BLOCKS");
    FLUPS_INFO("  - selfBlockN = %d", selfBlockN_);
    // FLUPS_INFO("  - nByBlock  = %d %d %d", nByBlock_[0], nByBlock_[1], nByBlock_[2]);
    FLUPS_INFO("  - inBlock = %d", inBlock_);
    FLUPS_INFO("  - onBlock = %d", onBlock_);
    FLUPS_INFO("------------------------------------------");
}

void SwitchTopo_nb_test() {
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {8, 8, 8};
    const int nproc[3] = {2, 2, 1};

    const int nglob_big[3] = {17, 8, 8};
    const int nproc_big[3] = {2, 2, 1};

    //===========================================================================
    // real numbers
    Topology* topo    = new Topology(0, 1, nglob, nproc, false,NULL,1, MPI_COMM_WORLD);
    Topology* topobig = new Topology(0, 1, nglob_big, nproc_big, false,NULL,1, MPI_COMM_WORLD);

    double* data = (double*)m_calloc(sizeof(double*) * std::max(topo->memsize(), topobig->memsize()));

    const int nmem[3] = {topo->nmem(0),topo->nmem(1),topo->nmem(2)};
    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localIndex(0,i0, i1, i2, 0,nmem,1, 0);
                data[id] = id;
            }
        }
    }
    // try the dump
#if (FLUPS_HDF5)
    hdf5_dump(topo, "test_real", data);
#endif
    const int fieldstart[3] = {0, 0, 0};
    // printf("\n=============================");
    SwitchTopo* switchtopo = new SwitchTopo_nb(topo, topobig, fieldstart, NULL);
    switchtopo->setup();
    switchtopo->disp();

    size_t         max_mem    = switchtopo->get_bufMemSize();
    opt_double_ptr send_buff  = (opt_double_ptr)m_calloc(max_mem * sizeof(double));
    opt_double_ptr recv_buff  = (opt_double_ptr)m_calloc(max_mem * sizeof(double));
    std::memset(send_buff, 0, max_mem * sizeof(double));
    std::memset(recv_buff, 0, max_mem * sizeof(double));
    // associate the buffer
    
    switchtopo->setup_buffers(send_buff, recv_buff);

    // printf("\n\n============ FORWARD =================");
    switchtopo->execute(data, FLUPS_FORWARD);

#if (FLUPS_HDF5)
    hdf5_dump(topobig, "test_real_padd", data);
#endif

    // printf("\n\n============ BACKWARD =================");
    switchtopo->execute(data, FLUPS_BACKWARD);

#if (FLUPS_HDF5)
    hdf5_dump(topo, "test_real_returned", data);
#endif

    m_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);

    //===========================================================================
    // complex numbers
    topo    = new Topology(0, 1, nglob, nproc, true,NULL,1, MPI_COMM_WORLD);
    topobig = new Topology(2, 1, nglob_big, nproc_big, true,NULL,1, MPI_COMM_WORLD);

    data = (double*)m_calloc(sizeof(double*) * topobig->memsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localIndex(0, i0, i1, i2, 0, nmem, 2, 0);
                data[id + 0] = 0;
                data[id + 1] = id;
            }
        }
    }

#if (FLUPS_HDF5) 
    // try the dump
    hdf5_dump(topo, "test_complex", data);
#endif

    // topobig->switch2complex();
    // printf("as complex: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));
    // topobig->switch2real();
    // printf("as real: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));

    const int fieldstart2[3] = {4, 0, 0};
    // printf("\n=============================");
    switchtopo = new SwitchTopo_nb(topo, topobig, fieldstart2, NULL);

    switchtopo->execute(data, FLUPS_FORWARD);

#if (FLUPS_HDF5)
    hdf5_dump(topobig, "test_complex_padd", data);
#endif 

    switchtopo->execute(data, FLUPS_BACKWARD);

#if (FLUPS_HDF5)
    hdf5_dump(topo, "test_complex_returned", data);
#endif

    m_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);
    END_FUNC;
}