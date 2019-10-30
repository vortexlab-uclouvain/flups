/**
 * @file SwitchTopo_nb.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-28
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
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

SwitchTopo_nb::SwitchTopo_nb(const Topology* topo_input, const Topology* topo_output, const int shift[3], Profiler* prof) {
    BEGIN_FUNC;

    FLUPS_CHECK(topo_input->isComplex() == topo_output->isComplex(), "both topologies have to be the same kind", LOCATION);

    _topo_in  = topo_input;
    _topo_out = topo_output;

    _inComm = _topo_in->get_comm();
    _outComm = _topo_out->get_comm();

    int comm_size;
    MPI_Comm_size(_inComm, &comm_size);

#ifdef PROF    
    _prof     = prof;
    for (int ip=0;ip<3;ip++){
        if(_topo_in->axis() == _topo_in->axproc(ip)){
            _iswitch  = ip;
        }
    }
#endif
    //-------------------------------------------------------------------------
    /** - compute block info */
    //-------------------------------------------------------------------------
    // get the blockshift
    _shift[0] = shift[0];
    _shift[1] = shift[1];
    _shift[2] = shift[2];

    _init_blockInfo(_topo_in, _topo_out);

    //-------------------------------------------------------------------------
    /** - initialize the profiler    */
    //-------------------------------------------------------------------------
#ifdef PROF
    if (_prof != NULL) {
        _prof->create("reorder","solve");

        _prof->create("switch0", "reorder");
        _prof->create("mem2buf0", "switch0");
        _prof->create("buf2mem0", "switch0");
        _prof->create("waiting0", "switch0");

        _prof->create("switch1", "reorder");
        _prof->create("mem2buf1", "switch1");
        _prof->create("buf2mem1", "switch1");
        _prof->create("waiting1", "switch1");

        _prof->create("switch2", "reorder");
        _prof->create("mem2buf2", "switch2");
        _prof->create("buf2mem2", "switch2");
        _prof->create("waiting2", "switch2");
    }
#endif
    END_FUNC;
}

/**
 * @brief initialize the blocks: compute their index, their number, their size and their source/destination
 * 
 */
void SwitchTopo_nb::_init_blockInfo(const Topology* topo_in, const Topology* topo_out){
    BEGIN_FUNC;
    
    int comm_size,ocomm_size;
    MPI_Comm_size(_inComm, &comm_size);
    MPI_Comm_size(_outComm, &ocomm_size);

    FLUPS_CHECK(ocomm_size==comm_size,"In and out communicators must have the same size.",LOCATION);


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

    // int  iblockIDStart[3];
    // int  oblockIDStart[3];
    // int* inBlockEachProc     = (int*)flups_malloc(comm_size * 3 * sizeof(int));
    // int* onBlockEachProc     = (int*)flups_malloc(comm_size * 3 * sizeof(int));
    // int* istartBlockEachProc = (int*)flups_malloc(comm_size * 3 * sizeof(int));
    // int* ostartBlockEachProc = (int*)flups_malloc(comm_size * 3 * sizeof(int));

    //-------------------------------------------------------------------------
    /** - Compute intersection ids */
    //-------------------------------------------------------------------------
    //recompute _start and _end
    topo_in->cmpt_intersect_id(_shift, topo_out, istart, iend);
    int mshift[3] = {-_shift[0], -_shift[1], -_shift[2]};
    topo_out->cmpt_intersect_id(mshift, topo_in, ostart, oend);

    //-------------------------------------------------------------------------
    /** - get the block size as the GCD of the memory among every process between send and receive */
    //-------------------------------------------------------------------------
    _cmpt_nByBlock(istart,iend,ostart,oend,nByBlock);

    // _cmpt_blockIndexes(istart, iend, nByBlock, topo_in, inBlockv, iblockIDStart, istartBlockEachProc, inBlockEachProc);
    // _cmpt_blockIndexes(ostart, oend, nByBlock, topo_out, onBlockv, oblockIDStart, ostartBlockEachProc, onBlockEachProc);
    _cmpt_blockIndexes(istart, iend, nByBlock, topo_in, inBlockv);
    _cmpt_blockIndexes(ostart, oend, nByBlock, topo_out, onBlockv);

    // // allocte the block size
    // for (int id = 0; id < 3; id++) {
    //     _iBlockSize[id] = (int*)flups_malloc(inBlockv[0] * inBlockv[1] * inBlockv[2] * sizeof(int));
    //     _oBlockSize[id] = (int*)flups_malloc(onBlockv[0] * onBlockv[1] * onBlockv[2] * sizeof(int));
    // }

    // allocate the destination ranks
    _i2o_destRank = (int*)flups_malloc(inBlockv[0] * inBlockv[1] * inBlockv[2] * sizeof(int));
    _o2i_destRank = (int*)flups_malloc(onBlockv[0] * onBlockv[1] * onBlockv[2] * sizeof(int));
    // // allocate the destination tags
    // _i2o_destTag = (int*)flups_malloc(inBlockv[0] * inBlockv[1] * inBlockv[2] * sizeof(int));
    // _o2i_destTag = (int*)flups_malloc(onBlockv[0] * onBlockv[1] * onBlockv[2] * sizeof(int));

    // // get the size of the blocks
    // _cmpt_blockSize(inBlockv, iblockIDStart, nByBlock, istart, iend, _iBlockSize);
    // _cmpt_blockSize(onBlockv, oblockIDStart, nByBlock, ostart, oend, _oBlockSize);

    // get the ranks
    _cmpt_blockDestRank(inBlockv,nByBlock,_shift,istart,topo_in,topo_out,_i2o_destRank);
    _cmpt_blockDestRank(onBlockv,nByBlock,mshift,ostart,topo_out,topo_in,_o2i_destRank);
    // _cmpt_blockDestRankAndTag(inBlockv, iblockIDStart, topo_out, ostartBlockEachProc, onBlockEachProc, _i2o_destRank, _i2o_destTag);
    // _cmpt_blockDestRankAndTag(onBlockv, oblockIDStart, topo_in, istartBlockEachProc, inBlockEachProc, _o2i_destRank,_o2i_destTag);

    // try to gather blocks together if possible, rewrittes the sizes, the blockistart, the number of blocks, the ranks and the tags
    _gather_blocks(topo_in, nByBlock, istart, iend, inBlockv, _iBlockSize, _iBlockiStart, &_inBlock, &_i2o_destRank);
    _gather_blocks(topo_out, nByBlock, ostart, oend, onBlockv, _oBlockSize, _oBlockiStart, &_onBlock, &_o2i_destRank);
    // get the tags for the gathered blocks
    _gather_tags(_inComm, _inBlock, _onBlock, _i2o_destRank, _o2i_destRank, &_i2o_destTag, &_o2i_destTag);

    // allocate the requests
    _i2o_sendRequest = (MPI_Request*)flups_malloc(_inBlock * sizeof(MPI_Request));
    _i2o_recvRequest = (MPI_Request*)flups_malloc(_onBlock * sizeof(MPI_Request));
    _o2i_sendRequest = (MPI_Request*)flups_malloc(_onBlock * sizeof(MPI_Request));
    _o2i_recvRequest = (MPI_Request*)flups_malloc(_inBlock * sizeof(MPI_Request));

    // free the temp arrays
    // flups_free(inBlockEachProc);
    // flups_free(onBlockEachProc);
    // flups_free(istartBlockEachProc);
    // flups_free(ostartBlockEachProc);

    END_FUNC;
}

void SwitchTopo_nb::_free_blockInfo(){

    if (_i2o_destRank != NULL) flups_free(_i2o_destRank);
    if (_o2i_destRank != NULL) flups_free(_o2i_destRank);
    if (_i2o_destTag != NULL) flups_free(_i2o_destTag);
    if (_o2i_destTag != NULL) flups_free(_o2i_destTag);

    _i2o_destRank = NULL;
    _o2i_destRank = NULL;
    _i2o_destTag  = NULL;
    _o2i_destTag  = NULL;

    for (int id = 0; id < 3; id++) {
        if (_iBlockSize[id] != NULL) flups_free(_iBlockSize[id]);
        if (_oBlockSize[id] != NULL) flups_free(_oBlockSize[id]);
        if (_iBlockiStart[id] != NULL) flups_free(_iBlockiStart[id]);
        if (_oBlockiStart[id] != NULL) flups_free(_oBlockiStart[id]);
        _iBlockSize[id]   = NULL;
        _oBlockSize[id]   = NULL;
        _iBlockiStart[id] = NULL;
        _oBlockiStart[id] = NULL;
    }

    if (_i2o_sendRequest != NULL) flups_free(_i2o_sendRequest);
    if (_i2o_recvRequest != NULL) flups_free(_i2o_recvRequest);
    if (_o2i_sendRequest != NULL) flups_free(_o2i_sendRequest);
    if (_o2i_recvRequest != NULL) flups_free(_o2i_recvRequest);

    _i2o_sendRequest = NULL;
    _i2o_recvRequest = NULL;
    _o2i_sendRequest = NULL;
    _o2i_recvRequest = NULL;

}

/**
 * @brief 
 * 
 */
void SwitchTopo_nb::setup(){
    BEGIN_FUNC;

    int rank, comm_size;
    MPI_Comm inComm = _topo_in->get_comm();
    MPI_Comm outComm = _topo_out->get_comm();

    MPI_Comm_rank(inComm, &rank);
    MPI_Comm_size(inComm, &comm_size);

    //Ensure that comms have not changed since init. Otherwise recompute the source/destination of blocks.
    int compIn, compOut;
    MPI_Comm_compare(inComm, _inComm, &compIn);
    MPI_Comm_compare(outComm, _outComm, &compOut);
    if( compIn != MPI_IDENT || compOut != MPI_IDENT){
        FLUPS_WARNING("The inComm and/or outComm have changed since this switchtopo was created. I will recompute the communication scheme.",LOCATION);

        _inComm = inComm;
        _outComm = outComm;

        //reinit the block information
        _free_blockInfo();

        //The input topo may have been reset to real, even if this switchtopo is a complex2complex. 
        //We create a tmp input topo which is complex if needed, for the computation of start and end.
        bool isC2C = _topo_out->isComplex();
        
        int tmp_nglob[3], tmp_nproc[3], tmp_axproc[3];
        for(int i = 0; i<3;i++){
            tmp_nglob[i] = _topo_in->nglob(i);
            tmp_nproc[i] = _topo_in->nproc(i);
            tmp_axproc[i] = _topo_in->axproc(i);
        }
        const Topology* topo_in_tmp = new Topology(_topo_in->axis(),tmp_nglob,tmp_nproc,isC2C,tmp_axproc,FLUPS_ALIGNMENT,_topo_in->get_comm());

        //recompute block info
        _init_blockInfo(topo_in_tmp, _topo_out);

        delete(topo_in_tmp);
    }

    //-------------------------------------------------------------------------
    /** - Setup subcomm */
    //-------------------------------------------------------------------------
    _cmpt_commSplit();
    // setup the dest rank, counts and starts
    _setup_subComm(_inBlock,_iBlockSize, _i2o_destRank,NULL,NULL);
    _setup_subComm(_onBlock,_oBlockSize, _o2i_destRank,NULL,NULL);

    //-------------------------------------------------------------------------
    /** - Compute the self blocks in the new comms   */
    //-------------------------------------------------------------------------
    int newrank;
    MPI_Comm_rank(_subcomm,&newrank);
    _selfBlockN = 0;
    for (int bid = 0; bid < _inBlock; bid++) {
        // for the send when doing input 2 output: send to rank i2o with tag _i2o_destTag[bid]
        if (_i2o_destRank[bid] == newrank) {
            _selfBlockN++;
        }
    }
    int temp = 0;
    for (int bid = 0; bid < _onBlock; bid++) {
        if (_o2i_destRank[bid] == newrank) {
            temp++;
        }
    }
    FLUPS_CHECK(temp == _selfBlockN, "the number of selfBlocks has to be the same in both TOPO!", LOCATION);
    _iselfBlockID = (int*)flups_malloc(_selfBlockN * sizeof(int));
    _oselfBlockID = (int*)flups_malloc(_selfBlockN * sizeof(int));
    //-------------------------------------------------------------------------
    /** - Display performance information if asked */
    //-------------------------------------------------------------------------
#ifdef PERF_VERBOSE
    int rankworld;
    MPI_Comm_rank(inComm, &rankworld);
    // we display important information for the performance
    string name = "./prof/SwitchTopo_" + std::to_string(_topo_in->axis()) + "to" + std::to_string(_topo_out->axis()) + "_rank" + std::to_string(rankworld) + ".txt";
    FILE* file = fopen(name.c_str(),"a+");
    if(file != NULL){
        fprintf(file,"============================================================\n");
+       fprintf(file,"NX = %d - rank = %d - threads = %d\n",_topo_in->nglob(0),comm_size,omp_get_max_threads());
        fprintf(file,"- non blocking and persistent communications\n");

        // int rlen;
        // char myname[MPI_MAX_OBJECT_NAME];
        // MPI_Comm_get_name(_subcomm, myname, &rlen);
        // fprintf(file,"- in subcom %s with rank %d/%d\n",myname,newrank,subsize);
        fprintf(file,"- nglob = %d %d %d to %d %d %d\n",_topo_in->nglob(0),_topo_in->nglob(1),_topo_in->nglob(2),_topo_out->nglob(0),_topo_out->nglob(1),_topo_out->nglob(2));
        fprintf(file,"- nproc = %d %d %d to %d %d %d\n",_topo_in->nproc(0),_topo_in->nproc(1),_topo_in->nproc(2),_topo_out->nproc(0),_topo_out->nproc(1),_topo_out->nproc(2));
        // fprintf(file,"- start = %d %d %d to %d %d %d\n",_istart[0],_istart[1],_istart[2],_ostart[0],_ostart[1],_ostart[2]);
        // fprintf(file,"- end = %d %d %d to %d %d %d\n",_iend[0],_iend[1],_iend[2],_oend[0],_oend[1],_oend[2]);
        // int totalsize = (_nByBlock[0]+_exSize[0]%2)*(_nByBlock[1]+_exSize[1]%2)*(_nByBlock[2]+_exSize[2]%2)*_topo_out->nf();
        // fprintf(file,"- nByBlock = %d %d %d, real size = %d %d %d, alignement padding? %d vs %d\n",_nByBlock[0],_nByBlock[1],_nByBlock[2],(_nByBlock[0] == 1)?1:_nByBlock[0]+_exSize[0]%2,(_nByBlock[1] == 1)?1:_nByBlock[1]+_exSize[1]%2,(_nByBlock[2] == 1)?1:_nByBlock[2]+_exSize[2]%2,totalsize,get_blockMemSize());

        fprintf(file,"--------------------------\n");
        fprintf(file,"%d inblock: %d\n",newrank,_inBlock);
        fprintf(file,"%d I2O:\n",newrank);
        for(int ib=0; ib<_inBlock; ib++){
            fprintf(file," block[%d]: size= %d %d %d - istart = %d %d %d - destrank = %d - desttag = %d\n",ib,_iBlockSize[0][ib],_iBlockSize[1][ib],_iBlockSize[2][ib],_iBlockiStart[0][ib],_iBlockiStart[1][ib],_iBlockiStart[2][ib],_i2o_destRank[ib],_i2o_destTag[ib]);
        }
        fprintf(file,"\n");
        fprintf(file,"--------------------------\n");
        fprintf(file,"%d onblock: %d\n",newrank,_onBlock);
        fprintf(file,"%d O2I:\n",newrank);
        for(int ib=0; ib<_onBlock; ib++){
            fprintf(file," block[%d]: size= %d %d %d - istart = %d %d %d - destrank = %d - desttag = %d \n",ib,_oBlockSize[0][ib],_oBlockSize[1][ib],_oBlockSize[2][ib],_oBlockiStart[0][ib],_oBlockiStart[1][ib],_oBlockiStart[2][ib],_o2i_destRank[ib],_o2i_destTag[ib]);
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
    const int nf = std::max(_topo_in->nf(),_topo_out->nf());

    int newrank;
    MPI_Comm_rank(_subcomm, &newrank);
    // allocate the second layer of buffers
    _sendBuf = (double**)flups_malloc(_inBlock * sizeof(double*));
    _recvBuf = (double**)flups_malloc(_onBlock * sizeof(double*));

    const bool doShuffle=(_topo_in->axis() != _topo_out->axis());
    
    if (doShuffle) {
        _i2o_shuffle = (fftw_plan*)flups_malloc(_onBlock * sizeof(fftw_plan));
        _o2i_shuffle = (fftw_plan*)flups_malloc(_inBlock * sizeof(fftw_plan));
    } else {
        _i2o_shuffle = NULL;
        _o2i_shuffle = NULL;
    }
    
    //-------------------------------------------------------------------------
    /** - for each block we associate the the data buffer and the MPI requests or associate it to NULL */
    //-------------------------------------------------------------------------
    // reset the counter to 0
    int selfcount = 0;
    for (int bid = 0; bid < _inBlock; bid++) {
        //associate the pointer to the correct block
        _sendBuf[bid] = sendData;
        sendData = sendData + get_blockMemSize(bid,nf,_iBlockSize);
        // for the send when doing input 2 output: send to rank i2o with tag _i2o_destTag[bid]
        if (_i2o_destRank[bid] == newrank) {
            // save the bid to the self block list
            _iselfBlockID[selfcount] = bid;
            // associate the request to NULL
            _i2o_sendRequest[bid] = MPI_REQUEST_NULL;
            _o2i_recvRequest[bid] = MPI_REQUEST_NULL;
            // increment the counter
            selfcount++;
        } else {
            // get the send size without padding
            const size_t sendSize = (size_t)_iBlockSize[0][bid] * (size_t)_iBlockSize[1][bid] * (size_t)_iBlockSize[2][bid] * (size_t)(nf);
            MPI_Send_init(_sendBuf[bid], sendSize, MPI_DOUBLE, _i2o_destRank[bid], _i2o_destTag[bid], _subcomm, &(_i2o_sendRequest[bid]));
            // for the send when doing output 2 input: send to rank o2i with tag o2i
            MPI_Recv_init(_sendBuf[bid], sendSize, MPI_DOUBLE, _i2o_destRank[bid], bid, _subcomm, &(_o2i_recvRequest[bid]));
        }

        // setup the suffle plan for the out 2 in transformation if needed
        if (doShuffle) {
            int tmp_size[3] = {_iBlockSize[0][bid], _iBlockSize[1][bid], _iBlockSize[2][bid]};
            _setup_shuffle(tmp_size, _topo_out, _topo_in, _sendBuf[bid], &_o2i_shuffle[bid]);
            FLUPS_INFO("doing a shuffle allocation");
        }
    }
    FLUPS_CHECK(selfcount == _selfBlockN,"the number of counted block has to match the allocation number: %d vs %d",selfcount,_selfBlockN,LOCATION);

    // reset the self count
    selfcount = 0;
    for (int bid = 0; bid < _onBlock; bid++) {
        //associate the pointer with the correct block
        _recvBuf[bid] = recvData;
        recvData = recvData+ get_blockMemSize(bid,nf,_oBlockSize);
        // create the request if needed
        if (_o2i_destRank[bid] == newrank) {
            // save the bid
            _oselfBlockID[selfcount] = bid;
            // associate the request to NULL
            _i2o_recvRequest[bid] = MPI_REQUEST_NULL;
            _o2i_sendRequest[bid] = MPI_REQUEST_NULL;
            // increment the counter
            selfcount++;
        } else {
            // get the receive size without padding
            const size_t recvSize = (size_t)_oBlockSize[0][bid] * (size_t)_oBlockSize[1][bid] * (size_t)_oBlockSize[2][bid] * (size_t)_topo_out->nf();
            // for the reception when doing input 2 output: receive from the rank o2i with tag bid
            MPI_Recv_init(_recvBuf[bid], recvSize, MPI_DOUBLE, _o2i_destRank[bid], bid, _subcomm, &(_i2o_recvRequest[bid]));
            // for the send when doing output 2 input: send to rank o2i with tag o2i
            MPI_Send_init(_recvBuf[bid], recvSize, MPI_DOUBLE, _o2i_destRank[bid], _o2i_destTag[bid], _subcomm, &(_o2i_sendRequest[bid]));
        }

        // setup the suffle plan for the in 2 out transformation
        if (doShuffle) {
            int tmp_size[3] = {_oBlockSize[0][bid], _oBlockSize[1][bid], _oBlockSize[2][bid]};
            _setup_shuffle(tmp_size,_topo_in,_topo_out, _recvBuf[bid], &_i2o_shuffle[bid]);
            FLUPS_INFO("doing a shuffle allocation");
        }
    }
    FLUPS_CHECK(selfcount == _selfBlockN,"the number of counted block has to match the allocation number: %d vs %d",selfcount,_selfBlockN,LOCATION);
    END_FUNC;
}

/**
 * @brief Destroy the Switch Topo
 * 
 */
SwitchTopo_nb::~SwitchTopo_nb() {
    BEGIN_FUNC;

    int  rlen, comp;
    MPI_Comm_compare(_subcomm,_inComm,&comp);
    if(comp!=MPI_IDENT){
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(_subcomm, myname, &rlen);
        FLUPS_INFO("freeing the comm %s",myname);
        MPI_Comm_free(&_subcomm);
    }

    for(int ib=0; ib< _inBlock; ib++){
        // if (_sendBuf[ib] != NULL) flups_free(_sendBuf[ib]);
        if (_i2o_sendRequest[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(_i2o_sendRequest[ib]));
        if (_o2i_recvRequest[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(_o2i_recvRequest[ib]));
    }
    for(int ib=0; ib< _onBlock; ib++){
        // if (_recvBuf[ib] != NULL) flups_free(_recvBuf[ib]);
        if (_i2o_recvRequest[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(_i2o_recvRequest[ib]));
        if (_o2i_sendRequest[ib] != MPI_REQUEST_NULL) MPI_Request_free(&(_o2i_sendRequest[ib]));
    }

    _free_blockInfo();

    if (_iselfBlockID != NULL) flups_free(_iselfBlockID);
    if (_oselfBlockID != NULL) flups_free(_oselfBlockID);

    if (_i2o_shuffle != NULL) {
        for (int ib = 0; ib < _onBlock; ib++) {
            fftw_destroy_plan(_i2o_shuffle[ib]);
        }
        flups_free(_i2o_shuffle);
    }
    if (_o2i_shuffle != NULL) {
        for (int ib = 0; ib < _inBlock; ib++) {
            fftw_destroy_plan(_o2i_shuffle[ib]);
        }
        flups_free(_o2i_shuffle);
    }

    flups_free((double*)_sendBuf);
    flups_free((double*)_recvBuf);
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

    FLUPS_CHECK(_topo_in->isComplex() == _topo_out->isComplex(),"both topologies have to be complex or real", LOCATION);
    FLUPS_CHECK(_topo_in->nf() <= 2, "the value of nf is not supported", LOCATION);

    PROF_START("reorder");
    int iswitch = _iswitch;

    //-------------------------------------------------------------------------
    /** - setup required memory arrays */
    //-------------------------------------------------------------------------

    const Topology* topo_in;
    const Topology* topo_out;

    MPI_Request* sendRequest;
    MPI_Request* recvRequest;

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
    
    int* oselfBlockID;
    int* destTag;

    // const int nByBlock[3] = {_nByBlock[0],_nByBlock[1],_nByBlock[2]};

    opt_double_ptr* sendBuf;
    opt_double_ptr* recvBuf;

    fftw_plan* shuffle = NULL;

    if (sign == FLUPS_FORWARD) {
        topo_in     = _topo_in;
        topo_out    = _topo_out;
        sendRequest = _i2o_sendRequest;
        recvRequest = _i2o_recvRequest;
        sendBuf     = _sendBuf;
        recvBuf     = _recvBuf;

        oselfBlockID = _oselfBlockID;
        destTag      = _i2o_destTag;

        shuffle = _i2o_shuffle;

        send_nBlock = _inBlock;
        recv_nBlock = _onBlock;

        for (int id = 0; id < 3; id++) {
            // send_nBlock[id] = _inBlock[id];
            // recv_nBlock[id] = _onBlock[id];
            // istart[id]      = _istart[id];
            // iend[id]        = _iend[id];
            // ostart[id]      = _ostart[id];
            // oend[id]        = _oend[id];
            inmem[id]       = _topo_in->nmem(id);
            onmem[id]       = _topo_out->nmem(id);
            iBlockSize[id]  = _iBlockSize[id];
            oBlockSize[id]  = _oBlockSize[id];
            iBlockiStart[id]  = _iBlockiStart[id];
            oBlockiStart[id]  = _oBlockiStart[id];
        }
    } else if (sign == FLUPS_BACKWARD) {
        topo_in     = _topo_out;
        topo_out    = _topo_in;
        sendRequest = _o2i_sendRequest;
        recvRequest = _o2i_recvRequest;
        sendBuf     = _recvBuf;
        recvBuf     = _sendBuf;

        oselfBlockID = _iselfBlockID;
        destTag      = _o2i_destTag;

        shuffle = _o2i_shuffle;

        send_nBlock = _onBlock;
        recv_nBlock = _inBlock;

        for (int id = 0; id < 3; id++) {
            // send_nBlock[id] = _onBlock[id];
            // recv_nBlock[id] = _inBlock[id];
            // istart[id]      = _ostart[id];
            // iend[id]        = _oend[id];
            // ostart[id]      = _istart[id];
            // oend[id]        = _iend[id];
            inmem[id]       = _topo_out->nmem(id);
            onmem[id]       = _topo_in->nmem(id);
            iBlockSize[id]  = _oBlockSize[id];
            oBlockSize[id]  = _iBlockSize[id];
            iBlockiStart[id]  = _oBlockiStart[id];
            oBlockiStart[id]  = _iBlockiStart[id];
            
        }
    } else {
        FLUPS_CHECK(false, "the sign is not FLUPS_FORWARD nor FLUPS_BACKWARD", LOCATION);
    }

    FLUPS_INFO("previous topo: %d,%d,%d axis=%d", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    FLUPS_INFO("new topo: %d,%d,%d  axis=%d", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());
    FLUPS_INFO("using %d blocks on send and %d on recv",send_nBlock,recv_nBlock);

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

    PROF_STARTi("switch",iswitch);
    PROF_STARTi("mem2buf",iswitch);
    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    // const int nblocks_send = send_nBlock[0] * send_nBlock[1] * send_nBlock[2];

#if defined(__INTEL_COMPILER)
#pragma omp parallel proc_bind(close) default(none) firstprivate(send_nBlock, v, sendBuf, recvBuf, destTag, iBlockSize,iBlockiStart, nf, inmem, iax0, iax1,iax2,sendRequest)
#elif defined(__GNUC__)
#pragma omp parallel proc_bind(close) default(none) shared(ompi_request_null) firstprivate(send_nBlock, v, sendBuf, recvBuf, destTag,iBlockSize,iBlockiStart, nf, inmem, iax0, iax1,iax2,sendRequest)
#endif
    for (int bid = 0; bid < send_nBlock; bid++) {
        // // get the split index
        // int ib[3];
        // localSplit(bid, send_nBlock, 0, ib, 1);
        // get the buffer data for this block
        double* data;
        if(sendRequest[bid] == MPI_REQUEST_NULL){
            // if we are doing a self block the data is the recv buff
            // the new block ID is given by destTag[bid]
            data = recvBuf[destTag[bid]];
        } else {
            // else we copy inside the sendbuffer
            data = sendBuf[bid];
        }
        // // get the starting index in the global memory
        // const int loci0 = istart[iax0] + ib[iax0] * nByBlock[iax0];
        // const int loci1 = istart[iax1] + ib[iax1] * nByBlock[iax1];
        // const int loci2 = istart[iax2] + ib[iax2] * nByBlock[iax2];

        // go inside the block
        const int id_max = iBlockSize[iax1][bid] * iBlockSize[iax2][bid];
        const size_t nmax = iBlockSize[iax0][bid] * nf;

        // the buffer is aligned if the starting id is aligned and if nmax is a multiple of the alignement
        const bool isBuffAligned = FLUPS_ISALIGNED(data) &&  nmax%FLUPS_ALIGNMENT == 0;
        // the data is aligned if the starting index is aligned AND if the gap between two entries, inmem[iax0] is a multiple of the alignment
        double*    my_v            = v + localIndex(iax0, iBlockiStart[iax0][bid], iBlockiStart[iax1][bid], iBlockiStart[iax2][bid], iax0, inmem, nf);
        const bool isVectorAligned = FLUPS_ISALIGNED(my_v) && inmem[iax0] % FLUPS_ALIGNMENT == 0;

        // we choose the best loop depending on the alignement
        if (isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
            for (int id = 0; id < id_max; id++) {
                // get the id from a small modulo
                const int i2 = id / iBlockSize[iax1][bid];
                const int i1 = id % iBlockSize[iax1][bid];

                // get the local starting location for the buffer and the field
                const opt_double_ptr vloc = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf);
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
                const double* __restrict vloc = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf);
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
                const opt_double_ptr vloc  = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf);
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
                const double* __restrict vloc  = my_v + localIndex(iax0, 0, i1, i2, iax0, inmem, nf);
                double* __restrict dataloc = data + id * nmax;

                // do the copy -> vectorized
                for (size_t i0 = 0; i0 < nmax; i0++) {
                    dataloc[i0] = vloc[i0];
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

    PROF_STOPi("mem2buf",iswitch);

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

#pragma omp parallel default(none) proc_bind(close) shared(status) firstprivate(recv_nBlock, oselfBlockID, v, recvBuf, oBlockSize, oBlockiStart, nf, onmem, oax0, oax1, oax2, recvRequest, iswitch, shuffle)
    for (int count = 0; count < recv_nBlock; count++) {
        // only the master receive the call
        int bid = -1;
        // if we are doing a self block
        if (count < _selfBlockN) {
            bid = oselfBlockID[count];
#pragma omp master
            {
                PROF_STARTi("buf2mem",iswitch);
                // only the master call the fftw_execute which is executed in multithreading
                if (shuffle != NULL) {
                    fftw_execute(shuffle[bid]);
                }
            }
#pragma omp barrier
        } else {
#pragma omp master
            {
                PROF_STARTi("waiting",iswitch);
                int request_index;
                MPI_Waitany(recv_nBlock, recvRequest, &request_index, &status);
                PROF_STOPi("waiting",iswitch);
                PROF_STARTi("buf2mem",iswitch);
                
                // bid is set for the master
                bid = status.MPI_TAG;
                // only the master call the fftw_execute which is executed in multithreading                
                if (shuffle != NULL) {
                    fftw_execute(shuffle[bid]);
                }
            }
            // make sure that the master has received the status before going further
            // there is no implicit barrier after
#pragma omp barrier
            // get the block id = the tag
            bid = status.MPI_TAG;
        }
        // add the bandwith info
        #pragma omp master
        {
#ifdef PROF            
            if (_prof != NULL) {
                _prof->addMem("waiting"+to_string(iswitch), get_blockMemSize(bid,nf,oBlockSize)*sizeof(double));
            }
#endif
        }
        
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
        const bool isBuffAligned = FLUPS_ISALIGNED(recvBuf[bid]) &&  nmax%FLUPS_ALIGNMENT == 0;
        // the data is aligned if the starting index is aligned AND if the gap between two entries, inmem[iax0] is a multiple of the alignment
        // double*    my_v            = v + localIndex(oax0, oBlockiStart[0][bid], oBlockiStart[1][bid], oBlockiStart[2][bid], oax0, onmem, nf);
        double*    my_v            = v + localIndex(oax0, oBlockiStart[oax0][bid], oBlockiStart[oax1][bid], oBlockiStart[oax2][bid], oax0, onmem, nf);
        const bool isVectorAligned = FLUPS_ISALIGNED(my_v) && onmem[oax0] % FLUPS_ALIGNMENT == 0;

        //choose the correct loop to improve the efficiency
        if (isBuffAligned && isVectorAligned) {
#pragma omp for schedule(static)
            for (int id = 0; id < id_max; id++) {
                // get the id from a small modulo
                const int i2 = id / oBlockSize[oax1][bid];
                const int i1 = id % oBlockSize[oax1][bid];

                // get the local starting id for the buffer and the data
                opt_double_ptr       vloc    = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf);
                const opt_double_ptr dataloc = recvBuf[bid] + id * nmax;
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
                double* __restrict vloc      = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf);
                const opt_double_ptr dataloc = recvBuf[bid] + id * nmax;
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
                opt_double_ptr vloc              = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf);
                const double* __restrict dataloc = recvBuf[bid] + id * nmax;
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
                double* __restrict vloc          = my_v + localIndex(oax0, 0, i1, i2, oax0, onmem, nf);
                const double* __restrict dataloc = recvBuf[bid] + id * nmax;
                // do the copy
                for (size_t i0 = 0; i0 < nmax; i0++) {
                    vloc[i0] = dataloc[i0];
                }
            }
        }

#pragma omp master
        {
            PROF_STOPi("buf2mem",iswitch);
        }
    }
    // now that we have received everything, close the send requests
    MPI_Waitall(send_nBlock, sendRequest,MPI_STATUSES_IGNORE);

    PROF_STOPi("switch",iswitch);
    PROF_STOP("reorder");
    END_FUNC;
}

void SwitchTopo_nb::disp() const {
    BEGIN_FUNC;
    FLUPS_INFO("------------------------------------------");
    FLUPS_INFO("## Topo Swticher MPI");
    FLUPS_INFO("--- INPUT");
    FLUPS_INFO("  - input axis = %d", _topo_in->axis());
    FLUPS_INFO("  - input local = %d %d %d", _topo_in->nloc(0), _topo_in->nloc(1), _topo_in->nloc(2));
    FLUPS_INFO("  - input global = %d %d %d", _topo_in->nglob(0), _topo_in->nglob(1), _topo_in->nglob(2));
    // FLUPS_INFO("  - istart = %d %d %d", _istart[0], _istart[1], _istart[2]);
    // FLUPS_INFO("  - iend = %d %d %d", _iend[0], _iend[1], _iend[2]);
    FLUPS_INFO("--- OUTPUT");
    FLUPS_INFO("  - output axis = %d", _topo_out->axis());
    FLUPS_INFO("  - output local = %d %d %d", _topo_out->nloc(0), _topo_out->nloc(1), _topo_out->nloc(2));
    FLUPS_INFO("  - output global = %d %d %d", _topo_out->nglob(0), _topo_out->nglob(1), _topo_out->nglob(2));
    // FLUPS_INFO("  - ostart = %d %d %d", _ostart[0], _ostart[1], _ostart[2]);
    // FLUPS_INFO("  - oend = %d %d %d", _oend[0], _oend[1], _oend[2]);
    FLUPS_INFO("--- BLOCKS");
    FLUPS_INFO("  - selfBlockN = %d", _selfBlockN);
    // FLUPS_INFO("  - nByBlock  = %d %d %d", _nByBlock[0], _nByBlock[1], _nByBlock[2]);
    FLUPS_INFO("  - inBlock = %d", _inBlock);
    FLUPS_INFO("  - onBlock = %d", _onBlock);
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
    Topology* topo    = new Topology(0, nglob, nproc, false,NULL,1, MPI_COMM_WORLD);
    Topology* topobig = new Topology(0, nglob_big, nproc_big, false,NULL,1, MPI_COMM_WORLD);

    double* data = (double*)flups_malloc(sizeof(double*) * std::max(topo->memsize(), topobig->memsize()));

    const int nmem[3] = {topo->nmem(0),topo->nmem(1),topo->nmem(2)};
    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localIndex(0,i0, i1, i2,0,nmem,1);
                data[id] = id;
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_real", data);

    const int fieldstart[3] = {0, 0, 0};
    // printf("\n=============================");
    SwitchTopo* switchtopo = new SwitchTopo_nb(topo, topobig, fieldstart, NULL);
    switchtopo->setup();
    switchtopo->disp();

    size_t         max_mem    = switchtopo->get_bufMemSize();
    opt_double_ptr send_buff  = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
    opt_double_ptr recv_buff  = (opt_double_ptr)flups_malloc(max_mem * sizeof(double));
    std::memset(send_buff, 0, max_mem * sizeof(double));
    std::memset(recv_buff, 0, max_mem * sizeof(double));
    // associate the buffer
    
    switchtopo->setup_buffers(send_buff, recv_buff);

    // printf("\n\n============ FORWARD =================");
    switchtopo->execute(data, FLUPS_FORWARD);

    hdf5_dump(topobig, "test_real_padd", data);

    // printf("\n\n============ BACKWARD =================");
    switchtopo->execute(data, FLUPS_BACKWARD);

    hdf5_dump(topo, "test_real_returned", data);

    flups_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);

    //===========================================================================
    // complex numbers
    topo    = new Topology(0, nglob, nproc, true,NULL,1, MPI_COMM_WORLD);
    topobig = new Topology(2, nglob_big, nproc_big, true,NULL,1, MPI_COMM_WORLD);

    data = (double*)flups_malloc(sizeof(double*) * topobig->memsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localIndex(0,i0, i1, i2,0,nmem,2);
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

    const int fieldstart2[3] = {4, 0, 0};
    // printf("\n=============================");
    switchtopo = new SwitchTopo_nb(topo, topobig, fieldstart2, NULL);

    switchtopo->execute(data, FLUPS_FORWARD);

    hdf5_dump(topobig, "test_complex_padd", data);

    switchtopo->execute(data, FLUPS_BACKWARD);

    hdf5_dump(topo, "test_complex_returned", data);

    flups_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);
    END_FUNC;
}