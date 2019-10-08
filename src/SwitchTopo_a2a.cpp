/**
 * @file SwitchTopo_a2a.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-25
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
using namespace FLUPS;

SwitchTopo_a2a::SwitchTopo_a2a(const Topology* topo_input, const Topology* topo_output, const int shift[3], Profiler* prof) {
    BEGIN_FUNC;

    FLUPS_CHECK(topo_input->isComplex() == topo_output->isComplex(), "both topologies have to be the same kind", LOCATION);

    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    _topo_in  = topo_input;
    _topo_out = topo_output;
    _prof     = prof;

    //-------------------------------------------------------------------------
    /** - get the starting and ending index of the shared zone */
    //-------------------------------------------------------------------------
    // get the blockshift
    _topo_in->cmpt_intersect_id(shift, _topo_out, _istart, _iend);
    int tmp[3] = {-shift[0], -shift[1], -shift[2]};
    _topo_out->cmpt_intersect_id(tmp, _topo_in, _ostart, _oend);

    //-------------------------------------------------------------------------
    /** - get the block size as the GCD of the memory among every process between send and receive */
    //-------------------------------------------------------------------------
    _cmpt_nByBlock();

    //-------------------------------------------------------------------------
    /** - get the number of blocks and for each block get the size and the destination rank */
    //-------------------------------------------------------------------------
    int  iblockIDStart[3];
    int  oblockIDStart[3];
    int* inBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));
    int* onBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));

    _cmpt_blockIndexes(_istart, _iend, _nByBlock, _topo_in, _inBlock, iblockIDStart, inBlockEachProc);
    _cmpt_blockIndexes(_ostart, _oend, _nByBlock, _topo_out, _onBlock, oblockIDStart, onBlockEachProc);

    // allocte the block size
    for (int id = 0; id < 3; id++) {
        _iBlockSize[id] = (int*)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
        _oBlockSize[id] = (int*)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));
    }

    // allocate the destination ranks
    _i2o_destRank = (opt_int_ptr)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
    _o2i_destRank = (opt_int_ptr)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));

    // get the send destination ranks in the ouput topo
    _cmpt_blockSize(_inBlock, iblockIDStart, _nByBlock, _istart, _iend, _iBlockSize);
    _cmpt_blockSize(_onBlock, oblockIDStart, _nByBlock, _ostart, _oend, _oBlockSize);

    _cmpt_blockDestRankAndTag(_inBlock, iblockIDStart, _topo_out, onBlockEachProc, _i2o_destRank,NULL);
    _cmpt_blockDestRankAndTag(_onBlock, oblockIDStart, _topo_in, inBlockEachProc, _o2i_destRank,NULL);

    // free the temp arrays
    fftw_free(inBlockEachProc);
    fftw_free(onBlockEachProc);

    //-------------------------------------------------------------------------
    /** - Setup subcomm */
    //-------------------------------------------------------------------------
    _cmpt_commSplit();
    // setup the dest rank, counts and starts
    _setup_subComm(_subcomm, _inBlock, _i2o_destRank, &_i2o_count, &_i2o_start);
    _setup_subComm(_subcomm, _onBlock, _o2i_destRank, &_o2i_count, &_o2i_start);

    //-------------------------------------------------------------------------
    /** - determine if we are all to all */
    //-------------------------------------------------------------------------
    int subsize;
    MPI_Comm_size(_subcomm, &subsize);

    int tmp_size = _i2o_count[0];
    _is_all2all  = (tmp_size != 0);
    _is_all2all  = _is_all2all && (tmp_size == _o2i_count[0]);
    for (int ir = 1; ir < subsize; ir++) {
        // if the count from and to the rank is the same, we can do an A2A
        _is_all2all = _is_all2all && (tmp_size == _i2o_count[ir]);
        _is_all2all = _is_all2all && (tmp_size == _o2i_count[ir]);
    }

    // if we are all to all, clean the start array
    if (_is_all2all) {
        if (_i2o_start != NULL) {
            fftw_free(_i2o_start);
            _i2o_start = NULL;
        }
        if (_o2i_start != NULL) {
            fftw_free(_o2i_start);
            _o2i_start = NULL;
        }
    }

    //-------------------------------------------------------------------------
    /** - initialize the profiler    */
    //-------------------------------------------------------------------------
    if (_prof != NULL) {
        _prof->create("reorder", "solve");

        _prof->create("switch0", "reorder");
        _prof->create("0mem2buf", "switch0");
        _prof->create("0buf2mem", "switch0");
        _prof->create("0all_2_all", "switch0");
        _prof->create("0all_2_all_v", "switch0");

        _prof->create("switch1", "reorder");
        _prof->create("1mem2buf", "switch1");
        _prof->create("1buf2mem", "switch1");
        _prof->create("1all_2_all", "switch1");
        _prof->create("1all_2_all_v", "switch1");

        _prof->create("switch2", "reorder");
        _prof->create("2mem2buf", "switch2");
        _prof->create("2buf2mem", "switch2");
        _prof->create("2all_2_all", "switch2");
        _prof->create("2all_2_all_v", "switch2");
    }

    //-------------------------------------------------------------------------
    /** - Display performance information if asked */
    //-------------------------------------------------------------------------
#ifdef PERF_VERBOSE
    // we display important information for the performance
    string name = "./prof/SwitchTopo_" + std::to_string(_topo_in->axis()) + "to" + std::to_string(_topo_out->axis()) + "_rank" + std::to_string(rank) + ".txt";
    FILE* file = fopen(name.c_str(),"a+");
    if(file != NULL){
        fprintf(file,"============================================================\n");
+       fprintf(file,"NX = %d - rank = %d - threads = %d\n",_topo_in->nglob(0),comm_size,omp_get_max_threads());
        if(_is_all2all) fprintf(file,"- is all to all\n");
        if(!_is_all2all) fprintf(file,"- is all to all VECTOR\n");
        
        int newrank;
        MPI_Comm_rank(_subcomm,&newrank);
        int rlen;
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(_subcomm, myname, &rlen);
        fprintf(file,"- in subcom %s with rank %d/%d\n",myname,newrank,subsize);
        fprintf(file,"- nglob = %d %d %d to %d %d %d\n",_topo_in->nglob(0),_topo_in->nglob(1),_topo_in->nglob(2),_topo_out->nglob(0),_topo_out->nglob(1),_topo_out->nglob(2));
        fprintf(file,"- nproc = %d %d %d to %d %d %d\n",_topo_in->nproc(0),_topo_in->nproc(1),_topo_in->nproc(2),_topo_out->nproc(0),_topo_out->nproc(1),_topo_out->nproc(2));
        int totalsize = (_nByBlock[0]+_exSize[0]%2)*(_nByBlock[1]+_exSize[1]%2)*(_nByBlock[2]+_exSize[2]%2);
        fprintf(file,"- nByBlock = %d %d %d, real size = %d %d %d, alignement padding? %d vs %d\n",_nByBlock[0],_nByBlock[1],_nByBlock[2],(_nByBlock[0] == 1)?1:_nByBlock[0]+_exSize[0]%2,(_nByBlock[1] == 1)?1:_nByBlock[1]+_exSize[1]%2,(_nByBlock[2] == 1)?1:_nByBlock[2]+_exSize[2]%2,totalsize,get_blockMemSize());

        fprintf(file,"--------------------------\n");
        fprintf(file,"%d SEND:",newrank);
        for(int ib=0; ib<_inBlock[0] * _inBlock[1] * _inBlock[2]; ib++){
            fprintf(file," %d ",_i2o_destRank[ib]);
        }
        fprintf(file,"\n");
        fprintf(file,"--------------------------\n");
        fprintf(file,"%d RECV:",newrank);
        for(int ib=0; ib<_onBlock[0] * _onBlock[1] * _onBlock[2]; ib++){
            fprintf(file," %d ",_o2i_destRank[ib]);
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

    int  rlen;
    char myname[MPI_MAX_OBJECT_NAME];
    MPI_Comm_get_name(_subcomm, myname, &rlen);
    FLUPS_INFO("freeing the comm %s",myname);

    MPI_Comm_free(&_subcomm);

    FLUPS_INFO("freeing the arrays");

    if (_i2o_destRank != NULL) fftw_free(_i2o_destRank);
    if (_o2i_destRank != NULL) fftw_free(_o2i_destRank);

    if (_i2o_count != NULL) fftw_free(_i2o_count);
    if (_o2i_count != NULL) fftw_free(_o2i_count);
    if (_i2o_start != NULL) fftw_free(_i2o_start);
    if (_o2i_start != NULL) fftw_free(_o2i_start);

    if (_sendBuf != NULL) fftw_free((double*)_sendBuf);
    if (_recvBuf != NULL) fftw_free((double*)_recvBuf);

    for (int id = 0; id < 3; id++) {
        if (_iBlockSize[id] != NULL) fftw_free(_iBlockSize[id]);
        if (_oBlockSize[id] != NULL) fftw_free(_oBlockSize[id]);
    }
    END_FUNC;
}


/**
 * @brief Setup the communication buffer for every block
 * 
 * We associate #_sendBuf[i_block] with the correct adress inside the raw buffer.
 * This way, we can use only #_sendBuf and #_recvBuf for each block without any additional computation inside the execute.
 * Moreover, asking the user to allocate the data reduces the memory footprint as it can be shared among several SwitchTopo
 * 
 * @param sendData the "raw" communication buffer allocated at least at the size returned by get_bufMemSize 
 * @param recvData the "raw" communication buffer allocated at least at the size returned by get_bufMemSize 
 */
void SwitchTopo_a2a::setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) {
    BEGIN_FUNC;
    int subsize;
    MPI_Comm_size(_subcomm, &subsize);
    // allocate the second layer of buffers
     _sendBuf = (double**)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(double*));
    _recvBuf = (double**)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(double*));
    
    // link the buff of every block to the data initialized
    int*      countPerRank = (int*)fftw_malloc(subsize * sizeof(int));
    const int blockSize    = get_blockMemSize();

    std::memset(countPerRank, 0, subsize * sizeof(int));
    for (int ib = 0; ib < _inBlock[0] * _inBlock[1] * _inBlock[2]; ib++) {
        // get the destination rank
        int destrank = _i2o_destRank[ib];
        // we count the number of blocks on the ranks bellow me
        size_t memblocks = 0;
        for (int ir = 0; ir < destrank; ir++) {
            memblocks += (size_t)_i2o_count[ir];
        }

        // place the block given the number of ranks bellow + the number of block already set to my rank
        _sendBuf[ib] = sendData + memblocks + countPerRank[destrank];
        // add the block size to the number of already added blocks
        countPerRank[destrank] += blockSize;
    }

    std::memset(countPerRank, 0, subsize * sizeof(int));
    for (int ib = 0; ib < _onBlock[0] * _onBlock[1] * _onBlock[2]; ib++) {
        // get the destination rank
        int destrank = _o2i_destRank[ib];
        // we count the number of blocks on the ranks bellow me
        size_t memblocks = 0;
        for (int ir = 0; ir < destrank; ir++) {
            memblocks += (size_t)_o2i_count[ir];
        }
        // place the block given the number of ranks bellow + the number of block already set to my rank
        _recvBuf[ib] = recvData + memblocks + countPerRank[destrank];
        // add the block size to the number of already added blocks
        countPerRank[destrank] += blockSize;
    }

    fftw_free(countPerRank);
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
void SwitchTopo_a2a::execute(opt_double_ptr v, const int sign, const int iswitch) const {
    BEGIN_FUNC;

    FLUPS_CHECK(_topo_in->isComplex() == _topo_out->isComplex(), "both topologies have to be complex or real", LOCATION);
    FLUPS_CHECK(_topo_in->nf() <= 2, "the value of nf is not supported", LOCATION);

    int comm_size;
    // MPI_Comm_rank(_subcomm, &rank);
    MPI_Comm_size(_subcomm, &comm_size);

    string profName;
    if (_prof != NULL) _prof->start("reorder");

    //-------------------------------------------------------------------------
    /** - setup required memory arrays */
    //-------------------------------------------------------------------------

    const Topology* topo_in;
    const Topology* topo_out;

    int send_nBlock[3];
    int recv_nBlock[3];

    int istart[3];
    int ostart[3];
    int iend[3];
    int oend[3];
    int inmem[3];
    int onmem[3];

    int* iBlockSize[3];
    int* oBlockSize[3];

    int* send_count;
    int* recv_count;
    int* send_start;
    int* recv_start;

    const int nByBlock[3] = {_nByBlock[0], _nByBlock[1], _nByBlock[2]};

    opt_double_ptr* sendBuf;
    opt_double_ptr* recvBuf;

    if (sign == FLUPS_FORWARD) {
        topo_in  = _topo_in;
        topo_out = _topo_out;
        sendBuf  = _sendBuf;
        recvBuf  = _recvBuf;

        send_count = _i2o_count;
        recv_count = _o2i_count;
        send_start = _i2o_start;
        recv_start = _o2i_start;

        for (int id = 0; id < 3; id++) {
            send_nBlock[id] = _inBlock[id];
            recv_nBlock[id] = _onBlock[id];
            istart[id]      = _istart[id];
            iend[id]        = _iend[id];
            ostart[id]      = _ostart[id];
            oend[id]        = _oend[id];
            inmem[id]       = _topo_in->nmem(id);
            onmem[id]       = _topo_out->nmem(id);
            iBlockSize[id]  = _iBlockSize[id];
            oBlockSize[id]  = _oBlockSize[id];
        }
    } else if (sign == FLUPS_BACKWARD) {
        topo_in  = _topo_out;
        topo_out = _topo_in;
        sendBuf  = _recvBuf;
        recvBuf  = _sendBuf;

        send_count = _o2i_count;
        recv_count = _i2o_count;
        send_start = _o2i_start;
        recv_start = _i2o_start;

        for (int id = 0; id < 3; id++) {
            send_nBlock[id] = _onBlock[id];
            recv_nBlock[id] = _inBlock[id];
            istart[id]      = _ostart[id];
            iend[id]        = _oend[id];
            ostart[id]      = _istart[id];
            oend[id]        = _iend[id];
            inmem[id]       = _topo_out->nmem(id);
            onmem[id]       = _topo_in->nmem(id);
            iBlockSize[id]  = _oBlockSize[id];
            oBlockSize[id]  = _iBlockSize[id];
        }
    } else {
        FLUPS_CHECK(false, "the sign is not FLUPS_FORWARD nor FLUPS_BACKWARD", LOCATION);
    }

    FLUPS_INFO("previous topo: %d,%d,%d axis=%d", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    FLUPS_INFO("new topo: %d,%d,%d  axis=%d", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());
    FLUPS_INFO("using %d blocks on send and %d on recv", send_nBlock[0] * send_nBlock[1] * send_nBlock[2], recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2]);

    //Testing if the topo in and out are exactly the same. In that case, we just return.
    if(_is_all2all){

        // printf("ALL INFOs: nblocks %d %d %d / %d %d %d\n",send_nBlock[0],send_nBlock[1],send_nBlock[2],recv_nBlock[0],recv_nBlock[1],recv_nBlock[2]);
        // printf("ALL INFOs: istart  %d %d %d / %d %d %d\n",istart[0],istart[1],istart[2],ostart[0],ostart[1],ostart[2]);
        // printf("ALL INFOs: iend    %d %d %d / %d %d %d\n",iend[0],iend[1],iend[2],oend[0],oend[1],oend[2]);
        // printf("ALL INFOs: inmem   %d %d %d / %d %d %d\n",inmem[0],inmem[1],inmem[2],onmem[0],onmem[1],onmem[2]);
        // printf("ALL INFOs: iBSize  %d %d %d / %d %d %d\n",iBlockSize[0],iBlockSize[1],iBlockSize[2],oBlockSize[0],oBlockSize[1],oBlockSize[2]);
        bool condition;
        for (int id=0;id<3;id++){
            condition &= (send_nBlock[id] == recv_nBlock[id]) && \
                         (istart[id] == ostart[id]) && \
                         (iend[id] == oend[id]) && \
                         (inmem[id] == onmem[id]) && \
                         (iBlockSize[id] == oBlockSize[id]);
        }

        if(condition){
            FLUPS_INFO("Skipping switchtopo: in and out topos are the same");
            if (_prof != NULL) {
                _prof->stop("reorder");
            }
            END_FUNC;
            return;
        }
    }


    // define important constants
    const int ax0 = topo_in->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nf  = topo_in->nf();

    if (_prof != NULL) {
        profName = "switch"+to_string(iswitch);
        _prof->start(profName);
        profName = to_string(iswitch) + "mem2buf";
        _prof->start(profName);
    }
    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    const int nblocks_send = send_nBlock[0] * send_nBlock[1] * send_nBlock[2];

#pragma omp parallel proc_bind(close) default(none) firstprivate(nblocks_send, send_nBlock, v, sendBuf, istart, nByBlock, iBlockSize, nf, inmem, ax0, ax1, ax2)
    for (int bid = 0; bid < nblocks_send; bid++) {
        // get the split index
        int ibv[3];
        localSplit(bid, send_nBlock, 0, ibv, 1);

        // get the starting index in the global memory using !!nByBlock!!
        // since only the last block may have a different size
        const int loci0 = istart[ax0] + ibv[ax0] * nByBlock[ax0];
        const int loci1 = istart[ax1] + ibv[ax1] * nByBlock[ax1];
        const int loci2 = istart[ax2] + ibv[ax2] * nByBlock[ax2];
        // get the memory to write to/from
        opt_double_ptr data = sendBuf[bid];
        double* __restrict my_v = v + localIndex(ax0, loci0, loci1, loci2, ax0, inmem, nf);

#pragma omp master
        {   // reset of the buffer, only done by the master
            std::memset(data, 0, iBlockSize[0][bid] * iBlockSize[1][bid] * iBlockSize[2][bid] * nf * sizeof(double));
        }
        // wait till the reset is over before doing the fill
#pragma omp barrier

        // go inside the block
        const int id_max = iBlockSize[ax1][bid] * iBlockSize[ax2][bid];
#pragma omp for schedule(static)
        for (int id = 0; id < id_max; id++) {
            // get the id from a small modulo
            const int i2 = id / iBlockSize[ax1][bid];
            const int i1 = id % iBlockSize[ax1][bid];
            // get the max counter
            const size_t nmax = iBlockSize[ax0][bid] * nf;
            // get the starting global id for the buffer and the field
            const size_t buf_idx = id * nmax;
            const size_t my_idx  = localIndex(ax0, 0, i1, i2, ax0, inmem, nf);

            // do the copy -> vectorized
            for (size_t i0 = 0; i0 < nmax; i0++) {
                data[buf_idx + i0] = my_v[my_idx + i0];
            }
        }
    }
    if (_prof != NULL) {
        _prof->stop(profName); //mem2buf
    }

    //-------------------------------------------------------------------------
    /** - Do the communication */
    //-------------------------------------------------------------------------
    if (_is_all2all) {
        if (_prof != NULL) {
            profName = to_string(iswitch) + "all_2_all";
            _prof->start(profName);
        }
        MPI_Alltoall(sendBuf[0], send_count[0], MPI_DOUBLE, recvBuf[0], recv_count[0], MPI_DOUBLE, _subcomm);
        if (_prof != NULL) {
            _prof->stop(profName);//all_2_all
            int loc_mem = send_count[0] * comm_size;
            _prof->addMem(profName, loc_mem*sizeof(double));
        }

    } else {
        if (_prof != NULL) {
            profName = to_string(iswitch) + "all_2_all_v";
            _prof->start(profName);
        }
        MPI_Alltoallv(sendBuf[0], send_count, send_start, MPI_DOUBLE, recvBuf[0], recv_count, recv_start, MPI_DOUBLE, _subcomm);
        if (_prof != NULL) {
            _prof->stop(profName); //all_2_all_v
            int loc_mem = 0;
            for (int ir = 0; ir < comm_size; ir++) {
                loc_mem += send_count[ir];
            }
            _prof->addMem(profName, loc_mem*sizeof(double));
        }
    }

    //-------------------------------------------------------------------------
    /** - reset the memory to 0 */
    //-------------------------------------------------------------------------
    // reset the memory to 0
    std::memset(v, 0, sizeof(double) * topo_out->memsize());

    //-------------------------------------------------------------------------
    /** - wait for a block and copy when it arrives */
    //-------------------------------------------------------------------------
    // get some counters
    const int nblocks_recv = recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2];
    const int out_axis     = topo_out->axis();
    // for each block
    if (_prof != NULL) {
        profName = to_string(iswitch) + "buf2mem";
        _prof->start(profName);
    }

#pragma omp parallel default(none) proc_bind(close) firstprivate(nblocks_recv, recv_nBlock, v, recvBuf, ostart, nByBlock, oBlockSize, nf, onmem, ax0, ax1, ax2)
    for (int bid = 0; bid < nblocks_recv; bid++) {
        // get the indexing of the block in 012-indexing
        int ibv[3];
        localSplit(bid, recv_nBlock, 0, ibv, 1);

        // get the starting index in the global memory using !!nByBlock!!
        // since only the last block may have a different size
        const int loci0 = ostart[ax0] + ibv[ax0] * nByBlock[ax0];
        const int loci1 = ostart[ax1] + ibv[ax1] * nByBlock[ax1];
        const int loci2 = ostart[ax2] + ibv[ax2] * nByBlock[ax2];
        // get the memory
        opt_double_ptr data = recvBuf[bid];
        double* __restrict my_v = v + localIndex(ax0, loci0, loci1, loci2, out_axis, onmem, nf);
        // get the stride
        const size_t stride = localIndex(ax0, 1, 0, 0, out_axis, onmem, nf);
        // get the max number of ids not aligned in ax0
        const size_t id_max = oBlockSize[ax1][bid] * oBlockSize[ax2][bid];

        if (nf == 1) {
#pragma omp for schedule(static)
            for (size_t id = 0; id < id_max; id++) {
                // get the id from a small modulo
                const int i2 = id / oBlockSize[ax1][bid];
                const int i1 = id % oBlockSize[ax1][bid];
                // get the starting global id for the buffer and the field
                const size_t buf_idx = id * oBlockSize[ax0][bid] * nf;
                const size_t my_idx  = localIndex(ax0, 0, i1, i2, out_axis, onmem, nf);
                // do the copy
#pragma vector always
                for (int i0 = 0; i0 < oBlockSize[ax0][bid]; i0++) {
                    my_v[my_idx + i0 * stride] = data[buf_idx + i0];
                }
            }
        } else if (nf == 2) {
#pragma omp for schedule(static)
            for (size_t id = 0; id < id_max; id++) {
                // get the id from a small modulo
                const int i2 = id / oBlockSize[ax1][bid];
                const int i1 = id % oBlockSize[ax1][bid];
                // get the starting global id for the buffer and the field
                const size_t buf_idx = id * oBlockSize[ax0][bid] * nf;
                const size_t my_idx  = localIndex(ax0, 0, i1, i2, out_axis, onmem, nf);
                // do the copy
                // This loop is the most painful ! no vectorization possible when stride>2, but substantial speedup if stride==; this is why we split.
		const int nmax = oBlockSize[ax0][bid];
		if(stride==2){
#pragma vector always
                	for (int i0 = 0; i0 < 2*nmax; i0++) {
                	    my_v[my_idx + i0] = data[buf_idx + i0];
			    //memcpy(&my_v[my_idx+i0*2],&data[buf_idx+i0*2],2*sizeof(double));
			}
		}else{
#pragma ivdep
                	for (int i0 = 0; i0 < nmax; i0++) {
                    	    my_v[my_idx + i0 * stride + 0] = data[buf_idx + i0 * 2 + 0];
                    	    my_v[my_idx + i0 * stride + 1] = data[buf_idx + i0 * 2 + 1];
                            // memcpy(&my_v[my_idx+i0*stride],&data[buf_idx+i0*2],2*sizeof(double));
                        }
		}
            }
        }
    }

    if (_prof != NULL) {
        _prof->stop(profName);//buf2mem
        profName = "switch"+to_string(iswitch);
        _prof->stop(profName);
        _prof->stop("reorder");
    }
    END_FUNC;
}

void SwitchTopo_a2a::disp() const {
    BEGIN_FUNC;
    FLUPS_INFO("------------------------------------------");
    if (_is_all2all) FLUPS_INFO("## Topo Swticher All to All !! MPI");
    if (!_is_all2all) FLUPS_INFO("## Topo Swticher All to All vector MPI");
    FLUPS_INFO("--- INPUT");
    FLUPS_INFO("  - input axis = %d", _topo_in->axis());
    FLUPS_INFO("  - input local = %d %d %d", _topo_in->nloc(0), _topo_in->nloc(1), _topo_in->nloc(2));
    FLUPS_INFO("  - input global = %d %d %d", _topo_in->nglob(0), _topo_in->nglob(1), _topo_in->nglob(2));
    FLUPS_INFO("  - istart = %d %d %d", _istart[0], _istart[1], _istart[2]);
    FLUPS_INFO("  - iend = %d %d %d", _iend[0], _iend[1], _iend[2]);
    FLUPS_INFO("--- OUTPUT");
    FLUPS_INFO("  - output axis = %d", _topo_out->axis());
    FLUPS_INFO("  - output local = %d %d %d", _topo_out->nloc(0), _topo_out->nloc(1), _topo_out->nloc(2));
    FLUPS_INFO("  - output global = %d %d %d", _topo_out->nglob(0), _topo_out->nglob(1), _topo_out->nglob(2));
    FLUPS_INFO("  - ostart = %d %d %d", _ostart[0], _ostart[1], _ostart[2]);
    FLUPS_INFO("  - oend = %d %d %d", _oend[0], _oend[1], _oend[2]);
    FLUPS_INFO("--- BLOCKS");
    FLUPS_INFO("  - nByBlock  = %d %d %d", _nByBlock[0], _nByBlock[1], _nByBlock[2]);
    FLUPS_INFO("  - inBlock = %d %d %d", _inBlock[0], _inBlock[1], _inBlock[2]);
    FLUPS_INFO("  - onBlock = %d %d %d", _onBlock[0], _onBlock[1], _onBlock[2]);
    FLUPS_INFO("------------------------------------------");
}

void SwitchTopo_a2a_test() {
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {8, 8, 8};
    const int nproc[3] = {2, 2, 1};

    const int nglob_big[3] = {17, 8, 8};
    const int nproc_big[3] = {2, 2, 1};

    //===========================================================================
    // real numbers
    Topology* topo    = new Topology(2, nglob, nproc, false, NULL,1);
    Topology* topobig = new Topology(0, nglob_big, nproc_big, false, NULL,1);

    double* data = (double*)fftw_malloc(sizeof(double*) * std::max(topo->memsize(), topobig->memsize()));

    const int nmem[3] = {topo->nmem(0),topo->nmem(1),topo->nmem(2)};
    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                const size_t id    = localIndex(0, i0, i1, i2, 0, nmem, 2);
                
                data[id + 0] = id;
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_real", data);

    const int fieldstart[3] = {0, 0, 0};
    // printf("\n=============================");
    SwitchTopo* switchtopo = new SwitchTopo_a2a(topo, topobig, fieldstart, NULL);
    size_t          max_mem    = switchtopo->get_bufMemSize();
    opt_double_ptr  send_buff  = (opt_double_ptr)fftw_malloc(max_mem * sizeof(double));
    opt_double_ptr  recv_buff  = (opt_double_ptr)fftw_malloc(max_mem * sizeof(double));
    std::memset(send_buff, 0, max_mem * sizeof(double));
    std::memset(recv_buff, 0, max_mem * sizeof(double));
    // associate the buffer
    switchtopo->setup_buffers(send_buff, recv_buff);

    switchtopo->disp();

    // printf("\n\n============ FORWARD =================");
    switchtopo->execute(data, FLUPS_FORWARD, 0);

    hdf5_dump(topobig, "test_real_padd", data);

    // printf("\n\n============ BACKWARD =================");
    switchtopo->execute(data, FLUPS_BACKWARD, 0);

    hdf5_dump(topo, "test_real_returned", data);

    fftw_free(data);
    fftw_free(send_buff);
    fftw_free(recv_buff);
    delete (switchtopo);
    delete (topo);
    delete (topobig);

    // //===========================================================================
    // // complex numbers
    // topo    = new Topology(0, nglob, nproc, true,NULL,1);
    // topobig = new Topology(2, nglob_big, nproc_big, true,NULL,1);

    // data = (double*)fftw_malloc(sizeof(double*) * topobig->memsize());

    // const int ax0 = topo->axis();
    // const int nmem[3] = {topo->nmem(0),topo->nmem(1),topo->nmem(2)};
    // for (int i2 = 0; i2 < topo->nloc(2); i2++) {
    //     for (int i1 = 0; i1 < topo->nloc(1); i1++) {
    //         for (int i0 = 0; i0 < topo->nloc(0); i0++) {
    //             size_t id    = localIndex(0,i0, i1, i2,0,nmem,2);
    //             data[id + 0] = 0;
    //             data[id + 1] = id;
    //         }
    //     }
    // }
    // // try the dump
    // hdf5_dump(topo, "test_complex", data);

    // // topobig->switch2complex();
    // // printf("as complex: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));
    // // topobig->switch2real();
    // // printf("as real: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));

    // const int fieldstart2[3] = {4, 0, 0};
    // // printf("\n=============================");
    // switchtopo = new SwitchTopo_a2a(topo, topobig, fieldstart2, NULL);

    // switchtopo->execute(data, FLUPS_FORWARD);

    // hdf5_dump(topobig, "test_complex_padd", data);

    // switchtopo->execute(data, FLUPS_BACKWARD);

    // hdf5_dump(topo, "test_complex_returned", data);

    // fftw_free(data);
    // delete (switchtopo);
    // delete (topo);
    // delete (topobig);
    END_FUNC;
}
