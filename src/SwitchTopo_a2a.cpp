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
    int* nperProc = (int*)fftw_malloc(comm_size * sizeof(int));
    for (int id = 0; id < 3; id++) {
        // get the gcd between send and receive
        int isend = (_iend[id] - _istart[id]);
        int osend = (_oend[id] - _ostart[id]);

        // compute the exchanged size same if from the input or output
        MPI_Allreduce(&isend, &_exSize[id], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        // we have summed the size nproc(id+1)*size nproc(id+2) * size, so we divide
        _exSize[id] /= _topo_in->nproc((id+1)%3) * _topo_in->nproc((id+2)%3);

        // if I am the last one, I decrease the blocksize by one if needed
        if (_topo_in->rankd(id) == (_topo_in->nproc(id) - 1)) {
            isend = isend - _exSize[id] % 2;
        }
        if (_topo_out->rankd(id) == (_topo_out->nproc(id) - 1)) {
            osend = osend - _exSize[id] % 2;
        }
        int npoints = gcd(isend, osend);
        // gather on each proc the gcd
        MPI_Allgather(&npoints, 1, MPI_INT, nperProc, 1, MPI_INT, MPI_COMM_WORLD);
        // get the Greatest Common Divider among every process
        int my_gcd = nperProc[0];
        for (int ip = 1; ip < comm_size; ip++) {
            my_gcd = gcd(my_gcd, nperProc[ip]);
        }
        // store it as the block dimension
        _nByBlock[id] = my_gcd;
    }
    fftw_free(nperProc);

    //-------------------------------------------------------------------------
    /** - get the number of blocks and for each block get the size and the destination rank */
    //-------------------------------------------------------------------------
    int  iblockIDStart[3];
    int  oblockIDStart[3];
    int* inBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));
    int* onBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));

    cmpt_blockIndexes(_istart, _iend, _nByBlock, _topo_in, _inBlock, iblockIDStart, inBlockEachProc);
    cmpt_blockIndexes(_ostart, _oend, _nByBlock, _topo_out, _onBlock, oblockIDStart, onBlockEachProc);

    // allocte the block size
    for (int id = 0; id < 3; id++) {
        _iBlockSize[id] = (int*)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
        _oBlockSize[id] = (int*)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));
    }

    // allocate the destination ranks
    _i2o_destRank = (opt_int_ptr)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
    _o2i_destRank = (opt_int_ptr)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));

    // get the send destination ranks in the ouput topo
    cmpt_blockSize(_inBlock, iblockIDStart, _nByBlock, _istart, _iend, _iBlockSize);
    cmpt_blockSize(_onBlock, oblockIDStart, _nByBlock, _ostart, _oend, _oBlockSize);

    cmpt_blockDestRank(_inBlock, iblockIDStart, _topo_out, onBlockEachProc, _i2o_destRank);
    cmpt_blockDestRank(_onBlock, oblockIDStart, _topo_in, inBlockEachProc, _o2i_destRank);

    // free the temp arrays
    fftw_free(inBlockEachProc);
    fftw_free(onBlockEachProc);

    //-------------------------------------------------------------------------
    /** - do the communication split */
    //-------------------------------------------------------------------------
    // compute the color among the proc I send to and I recv from
    FLUPS_INFO("Trying to determine the MPI communicators...");
    int   mycolor   = rank;
    int*  colors    = (int*)fftw_malloc(comm_size * sizeof(int));
    bool* inMyGroup = (bool*)fftw_malloc(comm_size * sizeof(bool));

    for (int ir = 0; ir < comm_size; ir++) {
        inMyGroup[ir] = false;
    }

    // do a first pass and give a color + who is in my group
    for (int ib = 0; ib < _inBlock[0] * _inBlock[1] * _inBlock[2]; ib++) {
        mycolor                      = std::min(mycolor, _i2o_destRank[ib]);
        inMyGroup[_i2o_destRank[ib]] = true;
    }
    for (int ib = 0; ib < _onBlock[0] * _onBlock[1] * _onBlock[2]; ib++) {
        mycolor                      = std::min(mycolor, _o2i_destRank[ib]);
        inMyGroup[_o2i_destRank[ib]] = true;
    }

    // count how much is should be in my group
    // by default we assume that nobody is in the same group
    int nleft = 0;
    for (int ir = 0; ir < comm_size; ir++) {
        if (inMyGroup[ir]) {
            nleft += 1;
        }
    }
    // continue while we haven't found a solution
    while (nleft > 0) {
        // gather the color info from everyone
        MPI_Allgather(&mycolor, 1, MPI_INT, colors, 1, MPI_INT, MPI_COMM_WORLD);
        // iterate on the proc
        int n_notInMyGroup = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            // if it is reachable and the color is not already the same
            if (inMyGroup[ir] && (colors[ir] != mycolor)) {
                // we first increment the counter flagging that one is missing
                n_notInMyGroup += 1;
                // then we solve the problem if we are able to do so....
                // remove 1 if we are able to solve the issue <=> my color > colors[ir]
                n_notInMyGroup = n_notInMyGroup - (colors[ir] < mycolor);
                // changing if possible
                mycolor = std::min(mycolor, colors[ir]);
            }
        }
        // compute among everybody, if we need to continue
        MPI_Allreduce(&n_notInMyGroup, &nleft, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        FLUPS_INFO("stil %d to find (@ my proc: %d)", nleft, n_notInMyGroup);
    }
    fftw_free(colors);
    fftw_free(inMyGroup);

    FLUPS_INFO("Group found: my color = %d", mycolor);

    //-------------------------------------------------------------------------
    /** - create the new communicator and assocatied data */
    //-------------------------------------------------------------------------
    MPI_Comm_split(MPI_COMM_WORLD, mycolor, rank, &_subcomm);
    // set a name
    std::string commname = "comm-" + std::to_string(mycolor);
    MPI_Comm_set_name(_subcomm, commname.c_str());
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
        _prof->create("mem2buf", "reorder");
        _prof->create("buf2mem", "reorder");
        _prof->create("all_2_all", "reorder");
        _prof->create("all_2_all_v", "reorder");
    }

    //-------------------------------------------------------------------------
    /** - Display performance information if asked */
    //-------------------------------------------------------------------------
#ifdef PERF_VERBOSE
    // we display important information for the performance
    string name = "./prof/SwitchTopo_" + std::to_string(_topo_in->axis()) + "to" + std::to_string(_topo_out->axis()) + "_rank" + std::to_string(rank) + ".txt";
    FILE* file = fopen(name.c_str(),"w+");
    if(file != NULL){
        if(_is_all2all) fprintf(file,"- is all to all\n");
        if(!_is_all2all) fprintf(file,"- is vectorized all to all\n");
        
        int newrank;
        MPI_Comm_rank(_subcomm,&newrank);
        int rlen;
        char myname[MPI_MAX_OBJECT_NAME];
        MPI_Comm_get_name(_subcomm, myname, &rlen);
        fprintf(file,"- in subcom %s with rank %d/%d\n",myname,newrank,subsize);
        fprintf(file,"- nglob = %d %d %d to %d %d %d\n",_topo_in->nglob(0),_topo_in->nglob(1),_topo_in->nglob(2),_topo_out->nglob(0),_topo_out->nglob(1),_topo_out->nglob(2));
        fprintf(file,"- nproc = %d %d %d to %d %d %d\n",_topo_in->nproc(0),_topo_in->nproc(1),_topo_in->nproc(2),_topo_out->nproc(0),_topo_out->nproc(1),_topo_out->nproc(2));
        fprintf(file,"- nByBlock = %d %d %d, real size = %d %d %d\n",_nByBlock[0],_nByBlock[1],_nByBlock[2],_nByBlock[0]+_exSize[0]%2,_nByBlock[1]+_exSize[1]%2,_nByBlock[2]+_exSize[2]%2);

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
}

/**
 * @brief setup the comm lists
 * 
 * We setup the following lists:
 * - destRank: transformed from the values in the world comm to the values in the new comm.
 * - count: the number of elements send to each proc form this proc
 * - start: the starting position of the data to send to each proc in the buffer
 * 
 * @param newcomm the new comm
 * @param nBlock the number of blocks
 * @param destRank the destination rank on the world comm. returns the new destination rank in the newcomm
 * @param count the number of information to send to each proc
 * @param start the id in the buffer where the information starts for each proc
 */
void SwitchTopo_a2a::_setup_subComm(MPI_Comm newcomm, const int nBlock[3], int* destRank, int** count, int** start) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the new destination ranks    */
    //-------------------------------------------------------------------------
    int newrank, worldsize;
    MPI_Comm_rank(newcomm, &newrank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);

    // get the new ranks from my old friends in the old communicator
    int* newRanks = (int*)fftw_malloc(worldsize * sizeof(int));
    MPI_Allgather(&newrank, 1, MPI_INT, newRanks, 1, MPI_INT, MPI_COMM_WORLD);
    // replace the old ranks by the newest ones
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        destRank[ib] = newRanks[destRank[ib]];
    }
    fftw_free(newRanks);

    //-------------------------------------------------------------------------
    /** - build the size vector of block to each procs    */
    //-------------------------------------------------------------------------
    int subsize;
    MPI_Comm_size(_subcomm, &subsize);
    // count the number of blocks to each process
    (*count) = (int*)fftw_malloc(subsize * sizeof(int));
    (*start) = (int*)fftw_malloc(subsize * sizeof(int));
    std::memset((*count), 0, subsize * sizeof(int));
    std::memset((*start), 0, subsize * sizeof(int));
    // get the size per block
    const int blockMem = get_blockMemSize();
    // count the number of blocks per rank
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        (*count)[destRank[ib]] += blockMem;
    }
    // compute the start indexes
    (*start)[0] = 0;
    for (int ir = 1; ir < subsize; ir++) {
        (*start)[ir] = (*start)[ir - 1] + (*count)[ir - 1];
    }
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
void SwitchTopo_a2a::execute(opt_double_ptr v, const int sign) {
    BEGIN_FUNC;

    FLUPS_CHECK(_topo_in->isComplex() == _topo_out->isComplex(), "both topologies have to be complex or real", LOCATION);
    FLUPS_CHECK(_topo_in->nf() <= 2, "the value of nf is not supported", LOCATION);

    int comm_size;
    // MPI_Comm_rank(_subcomm, &rank);
    MPI_Comm_size(_subcomm, &comm_size);

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
    int inloc[3];
    int onloc[3];

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
            inloc[id]       = _topo_in->nloc(id);
            onloc[id]       = _topo_out->nloc(id);
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
            inloc[id]       = _topo_out->nloc(id);
            onloc[id]       = _topo_in->nloc(id);
            iBlockSize[id]  = _oBlockSize[id];
            oBlockSize[id]  = _iBlockSize[id];
        }
    } else {
        FLUPS_CHECK(false, "the sign is not FLUPS_FORWARD nor FLUPS_BACKWARD", LOCATION);
    }

    FLUPS_INFO("previous topo: %d,%d,%d axis=%d", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    FLUPS_INFO("new topo: %d,%d,%d  axis=%d", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());
    FLUPS_INFO("using %d blocks on send and %d on recv", send_nBlock[0] * send_nBlock[1] * send_nBlock[2], recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2]);

    // define important constants
    const int ax0 = topo_in->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nf  = topo_in->nf();

    if (_prof != NULL) {
        _prof->start("mem2buf");
    }
    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    const int nblocks_send = send_nBlock[0] * send_nBlock[1] * send_nBlock[2];

#pragma omp parallel proc_bind(close) default(none) firstprivate(nblocks_send, send_nBlock, v, sendBuf, istart, nByBlock, iBlockSize, nf, inloc, ax0, ax1, ax2)
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
        double* __restrict data = sendBuf[bid];
        double* __restrict my_v = v + localIndex(ax0, loci0, loci1, loci2, ax0, inloc, nf);

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
            const size_t my_idx  = localIndex(ax0, 0, i1, i2, ax0, inloc, nf);

            // do the copy -> vectorized
            for (size_t i0 = 0; i0 < nmax; i0++) {
                data[buf_idx + i0] = my_v[my_idx + i0];
            }
        }
    }
    if (_prof != NULL) {
        _prof->stop("mem2buf");
    }

    //-------------------------------------------------------------------------
    /** - Do the communication */
    //-------------------------------------------------------------------------
    if (_is_all2all) {
        if (_prof != NULL) {
            _prof->start("all_2_all");
        }
        MPI_Alltoall(sendBuf[0], send_count[0], MPI_DOUBLE, recvBuf[0], recv_count[0], MPI_DOUBLE, _subcomm);
        if (_prof != NULL) {
            _prof->stop("all_2_all");
            int loc_mem = send_count[0] * comm_size;
            _prof->addMem("all_2_all", loc_mem);
        }

    } else {
        if (_prof != NULL) {
            _prof->start("all_2_all_v");
        }
        MPI_Alltoallv(sendBuf[0], send_count, send_start, MPI_DOUBLE, recvBuf[0], recv_count, recv_start, MPI_DOUBLE, _subcomm);
        if (_prof != NULL) {
            _prof->stop("all_2_all_v");
            int loc_mem = 0;
            for (int ir = 0; ir < comm_size; ir++) {
                loc_mem += send_count[ir];
            }
            _prof->addMem("all_2_all_v", loc_mem);
        }
    }

    //-------------------------------------------------------------------------
    /** - reset the memory to 0 */
    //-------------------------------------------------------------------------
    // reset the memory to 0
    std::memset(v, 0, sizeof(double) * topo_out->locmemsize());

    //-------------------------------------------------------------------------
    /** - wait for a block and copy when it arrives */
    //-------------------------------------------------------------------------
    // get some counters
    const int nblocks_recv = recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2];
    const int out_axis     = topo_out->axis();
    // for each block
    if (_prof != NULL) {
        _prof->start("buf2mem");
    }

#pragma omp parallel default(none) proc_bind(close) firstprivate(nblocks_recv, recv_nBlock, v, recvBuf, ostart, nByBlock, oBlockSize, nf, onloc, ax0, ax1, ax2, out_axis)
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
        double* __restrict data = recvBuf[bid];
        double* __restrict my_v = v + localIndex(ax0, loci0, loci1, loci2, out_axis, onloc, nf);
        // get the stride
        const size_t stride = localIndex(ax0, 1, 0, 0, out_axis, onloc, nf);
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
                const size_t my_idx  = localIndex(ax0, 0, i1, i2, out_axis, onloc, nf);
                // do the copy
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
                const size_t my_idx  = localIndex(ax0, 0, i1, i2, out_axis, onloc, nf);
                // do the copy
                for (int i0 = 0; i0 < oBlockSize[ax0][bid]; i0++) {
                    my_v[my_idx + i0 * stride + 0] = data[buf_idx + i0 * 2 + 0];
                    my_v[my_idx + i0 * stride + 1] = data[buf_idx + i0 * 2 + 1];
                }
            }
        }
    }

    if (_prof != NULL) {
        _prof->stop("buf2mem");
        _prof->stop("reorder");
    }
}

void SwitchTopo_a2a::disp() {
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
    Topology* topo    = new Topology(2, nglob, nproc, false, NULL);
    Topology* topobig = new Topology(0, nglob_big, nproc_big, false, NULL);

    double* data = (double*)fftw_malloc(sizeof(double*) * std::max(topo->locmemsize(), topobig->locmemsize()));

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localindex_xyz(i0, i1, i2, topo);
                data[id + 0] = id;
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_real", data);

    const int fieldstart[3] = {0, 0, 0};
    // printf("\n=============================");
    SwitchTopo_a2a* switchtopo = new SwitchTopo_a2a(topo, topobig, fieldstart, NULL);
    size_t          max_mem    = switchtopo->get_bufMemSize();
    opt_double_ptr  send_buff  = (opt_double_ptr)fftw_malloc(max_mem * sizeof(double));
    opt_double_ptr  recv_buff  = (opt_double_ptr)fftw_malloc(max_mem * sizeof(double));
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

    fftw_free(data);
    fftw_free(send_buff);
    fftw_free(recv_buff);
    delete (switchtopo);
    delete (topo);
    delete (topobig);

    // //===========================================================================
    // // complex numbers
    // topo    = new Topology(0, nglob, nproc, true,NULL);
    // topobig = new Topology(2, nglob_big, nproc_big, true,NULL);

    // data = (double*)fftw_malloc(sizeof(double*) * topobig->locmemsize());

    // for (int i2 = 0; i2 < topo->nloc(2); i2++) {
    //     for (int i1 = 0; i1 < topo->nloc(1); i1++) {
    //         for (int i0 = 0; i0 < topo->nloc(0); i0++) {
    //             size_t id    = localindex_xyz(i0, i1, i2, topo);
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
}