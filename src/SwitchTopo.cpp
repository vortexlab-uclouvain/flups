/**
 * @file SwitchTopo.cpp
 * @author Thomas Gillis
 * @brief
 * @version
 * @date 2019-07-25
 *
 * @copyright Copyright © UCLouvain 2019
 *
 */

#include "SwitchTopo.hpp"

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
 * 
 * @param topo_input the input topology
 * @param topo_output the output topology 
 * @param shift the shift is the position of the (0,0,0) of topo_input in the topo_output indexing (in XYZ-indexing)
 */
SwitchTopo::SwitchTopo(const Topology* topo_input, const Topology* topo_output,
                       const int shift[3]) {
    BEGIN_FUNC

    UP_CHECK0(topo_input->isComplex() == topo_output->isComplex(),
              "both topologies have to be the same kind");

    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    _topo_in  = topo_input;
    _topo_out = topo_output;

    //-------------------------------------------------------------------------
    /** - allocate the size arrays */
    //-------------------------------------------------------------------------
    _nsend = (int*)fftw_malloc(comm_size * sizeof(int));
    _nrecv = (int*)fftw_malloc(comm_size * sizeof(int));
    _ssend = (int*)fftw_malloc(comm_size * sizeof(int));
    _srecv = (int*)fftw_malloc(comm_size * sizeof(int));
    _count = (int*)fftw_malloc(comm_size * sizeof(int));
    // n axis arrays
    _naxis_i2o_send = (int*)fftw_malloc(comm_size * sizeof(int));
    _naxis_o2i_send = (int*)fftw_malloc(comm_size * sizeof(int));
    _naxis_i2o_recv = (int*)fftw_malloc(comm_size * sizeof(int));
    _naxis_o2i_recv = (int*)fftw_malloc(comm_size * sizeof(int));

    UP_CHECK1(UP_ISALIGNED(_nsend), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_nrecv), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_ssend), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_srecv), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_count), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_naxis_i2o_send), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_naxis_o2i_send), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_naxis_i2o_recv), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(_naxis_o2i_recv), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);

    std::memset(_nsend, 0, sizeof(int) * comm_size);
    std::memset(_nrecv, 0, sizeof(int) * comm_size);
    std::memset(_ssend, 0, sizeof(int) * comm_size);
    std::memset(_srecv, 0, sizeof(int) * comm_size);
    std::memset(_count, 0, sizeof(int) * comm_size);

    //-------------------------------------------------------------------------
    /** - Compute the naxis for loop go through */
    //-------------------------------------------------------------------------
    // const int ax0 = _topo_in->axis();
    // const int ax1 = (ax0 + 1) % 3;
    // const int ax2 = (ax0 + 2) % 3;

    //-------------------------------------------------------------------------
    /** - for each dimension, get the shared zone */
    //-------------------------------------------------------------------------
    // for (int id = 0; id < 3; id++) {
    //     _ishift[id] = shift[id];
    //     _oshift[id] = -shift[id];
    // }

    // _topo_in->cmpt_intersect_id(_ishift, _topo_out, istart, _iend);
    // _topo_out->cmpt_intersect_id(_oshift, _topo_in, ostart, _oend);

    //-------------------------------------------------------------------------
    /** - get the starting and ending index of the shared zone */
    //-------------------------------------------------------------------------
    // get the blockshift
    for (int id = 0; id < 3; id++) {
        _ib2o_shift[id] = shift[id];
        _ob2i_shift[id] = -shift[id];
    }
    // get how much pints we send/recv in each direction
    _topo_in->cmpt_intersect_id(_ib2o_shift, _topo_out, _istart, _iend);
    _topo_out->cmpt_intersect_id(_ob2i_shift, _topo_in, _ostart, _oend);

    
    
    //-------------------------------------------------------------------------
    /** - get the block size as the GCD among every process between send and receive*/
    //-------------------------------------------------------------------------
    int* onProc    = (int*)fftw_malloc(comm_size * sizeof(int));
    for(int id=0; id<3; id++){
        // get the gcd between send and receive
        int  npoints = gcd(_iend[id] - _istart[id], _oend[id] - _ostart[id]);
        // printf("npoints = %d from in: %d %d - out: %d %d\n",npoints,_istart[id],_iend[id],_ostart[id],_oend[id]);
        // gather on each proc the gcd
        MPI_Allgather(&npoints, 1, MPI_INT, onProc, 1, MPI_INT, MPI_COMM_WORLD);
        // get the Greatest Common Divider among every process
        int my_gcd=onProc[0];
        for(int ip=1; ip<comm_size; ip++){
            my_gcd = gcd(my_gcd,onProc[ip]);
        }
        // store it as the block dimension
        _nByBlock[id] = my_gcd;

        // // store the number of blocks
        // _nBlock[id] = std::min(_topo_in->nglob(id),_topo_out->nglob(id))/_nByBlock[id];
    }
    fftw_free(onProc);
    
    // printf("send topo: nloc = %d %d %d \n",_topo_in->nloc(0),_topo_in->nloc(1),_topo_in->nloc(2));
    // printf("recv topo: nloc = %d %d %d \n",_topo_out->nloc(0),_topo_out->nloc(1),_topo_out->nloc(2));
    // printf("\nnByBlock = %d %d %d\n",_nByBlock[0],_nByBlock[1],_nByBlock[2]);

    //-------------------------------------------------------------------------
    /** - get the starting index of the block 0,0,0 for input and output */
    //-------------------------------------------------------------------------
    int* inBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));
    int* onBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));

    cmpt_blockIndexes(_istart,_iend,_nByBlock,_topo_in,_inBlock,_iblockIDStart,inBlockEachProc);
    cmpt_blockIndexes(_ostart,_oend,_nByBlock,_topo_out,_onBlock,_oblockIDStart,onBlockEachProc);
    

    // printf("the proc has %d %d %d blocks in INPUT conf\n",_inBlock[0],_inBlock[1],_inBlock[2]);
    // printf("the proc has %d %d %d blocks in OUTPUT conf\n",_onBlock[0],_onBlock[1],_onBlock[2]);

    //-------------------------------------------------------------------------
    /** - allocate the arrays */
    //-------------------------------------------------------------------------
    // allocate the buffers
    _sendBuf = (opt_double_ptr*)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(double*));
    _recvBuf = (opt_double_ptr*)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(double*));
    // allocate the requests
    _sendRequest = (MPI_Request*)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(MPI_Request));
    _recvRequest = (MPI_Request*)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(MPI_Request));
    // allocate the destination ranks
    _i2o_destRank = (int*)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
    _o2i_destRank = (int*)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));
    // allocate the destination tags
    _i2o_destTag = (int*)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
    _o2i_destTag = (int*)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));

    //-------------------------------------------------------------------------
    /** - for each block, get the destination rank */
    //-------------------------------------------------------------------------
    // send destination ranks in the ouput topo
    cmpt_blockDestRankAndTag(_inBlock,_iblockIDStart,_topo_out,onBlockEachProc,_i2o_destRank,_i2o_destTag);
    cmpt_blockDestRankAndTag(_onBlock,_oblockIDStart,_topo_in, inBlockEachProc,_o2i_destRank,_o2i_destTag);

    // free the temp arrays
    fftw_free(inBlockEachProc);
    fftw_free(onBlockEachProc);
    
    //-------------------------------------------------------------------------
    /** - for each block allocate the data buffer */
    //-------------------------------------------------------------------------
    for (int ib2 = 0; ib2 < _inBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < _inBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < _inBlock[0]; ib0++) {
                
                // get the local block index
                const size_t send_bid = blockID(ib0, ib1, ib2, _inBlock);
                // printf("sending block %d to proc %d with tag %d\n",send_bid,_i2o_destRank[send_bid],_i2o_destTag[send_bid]);
                // allocate at the correct size
                _sendBuf[send_bid] = (double*) fftw_malloc(_nByBlock[0] * _nByBlock[1] * _nByBlock[2] * sizeof(double) * _topo_in->nf());
                UP_CHECK1(UP_ISALIGNED(_sendBuf[send_bid]), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
            }
        }
    }
    for (int ib2 = 0; ib2 < _onBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < _onBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < _onBlock[0]; ib0++) {
                // get the local block index
                const size_t recv_bid = blockID(ib0, ib1, ib2, _onBlock);
                // printf("recving block %d from proc %d with tag %d\n",recv_bid,_o2i_destRank[recv_bid],_o2i_destTag[recv_bid]);
                // allocate at the correct size
                _recvBuf[recv_bid] = (double*) fftw_malloc(_nByBlock[0] * _nByBlock[1] * _nByBlock[2] * sizeof(double) * _topo_in->nf());
                UP_CHECK1(UP_ISALIGNED(_recvBuf[recv_bid]), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
            }
        }
    }

                

    // //-------------------------------------------------------------------------
    // /** - allocate the arrays */
    // //-------------------------------------------------------------------------
    // // allocate the buffers
    // _sendbuf        =(double**) fftw_malloc(_send_nBlockByProc[0] * _send_nBlockByProc[1] * _send_nBlockByProc[2] * sizeof(double*));
    // _recvbuf        =(double**) fftw_malloc(_recv_nBlockByProc[0] * _recv_nBlockByProc[1] * _recv_nBlockByProc[2] * sizeof(double*));
    // _send_blockSize =(bcoord*) fftw_malloc(_send_nBlockByProc[0] * _send_nBlockByProc[1] * _send_nBlockByProc[2] * sizeof(bcoord));
    // _recv_blockSize =(bcoord*) fftw_malloc(_recv_nBlockByProc[0] * _recv_nBlockByProc[1] * _recv_nBlockByProc[2] * sizeof(bcoord));
    // _send_request   =(MPI_Request*) fftw_malloc(_send_nBlockByProc[0] * _send_nBlockByProc[1] * _send_nBlockByProc[2] * sizeof(MPI_Request));
    // _recv_request   =(MPI_Request*) fftw_malloc(_recv_nBlockByProc[0] * _recv_nBlockByProc[1] * _recv_nBlockByProc[2] * sizeof(MPI_Request));
    // UP_CHECK1(UP_ISALIGNED(_send_blockSize), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    // UP_CHECK1(UP_ISALIGNED(_sendbuf), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    // UP_CHECK1(UP_ISALIGNED(_recv_blockSize), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    // UP_CHECK1(UP_ISALIGNED(_recvbuf), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    // UP_CHECK1(UP_ISALIGNED(_send_request), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    // UP_CHECK1(UP_ISALIGNED(_recv_request), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);

    // // int send_bid_start[3] = {topo_in->rankd(0) * _nBlockByProc[0], topo_in->rankd(1) * _nBlockByProc[1], topo_in->rankd(2) * _nBlockByProc[2]};
    // // int recv_bid_start[3] = {topo_out->rankd(0) * _nBlockByProc[0], topo_out->rankd(1) * _nBlockByProc[1], topo_out->rankd(2) * _nBlockByProc[2]};
    // // ---------- SEND -------------------
    // for (int ib2 = 0; ib2 < _send_nBlockByProc[2]; ib2++) {
    //     for (int ib1 = 0; ib1 < _send_nBlockByProc[1]; ib1++) {
    //         for (int ib0 = 0; ib0 < _send_nBlockByProc[0]; ib0++) {
    //             // get the local block index
    //             const int    send_bid[3] = {ib0, ib1, ib2};
    //             const size_t send_bindex = blockID(ib0, ib1, ib2, _send_nBlockByProc);
    //             // get the size of the block in 3D
    //             cmpt_blockSize(_topo_in, _topo_out, _send_nBlockByProc, _nByBlock, send_bid, _send_blockSize[send_bindex]);
    //             // allocate the buffer for the block
    //             _sendbuf[send_bindex] = (double*) fftw_malloc(_send_blockSize[send_bindex][0] * _send_blockSize[send_bindex][1] * _send_blockSize[send_bindex][2] * sizeof(double) * _topo_in->nf());
    //             UP_CHECK1(UP_ISALIGNED(_sendbuf[send_bindex]), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    //             // get the rank
    //             const int destrank = cmpt_blockRank(_topo_in, _topo_out, _ib2o_shift, _send_nBlockByProc,_nByBlock, send_bid);
    //             // add one block
    //             _nsend[destrank] = _nsend[destrank] + 1;
    //         }
    //     }
    // }
    // // ---------- RECEIVE -------------------
    // for (int ib2 = 0; ib2 < _recv_nBlockByProc[2]; ib2++) {
    //     for (int ib1 = 0; ib1 < _recv_nBlockByProc[1]; ib1++) {
    //         for (int ib0 = 0; ib0 < _recv_nBlockByProc[0]; ib0++) {
                
    //             // get global block index
    //             const int    recv_bid[3] = {ib0, ib1, ib2};
    //             const size_t recv_bindex = blockID(ib0, ib1, ib2, _recv_nBlockByProc);
    //             // get the size of the block in 3D
    //             cmpt_blockSize(_topo_out, _topo_in, _recv_nBlockByProc, _nByBlock, recv_bid, _recv_blockSize[recv_bindex]);
    //             // allocate the buffer for the block
    //             _recvbuf[recv_bindex] = (double*) fftw_malloc(_recv_blockSize[recv_bindex][0] * _recv_blockSize[recv_bindex][1] * _recv_blockSize[recv_bindex][2] * sizeof(double) * _topo_in->nf());
    //             UP_CHECK1(UP_ISALIGNED(_recvbuf[recv_bindex]), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    //             // get the rank
    //             const int origrank = cmpt_blockRank(_topo_out, _topo_in, _ob2i_shift,_recv_nBlockByProc, _nByBlock, recv_bid);
    //             // add one block
    //             _nrecv[origrank] = _nrecv[origrank] + 1;
    //         }
    //     }
    // }
    // //-------------------------------------------------------------------------
    // /** - for each data get its destination rank */
    // //-------------------------------------------------------------------------
    // int dest_rankd[3];
    // for (int i2 = _istart[ax2]; i2 < _iend[ax2]; ++i2) {
    //     dest_rankd[ax2] = _topo_in->cmpt_matchrank(ax2, _topo_out, i2 + _ishift[ax2]);
    //     for (int i1 = _istart[ax1]; i1 < _iend[ax1]; ++i1) {
    //         dest_rankd[ax1] = _topo_in->cmpt_matchrank(ax1, _topo_out, i1 + _ishift[ax1]);
    //         for (int i0 = _istart[ax0]; i0 < _iend[ax0]; ++i0) {
    //             dest_rankd[ax0] = _topo_in->cmpt_matchrank(ax0, _topo_out, i0 + _ishift[ax0]);
    //             // add one info to the destination index
    //             _nsend[rankindex(dest_rankd, _topo_out)] += _topo_in->nf();
    //         }
    //     }
    // }

    //-------------------------------------------------------------------------
    /** - send the number of block to each proc and store it */
    //-------------------------------------------------------------------------
    // MPI_Alltoall(_nsend, 1, MPI_INT, _nrecv, 1, MPI_INT, MPI_COMM_WORLD);

    // for (int i = 1; i < comm_size; ++i) {
    //     _ssend[i] = _ssend[i - 1] + _nsend[i - 1];
    //     _srecv[i] = _srecv[i - 1] + _nrecv[i - 1];
    // }
    // MPI_Alltoall(_naxis_i2o_send, 1, MPI_INT, _naxis_i2o_recv, 1, MPI_INT, MPI_COMM_WORLD);
    // MPI_Alltoall(_naxis_o2i_send, 1, MPI_INT, _naxis_o2i_recv, 1, MPI_INT, MPI_COMM_WORLD);

    //-------------------------------------------------------------------------
    /** - allocate the buffer if needed */
    //-------------------------------------------------------------------------
    // _bufsend = (double*)fftw_malloc(sizeof(double) * _topo_in->locmemsize());
    // _bufrecv = (double*)fftw_malloc(sizeof(double) * _topo_out->locmemsize());

    // UP_CHECK1(UP_ISALIGNED(_bufsend), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    // UP_CHECK1(UP_ISALIGNED(_bufrecv), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
}

SwitchTopo::~SwitchTopo() {
    if (_nsend != NULL) fftw_free(_nsend);
    if (_nrecv != NULL) fftw_free(_nrecv);
    if (_ssend != NULL) fftw_free(_ssend);
    if (_srecv != NULL) fftw_free(_srecv);
    if (_srecv != NULL) fftw_free(_count);
    if (_bufsend != NULL) fftw_free(_bufsend);
    if (_bufrecv != NULL) fftw_free(_bufrecv);

    if (_naxis_i2o_send != NULL) fftw_free(_naxis_i2o_send);
    if (_naxis_o2i_send != NULL) fftw_free(_naxis_o2i_send);
    if (_naxis_i2o_recv != NULL) fftw_free(_naxis_i2o_recv);
    if (_naxis_o2i_recv != NULL) fftw_free(_naxis_o2i_recv);

    // if (_send_blockSize != NULL) fftw_free(_send_blockSize);
    // if (_sendbuf != NULL) fftw_free(_sendbuf);
    // if (_recv_blockSize != NULL) fftw_free(_recv_blockSize);
    // if (_recvbuf != NULL) fftw_free(_recvbuf);

    // for (int ib2 = 0; ib2 < _send_nBlockByProc[2]; ib2++) {
    //     for (int ib1 = 0; ib1 < _send_nBlockByProc[1]; ib1++) {
    //         for (int ib0 = 0; ib0 < _send_nBlockByProc[0]; ib0++) {
    //             // send
    //             const int    send_bid[3] = {ib0 , ib1 , ib2 };
    //             const size_t send_bindex = blockID(ib0, ib1, ib2, _send_nBlockByProc);
    //             if (_sendbuf[send_bindex] != NULL) fftw_free(_sendbuf[send_bindex]);
    //         }
    //     }
    // }
    // for (int ib2 = 0; ib2 < _recv_nBlockByProc[2]; ib2++) {
    //     for (int ib1 = 0; ib1 < _recv_nBlockByProc[1]; ib1++) {
    //         for (int ib0 = 0; ib0 < _recv_nBlockByProc[0]; ib0++) {
    //             // recv
    //             const int    recv_bid[3] = {ib0 , ib1 , ib2};
    //             const size_t recv_bindex = blockID(ib0, ib1, ib2, _recv_nBlockByProc);
    //             if (_recvbuf[recv_bindex] != NULL) fftw_free(_recvbuf[recv_bindex]);
    //         }
    //     }
    // }
}


/**
 * @brief execute the switch from one topo to another
 * 
 * #### Buffer writting
 * The buffer memory writting is done according to the axis of the input topologies.
 * This allows to have a continuous memory access while filling the buffer.
 * 
 * Moreover, using naxis, we know that one a data has to be send to a proc P,
 * the following (naxis-1) have the same proc P as destination.
 * So, we are able to write chunks of size naxis-1
 * 
 * Since the writting of buffers is aligned with the topo_in axis, the loops are continuous in memory and fully vectorizable.
 * 
 * #### Buffer reading
 * The buffer reading has to follow the same order as in the buffer writting, so the axis of the topo_in in the inner loop.
 * Similarly to the writting case, when one data comes from proc P, the naxis-1 following data comes from the same proc.
 * 
 * The reading of the buffer is hence continuous but the writting inside the memory has an apriori unkown stride.
 * The stride may be computed using the difference of axis between the two topologies.
 * Hence the reading will be a bit slower since the writting due to memory discontinuities
 * 
 * 
 * @param v the memory to switch from one topo to another. It has to be large enough to contain both local data's
 * @param sign if the switch is forward (UP_FORWARD) or backward (UP_BACKWARD) w.r.t. the order defined at init.
 * 
 * -----------------------------------------------
 * We do the following:
 */
void SwitchTopo::execute(opt_double_ptr v, const int sign) {
    BEGIN_FUNC

    UP_CHECK0(_topo_in->isComplex() == _topo_out->isComplex(),
              "both topologies have to be complex or real");

    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // opt_int_ptr count = _count;

    //-------------------------------------------------------------------------
    /** - set the size counters*/
    //-------------------------------------------------------------------------
    // int istart[3];
    // int ostart[3];
    // int iend[3];
    // int oend[3];
    // int ishift[3];
    // int oshift[3];
    // int send_nBlockByProc[3];
    // int recv_nBlockByProc[3];

    const Topology* topo_in;
    const Topology* topo_out;

    MPI_Request* sendRequest;
    MPI_Request* recvRequest;

    int *destRank;
    int *destTag;
    int *origRank;

    int send_nBlock[3];
    int recv_nBlock[3];

    int istart[3];
    int ostart[3];
    int iend[3];
    int oend[3];
    int inloc[3];
    int onloc[3];

    // bcoord* send_blockSize;
    // bcoord* recv_blockSize;

    

    // opt_int_ptr ib2o_shift;
    // opt_int_ptr ob2i_shift;

    // opt_int_ptr    ssend;
    // opt_int_ptr    srecv;
    // opt_int_ptr    nsend;
    // opt_int_ptr    nrecv;
    opt_double_ptr* sendBuf;
    opt_double_ptr* recvBuf;
    

    // opt_int_ptr naxis_send;
    // opt_int_ptr naxis_recv;

    if (sign == UP_FORWARD) {
        topo_in  = _topo_in;
        topo_out = _topo_out;

        sendRequest = _sendRequest;
        recvRequest = _recvRequest;

        destRank = _i2o_destRank;
        destTag = _i2o_destTag;
        origRank = _o2i_destRank;

        sendBuf = _sendBuf;
        recvBuf = _recvBuf;

        for (int id = 0; id < 3; id++) {

            send_nBlock[id] = _inBlock[id];
            recv_nBlock[id] = _onBlock[id];
            // input parameters
            istart[id] = _istart[id];
            iend[id]   = _iend[id];
            // ishift[id] = _ishift[id];
            // output params
            ostart[id] = _ostart[id];
            oend[id]   = _oend[id];
            // oshift[id] = _oshift[id];

            inloc[id] = _topo_in->nloc(id);
            onloc[id] = _topo_out->nloc(id);

            // send_nBlockByProc[id] = _send_nBlockByProc[id];
            // recv_nBlockByProc[id] = _recv_nBlockByProc[id];
        }
        // ssend      = _ssend;
        // srecv      = _srecv;
        // nsend      = _nsend;
        // nrecv      = _nrecv;
        // sendbuf    = _sendbuf;
        // recvbuf    = _recvbuf;
        // naxis_send = _naxis_i2o_send;
        // naxis_recv = _naxis_i2o_recv;

        // send_blockSize = _send_blockSize;
        // recv_blockSize = _recv_blockSize;

        // send_request = _send_request;
        // recv_request = _recv_request;

        

        // ib2o_shift = _ib2o_shift;
        // ob2i_shift = _ob2i_shift;
        

    } else if (sign == UP_BACKWARD) {
        topo_in  = _topo_out;
        topo_out = _topo_in;

        sendRequest = _recvRequest;
        recvRequest = _sendRequest;

        destRank = _o2i_destRank;
        destTag  = _o2i_destTag;
        origRank = _i2o_destRank;

        sendBuf = _recvBuf;
        recvBuf = _sendBuf;

        for (int id = 0; id < 3; id++) {

            send_nBlock[id] = _onBlock[id];
            recv_nBlock[id] = _inBlock[id];
            // // input parameters
            istart[id] = _ostart[id];
            iend[id]   = _oend[id];
            // ishift[id] = _oshift[id];
            // // output params
            ostart[id] = _istart[id];
            oend[id]   = _iend[id];

            inloc[id] = _topo_out->nloc(id);
            onloc[id] = _topo_in->nloc(id);
            // oshift[id] = _ishift[id];

            // send_nBlockByProc[id] = _recv_nBlockByProc[id];
            // recv_nBlockByProc[id] = _send_nBlockByProc[id];
        }
        // ssend      = _srecv;
        // srecv      = _ssend;
        // nsend      = _nrecv;
        // nrecv      = _nsend;
        // sendbuf    = _recvBuf;
        // recvbuf    = _sendBuf;
        // naxis_send = _naxis_o2i_send;
        // naxis_recv = _naxis_o2i_recv;

        // send_blockSize = _recv_blockSize;
        // recv_blockSize = _send_blockSize;

        // send_request = _recv_request;
        // recv_request = _send_request;

        // ib2o_shift = _ob2i_shift;
        // ob2i_shift = _ib2o_shift;
        
    } else {
        UP_CHECK0(false, "the sign is not UP_FORWARD nor UP_BACKWARD");
    }

    INFOLOG5("previous topo: %d,%d,%d axis=%d\n", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    INFOLOG5("new topo: %d,%d,%d  axis=%d\n", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());

    const int ax0 = topo_in->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nf  = topo_in->nf();

    // printf("-------- BEGIN ----------\n");

    // printf("sending %d blocks\n",send_nBlock[0]*send_nBlock[1]*send_nBlock[2]);
    // printf("recving %d blocks\n",recv_nBlock[0]*recv_nBlock[1]*recv_nBlock[2]);
    

    //-------------------------------------------------------------------------
    /** - generate the reception requests so we are ready to receive */
    //-------------------------------------------------------------------------
    for (int ib2 = 0; ib2 < recv_nBlock[ax2]; ib2++) {
        for (int ib1 = 0; ib1 < recv_nBlock[ax1]; ib1++) {
            for (int ib0 = 0; ib0 < recv_nBlock[ax0]; ib0++) {
                // get the block ID
                const int      bid      = localIndex(ax0,ib0,ib1,ib2,0,recv_nBlock,1);
                opt_double_ptr data     = recvBuf[bid];
                const int      datasize = _nByBlock[0] * _nByBlock[1] * _nByBlock[2] * nf;
                // generate the request
                MPI_Irecv(data, datasize, MPI_DOUBLE, origRank[bid], MPI_ANY_TAG, MPI_COMM_WORLD, &(recvRequest[bid]));
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    // std::memset(count, 0, sizeof(int) * comm_size);

    for (int ib2 = 0; ib2 < send_nBlock[ax2]; ib2++) {
        for (int ib1 = 0; ib1 < send_nBlock[ax1]; ib1++) {
            for (int ib0 = 0; ib0 < send_nBlock[ax0]; ib0++) {
                const int      bid  = localIndex(ax0, ib0, ib1, ib2, 0, send_nBlock, 1);
                opt_double_ptr data = sendBuf[bid];

                // go inside the block
                for (int i2 = 0; i2 < _nByBlock[ax2]; i2++) {
                    for (int i1 = 0; i1 < _nByBlock[ax1]; i1++) {
                        // get the starting id for the buffer
                        const size_t buf_idx = localIndex(ax0, 0, i1, i2, ax0, _nByBlock, nf);
                        // get the starting id for the domain
                        const int    loci0  = istart[ax0] + ib0 * _nByBlock[ax0] + 0;
                        const int    loci1  = istart[ax1] + ib1 * _nByBlock[ax1] + i1;
                        const int    loci2  = istart[ax2] + ib2 * _nByBlock[ax2] + i2;
                        const size_t my_idx = localIndex(ax0, loci0, loci1, loci2, ax0, inloc, nf);
                        // get the max counter
                        const size_t nmax = _nByBlock[ax0] * nf;
                        // do the copy
                        for (size_t i0 = 0; i0 < nmax; i0++) {
                            data[buf_idx + i0] = v[my_idx + i0];
                        }
                    }
                }
                // printf("sending block %d to proc %d with tag %d\n",bid,destRank[bid],destTag[bid]);
                // send the block and continue
                const int datasize = _nByBlock[0] * _nByBlock[1] * _nByBlock[2] * topo_in->nf();
                MPI_Isend(data, datasize, MPI_DOUBLE, destRank[bid], destTag[bid], MPI_COMM_WORLD, &(sendRequest[bid]));
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - wait for a block and copy when it arrives */
    //-------------------------------------------------------------------------
    // reset the memory to 0
    std::memset(v, 0, sizeof(double) * topo_out->locmemsize());

    const int nblocks = recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2];
    int request_index;
    MPI_Status status;
    
    for(int count=0; count<nblocks; count++){
        // wait for a block
        MPI_Waitany(nblocks, recvRequest, &request_index,&status);
        
        // get the block id = the tag
        int bid = status.MPI_TAG;
        // get the indexing of the block in 012-indexing
        int ibv[3];
        localSplit(bid,recv_nBlock,0,ibv);
        // printf("doing split with size = %d %d %d\n",recv_nBlock[0],recv_nBlock[1],recv_nBlock[2]);
        // printf("block ID %d become %d %d %d\n",bid,ibv[0],ibv[1],ibv[2]);
        // get the associated data
        opt_double_ptr data = recvBuf[bid];
        // // get the block size
        // const bcoord blockSize = {recv_blockSize[recv_bid][0],recv_blockSize[recv_bid][1],recv_blockSize[recv_bid][2]};
        // get the starting id in memory
        // const size_t start_idx  = localIndex_withBlock_ao(0, 0, 0, topo_out, recv_ibv, _nByBlock);

        // printf("doing block %d with nf=%d\n",bid,nf);

        // go inside the block using the axis of the topo_in
        for (int i2 = 0; i2 < _nByBlock[ax2]; i2++) {
            for (int i1 = 0; i1 < _nByBlock[ax1]; i1++) {
                // get the starting id for the buffer
                const size_t buf_idx = localIndex(ax0, 0, i1, i2, ax0, _nByBlock, nf);
                // get the starting id for the domain
                const int    loci0  = ostart[ax0] + ibv[ax0] * _nByBlock[ax0] + 0;
                const int    loci1  = ostart[ax1] + ibv[ax1] * _nByBlock[ax1] + i1;
                const int    loci2  = ostart[ax2] + ibv[ax2] * _nByBlock[ax2] + i2;
                const size_t my_idx = localIndex(ax0, loci0, loci1, loci2, topo_out->axis(), onloc, nf);
                const size_t stride = localIndex(ax0, 1, 0, 0, topo_out->axis(), onloc, nf);


                // UP_CHECK2((my_idx < onloc[0]*onloc[1]*onloc[2]),"myidx out of scope (%d/%d)",my_idx, _nByBlock[0]*_nByBlock[1]*_nByBlock[2]);

                // printf("doing the copy with my_idx = %d and stride = %d\n",my_idx,stride);
                // do the copy
                if (nf == 1) {
                    for (int i0 = 0; i0 < _nByBlock[ax0]; i0++) {
                        // UP_CHECK4((my_idx + i0 * stride >= 0),"index out of scope (%d/%d): i0=%d, stride=%d",my_idx + i0 * stride , 0,i0,stride);
                        // UP_CHECK4((my_idx + i0 * stride < onloc[0]*onloc[1]*onloc[2]),"index out of scope (%d/%d): i0=%d, stride=%d",my_idx + i0 * stride , onloc[0]*onloc[1]*onloc[2],i0,stride);
                        // UP_CHECK3((buf_idx + i0 < _nByBlock[0]*_nByBlock[1]*_nByBlock[2]),"index out of scope (%d/%d): i0=%d",buf_idx + i0 , _nByBlock[0]*_nByBlock[1]*_nByBlock[2],i0);
                        v[my_idx + i0 * stride] = data[buf_idx + i0];
                    }
                } else if (nf == 2) {
                    for (int i0 = 0; i0 < _nByBlock[ax0]; i0++) {
                        v[my_idx + i0 * stride + 0] = data[buf_idx + i0 * 2 + 0];
                        v[my_idx + i0 * stride + 1] = data[buf_idx + i0 * 2 + 1];
                    }
                }
                else{
                    UP_CHECK0(false,"the value of nf is not supported")
                }
            }
        }
    }


    // MPI_Barrier(MPI_COMM_WORLD);

    // int dest_rankd[3];
    // for (int i2 = istart[ax2]; i2 < iend[ax2]; ++i2) {
    //     dest_rankd[ax2] = topo_in->cmpt_matchrank(ax2, topo_out, i2 + ishift[ax2]);
    //     for (int i1 = istart[ax1]; i1 < iend[ax1]; ++i1) {
    //         dest_rankd[ax1] = topo_in->cmpt_matchrank(ax1, topo_out, i1 + ishift[ax1]);
    //         for (int i0 = istart[ax0]; i0 < iend[ax0]; ++i0) {
    //             dest_rankd[ax0] = topo_in->cmpt_matchrank(ax0, topo_out, i0 + ishift[ax0]);
    //             // add one info to the destination index
    //             const int    dest_rank = rankindex(dest_rankd, topo_out);
    //             const int    buf_idx   = ssend[dest_rank] + count[dest_rank];
    //             const size_t my_idx    = localindex_ao(i0, i1, i2, topo_in);

    //             for (int ib = 0; ib < naxis_send[dest_rank]*topo_in->nf(); ib++) {
    //                 sendbuf[buf_idx + ib] = v[my_idx + ib];
    //             }
    //             count[dest_rank] += naxis_send[dest_rank]*topo_in->nf();
    //             i0 = i0 + naxis_send[dest_rank] - 1;
    //         }
    //     }
    // }

    // //-------------------------------------------------------------------------
    // /** - Send it */
    // //-------------------------------------------------------------------------
    // MPI_Alltoallv(sendbuf, nsend, ssend, MPI_DOUBLE, recvbuf, nrecv, srecv, MPI_DOUBLE, MPI_COMM_WORLD);
    // // for(int ip=0; ip<comm_size; ip++){
    // //     MPI_Isend(sendbuf[ssend[ip]],nsend[ip],MPI_DOUBLE,ip)
    // // }

    //-------------------------------------------------------------------------
    /** - Fill the memory */
    //-------------------------------------------------------------------------
    

    

    // int orig_rankd[3];
    // if (topo_out->nf() == 1) {
    //     for (int i2 = ostart[ax2]; i2 < oend[ax2]; ++i2) {
    //         orig_rankd[ax2] = topo_out->cmpt_matchrank(ax2, topo_in, i2 + oshift[ax2]);
    //         for (int i1 = ostart[ax1]; i1 < oend[ax1]; ++i1) {
    //             orig_rankd[ax1] = topo_out->cmpt_matchrank(ax1, topo_in, i1 + oshift[ax1]);
    //             for (int i0 = ostart[ax0]; i0 < oend[ax0]; ++i0) {
    //                 orig_rankd[ax0] = topo_out->cmpt_matchrank(ax0, topo_in, i0 + oshift[ax0]);

    //                 const int    orig_rank = rankindex(orig_rankd, topo_in);
    //                 const size_t buf_idx   = srecv[orig_rank] + count[orig_rank];
    //                 // get the starting index and the stride obtained when adding one point in i0
    //                 const size_t my_idx = localindex(topo_in->axis(), i0, i1, i2, topo_out);
    //                 const size_t stride = localindex(topo_in->axis(), 1, 0, 0, topo_out);

    //                 for (int ib = 0; ib < naxis_recv[orig_rank]; ib++) {
    //                     v[my_idx + ib * stride] = recvbuf[buf_idx + ib];
    //                 }
    //                 count[orig_rank] += naxis_recv[orig_rank];

    //                 i0 = i0 + naxis_recv[orig_rank] - 1;
    //             }
    //         }
    //     }
    // } else if (topo_out->nf() == 2) {
    //     for (int i2 = ostart[ax2]; i2 < oend[ax2]; ++i2) {
    //         orig_rankd[ax2] = topo_out->cmpt_matchrank(ax2, topo_in, i2 + oshift[ax2]);
    //         for (int i1 = ostart[ax1]; i1 < oend[ax1]; ++i1) {
    //             orig_rankd[ax1] = topo_out->cmpt_matchrank(ax1, topo_in, i1 + oshift[ax1]);
    //             for (int i0 = ostart[ax0]; i0 < oend[ax0]; ++i0) {
    //                 orig_rankd[ax0] = topo_out->cmpt_matchrank(ax0, topo_in, i0 + oshift[ax0]);

    //                 // get the rank of origin
    //                 const int orig_rank = rankindex(orig_rankd, topo_in);
    //                 // get the starting buffer index
    //                 const size_t buf_idx = srecv[orig_rank] + count[orig_rank];
    //                 // get the starting index and the stride obtained when adding one point in i0
    //                 const size_t my_idx = localindex(topo_in->axis(), i0, i1, i2, topo_out);
    //                 const size_t stride = localindex(topo_in->axis(), 1, 0, 0, topo_out);

    //                 // for every element
    //                 for (int ib = 0; ib < naxis_recv[orig_rank]; ib++) {
    //                     v[my_idx + ib * stride + 0] = recvbuf[buf_idx + ib * 2 + 0];
    //                     v[my_idx + ib * stride + 1] = recvbuf[buf_idx + ib * 2 + 1];
    //                 }
    //                 count[orig_rank] += naxis_recv[orig_rank] * 2;
    //                 i0 = i0 + naxis_recv[orig_rank] - 1;
    //             }
    //         }
    //     }
    // } else {
    //     UP_CHECK0(false, "size of Topological nf not supported");
    // }
}

void SwitchTopo::disp() {
    BEGIN_FUNC
    INFO("------------------------------------------\n");
    INFO("## Topo Swticher MPI\n");
    INFO("--- INPUT\n");
    INFO2("  - input axis = %d\n", _topo_in->axis());
    INFO4("  - input local = %d %d %d\n", _topo_in->nloc(0), _topo_in->nloc(1), _topo_in->nloc(2));
    INFO4("  - input global = %d %d %d\n", _topo_in->nglob(0), _topo_in->nglob(1), _topo_in->nglob(2));
    INFO4("  - ishift = %d %d %d\n", _ishift[0], _ishift[1], _ishift[2]);
    INFO4("  - istart = %d %d %d\n", _istart[0], _istart[1], _istart[2]);
    INFO4("  - iend = %d %d %d\n", _iend[0], _iend[1], _iend[2]);
    INFO("--- OUTPUT\n");
    INFO2("  - output axis = %d\n", _topo_out->axis());
    INFO4("  - output local = %d %d %d\n", _topo_out->nloc(0), _topo_out->nloc(1), _topo_out->nloc(2));
    INFO4("  - output global = %d %d %d\n", _topo_out->nglob(0), _topo_out->nglob(1), _topo_out->nglob(2));
    INFO4("  - oshift = %d %d %d\n", _oshift[0], _oshift[1], _oshift[2]);
    INFO4("  - ostart = %d %d %d\n", _ostart[0], _ostart[1], _ostart[2]);
    INFO4("  - oend = %d %d %d\n", _oend[0], _oend[1], _oend[2]);
    INFO("------------------------------------------\n");
}

void SwitchTopo_test() {
    BEGIN_FUNC

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {8, 8, 8};
    // const int nproc[3] = {comm_size, 1, 1};
    const int nproc[3] = {2, 2, 1};

    const int nglob_big[3] = {17, 8, 8};
    const int nproc_big[3] = {2, 2, 1};

    //===========================================================================
    // real numbers
    Topology* topo    = new Topology(0, nglob, nproc, false);
    Topology* topobig = new Topology(0, nglob_big, nproc_big, false);

    // printf("nloc topo = %d %d %d\n",topo->nloc(0),topo->nloc(1),topo->nloc(2));
    // printf("nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));

    double* data = (double*)fftw_malloc(sizeof(double*) * std::max(topo->locmemsize(),topobig->locmemsize()));

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

    // topobig->switch2complex();
    // printf("as complex: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));
    // topobig->switch2real();
    // printf("as real: nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));

    const int   fieldstart[3] = {0, 0, 0};
    printf("\n\n=============================\n");
    SwitchTopo* switchtopo    = new SwitchTopo(topo, topobig, fieldstart);

    printf("\n\n============ FORWARD =================\n");
    switchtopo->execute(data, UP_FORWARD);

    hdf5_dump(topobig, "test_real_padd", data);

    printf("\n\n============ BACKWARD =================\n");
    switchtopo->execute(data, UP_BACKWARD);

    hdf5_dump(topo, "test_real_returned", data);

    fftw_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);
    

    //===========================================================================
    // complex numbers
    topo    = new Topology(0, nglob, nproc, true);
    topobig = new Topology(2, nglob_big, nproc_big, true);

    // printf("nloc topo = %d %d %d\n",topo->nloc(0),topo->nloc(1),topo->nloc(2));
    // printf("nloc topobig = %d %d %d\n",topobig->nloc(0),topobig->nloc(1),topobig->nloc(2));

    data = (double*)fftw_malloc(sizeof(double*) * topobig->locmemsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localindex_xyz(i0, i1, i2, topo);
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

    const int   fieldstart2[3] = {4, 0, 0};
    printf("\n=============================\n");
    switchtopo    = new SwitchTopo(topo, topobig, fieldstart2);

    switchtopo->execute(data, UP_FORWARD);

    hdf5_dump(topobig, "test_complex_padd", data);

    switchtopo->execute(data, UP_BACKWARD);

    hdf5_dump(topo, "test_complex_returned", data);

    fftw_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);
    
}