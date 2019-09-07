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
SwitchTopo::SwitchTopo(const Topology* topo_input, const Topology* topo_output, const int shift[3], Profiler* prof) {
    BEGIN_FUNC

    FLUPS_CHECK0(topo_input->isComplex() == topo_output->isComplex(), "both topologies have to be the same kind");

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
    for (int id = 0; id < 3; id++) {
        _ib2o_shift[id] = shift[id];
        _ob2i_shift[id] = -shift[id];
    }
    // get how much pints we send/recv in each direction
    _topo_in->cmpt_intersect_id(_ib2o_shift, _topo_out, _istart, _iend);
    _topo_out->cmpt_intersect_id(_ob2i_shift, _topo_in, _ostart, _oend);

    //-------------------------------------------------------------------------
    /** - get the block size as the GCD of the memory among every process between send and receive */
    //-------------------------------------------------------------------------
    // We use the greatest common divisor, because it is possible that the last proc
    // in a given direction has a bit more data than the others (which all have
    // _nbyproc points).
    int* onProc = (int*)fftw_malloc(comm_size * sizeof(int));
    for (int id = 0; id < 3; id++) {
        // get the gcd between send and receive
        int npoints = gcd(_iend[id] - _istart[id], _oend[id] - _ostart[id]);
        // gather on each proc the gcd
        MPI_Allgather(&npoints, 1, MPI_INT, onProc, 1, MPI_INT, MPI_COMM_WORLD);
        // get the Greatest Common Divider among every process
        int my_gcd = onProc[0];
        for (int ip = 1; ip < comm_size; ip++) {
            my_gcd = gcd(my_gcd, onProc[ip]);
        }
        // store it as the block dimension
        _nByBlock[id] = my_gcd;
    }
    fftw_free(onProc);

    //-------------------------------------------------------------------------
    /** - get the starting index of the block 0,0,0 for input and output */
    //-------------------------------------------------------------------------
    int* inBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));
    int* onBlockEachProc = (int*)fftw_malloc(comm_size * 3 * sizeof(int));

    cmpt_blockIndexes(_istart, _iend, _nByBlock, _topo_in, _inBlock, _iblockIDStart, inBlockEachProc);
    cmpt_blockIndexes(_ostart, _oend, _nByBlock, _topo_out, _onBlock, _oblockIDStart, onBlockEachProc);

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
    _i2o_destRank = (opt_int_ptr)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
    _o2i_destRank = (opt_int_ptr)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));
    // allocate the destination tags
    _i2o_destTag = (opt_int_ptr)fftw_malloc(_inBlock[0] * _inBlock[1] * _inBlock[2] * sizeof(int));
    _o2i_destTag = (opt_int_ptr)fftw_malloc(_onBlock[0] * _onBlock[1] * _onBlock[2] * sizeof(int));

    //-------------------------------------------------------------------------
    /** - for each block, get the destination rank */
    //-------------------------------------------------------------------------
    // send destination ranks in the ouput topo
    cmpt_blockDestRankAndTag(_inBlock, _iblockIDStart, _topo_out, onBlockEachProc, _i2o_destRank, _i2o_destTag);
    cmpt_blockDestRankAndTag(_onBlock, _oblockIDStart, _topo_in, inBlockEachProc, _o2i_destRank, _o2i_destTag);

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
                const size_t send_bid = localIndex(0, ib0, ib1, ib2, 0, _inBlock, 1);
                // allocate at the correct size
                _sendBuf[send_bid] = (double*)fftw_malloc(_nByBlock[0] * _nByBlock[1] * _nByBlock[2] * sizeof(double) * _topo_in->nf());
                FLUPS_CHECK1(FLUPS_ISALIGNED(_sendBuf[send_bid]), "FFTW alignement not compatible with FLUPS_ALIGNMENT (=%d)", FLUPS_ALIGNMENT);
            }
        }
    }
    for (int ib2 = 0; ib2 < _onBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < _onBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < _onBlock[0]; ib0++) {
                // get the local block index
                const size_t recv_bid = localIndex(0, ib0, ib1, ib2, 0, _onBlock, 1);
                // allocate at the correct size
                _recvBuf[recv_bid] = (double*)fftw_malloc(_nByBlock[0] * _nByBlock[1] * _nByBlock[2] * sizeof(double) * _topo_in->nf());
                FLUPS_CHECK1(FLUPS_ISALIGNED(_recvBuf[recv_bid]), "FFTW alignement not compatible with FLUPS_ALIGNMENT (=%d)", FLUPS_ALIGNMENT);
            }
        }
    }
    //-------------------------------------------------------------------------
    /** - setup the profiler    */
    //-------------------------------------------------------------------------
    if (_prof != NULL) {
        _prof->create("reorder_mem2buf");
        _prof->create("reorder_buf2mem");
        _prof->create("reorder_waiting");
    }
}

/**
 * @brief Destroy the Switch Topo
 * 
 */
SwitchTopo::~SwitchTopo() {
    if (_i2o_destRank != NULL) fftw_free(_i2o_destRank);
    if (_o2i_destRank != NULL) fftw_free(_o2i_destRank);
    if (_i2o_destTag != NULL) fftw_free(_i2o_destTag);
    if (_o2i_destTag != NULL) fftw_free(_o2i_destTag);

    if (_sendRequest != NULL) fftw_free(_sendRequest);
    if (_recvRequest != NULL) fftw_free(_recvRequest);

    for (int ib2 = 0; ib2 < _inBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < _inBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < _inBlock[0]; ib0++) {
                // get the local block index
                const size_t send_bid = localIndex(0, ib0, ib1, ib2, 0, _inBlock, 1);
                if (_sendBuf[send_bid] != NULL) fftw_free(_sendBuf[send_bid]);
            }
        }
    }
    for (int ib2 = 0; ib2 < _onBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < _onBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < _onBlock[0]; ib0++) {
                // get the local block index
                const size_t recv_bid = localIndex(0, ib0, ib1, ib2, 0, _onBlock, 1);
                if (_recvBuf[recv_bid] != NULL) fftw_free(_recvBuf[recv_bid]);
            }
        }
    }
    fftw_free((double**)_sendBuf);
    fftw_free((double**)_recvBuf);
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
void SwitchTopo::execute(opt_double_ptr v, const int sign) {
    BEGIN_FUNC

    FLUPS_CHECK0(_topo_in->isComplex() == _topo_out->isComplex(),
              "both topologies have to be complex or real");

    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //-------------------------------------------------------------------------
    /** - setup required memory arrays */
    //-------------------------------------------------------------------------

    const Topology* topo_in;
    const Topology* topo_out;

    MPI_Request* sendRequest;
    MPI_Request* recvRequest;

    opt_int_ptr destRank;
    opt_int_ptr destTag;
    opt_int_ptr origRank;

    int send_nBlock[3];
    int recv_nBlock[3];

    int istart[3];
    int ostart[3];
    int iend[3];
    int oend[3];
    int inloc[3];
    int onloc[3];

    opt_double_ptr* sendBuf;
    opt_double_ptr* recvBuf;

    if (sign == FLUPS_FORWARD) {
        topo_in     = _topo_in;
        topo_out    = _topo_out;
        sendRequest = _sendRequest;
        recvRequest = _recvRequest;
        destRank    = _i2o_destRank;
        destTag     = _i2o_destTag;
        origRank    = _o2i_destRank;
        sendBuf     = _sendBuf;
        recvBuf     = _recvBuf;

        for (int id = 0; id < 3; id++) {
            send_nBlock[id] = _inBlock[id];
            recv_nBlock[id] = _onBlock[id];
            istart[id]      = _istart[id];
            iend[id]        = _iend[id];
            ostart[id]      = _ostart[id];
            oend[id]        = _oend[id];
            inloc[id]       = _topo_in->nloc(id);
            onloc[id]       = _topo_out->nloc(id);
        }
    } else if (sign == FLUPS_BACKWARD) {
        topo_in     = _topo_out;
        topo_out    = _topo_in;
        sendRequest = _recvRequest;
        recvRequest = _sendRequest;
        destRank    = _o2i_destRank;
        destTag     = _o2i_destTag;
        origRank    = _i2o_destRank;
        sendBuf     = _recvBuf;
        recvBuf     = _sendBuf;

        for (int id = 0; id < 3; id++) {
            send_nBlock[id] = _onBlock[id];
            recv_nBlock[id] = _inBlock[id];
            istart[id]      = _ostart[id];
            iend[id]        = _oend[id];
            ostart[id]      = _istart[id];
            oend[id]        = _iend[id];
            inloc[id]       = _topo_out->nloc(id);
            onloc[id]       = _topo_in->nloc(id);
        }
    } else {
        FLUPS_CHECK0(false, "the sign is not FLUPS_FORWARD nor FLUPS_BACKWARD");
    }

    INFOLOG5("previous topo: %d,%d,%d axis=%d\n", topo_in->nglob(0), topo_in->nglob(1), topo_in->nglob(2), topo_in->axis());
    INFOLOG5("new topo: %d,%d,%d  axis=%d\n", topo_out->nglob(0), topo_out->nglob(1), topo_out->nglob(2), topo_out->axis());

    const int ax0 = topo_in->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nf  = topo_in->nf();

    //-------------------------------------------------------------------------
    /** - generate the reception requests so we are ready to receive */
    //-------------------------------------------------------------------------
    for (int ib2 = 0; ib2 < recv_nBlock[ax2]; ib2++) {
        for (int ib1 = 0; ib1 < recv_nBlock[ax1]; ib1++) {
            for (int ib0 = 0; ib0 < recv_nBlock[ax0]; ib0++) {
                // get the block ID
                const int      bid      = localIndex(ax0, ib0, ib1, ib2, 0, recv_nBlock, 1);
                opt_double_ptr data     = recvBuf[bid];
                const int      datasize = _nByBlock[0] * _nByBlock[1] * _nByBlock[2] * nf;
                // generate the request
                MPI_Irecv(data, datasize, MPI_DOUBLE, origRank[bid], MPI_ANY_TAG, MPI_COMM_WORLD, &(recvRequest[bid]));
            }
        }
    }


    if (_prof != NULL) {
        _prof->start("reorder_mem2buf");
    }
    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
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
                // send the block and continue
                const int datasize = _nByBlock[0] * _nByBlock[1] * _nByBlock[2] * topo_in->nf();
                MPI_Isend(data, datasize, MPI_DOUBLE, destRank[bid], destTag[bid], MPI_COMM_WORLD, &(sendRequest[bid]));
            }
        }
    }
    if (_prof != NULL) {
        _prof->stop("reorder_mem2buf");
    }

    //-------------------------------------------------------------------------
    /** - wait for a block and copy when it arrives */
    //-------------------------------------------------------------------------
    // reset the memory to 0
    std::memset(v, 0, sizeof(double) * topo_out->locmemsize());
    // get some counters
    const int  nblocks = recv_nBlock[0] * recv_nBlock[1] * recv_nBlock[2];
    int        request_index;
    MPI_Status status;
    // for each block
    if (_prof != NULL) {
        _prof->start("reorder_buf2mem");
    }
    for (int count = 0; count < nblocks; count++) {
        // wait for a block
        if (_prof != NULL) {
            _prof->start("reorder_waiting");
        }
        MPI_Waitany(nblocks, recvRequest, &request_index, &status);
        if (_prof != NULL) {
            _prof->stop("reorder_waiting");
        }
        // get the block id = the tag
        int bid = status.MPI_TAG;
        // get the indexing of the block in 012-indexing
        int ibv[3];
        localSplit(bid, recv_nBlock, 0, ibv);
        // get the associated data
        opt_double_ptr data = recvBuf[bid];
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
                // do the copy
                if (nf == 1) {
                    for (int i0 = 0; i0 < _nByBlock[ax0]; i0++) {
                        v[my_idx + i0 * stride] = data[buf_idx + i0];
                    }
                } else if (nf == 2) {
                    for (int i0 = 0; i0 < _nByBlock[ax0]; i0++) {
                        v[my_idx + i0 * stride + 0] = data[buf_idx + i0 * 2 + 0];
                        v[my_idx + i0 * stride + 1] = data[buf_idx + i0 * 2 + 1];
                    }
                } else {
                    FLUPS_CHECK0(false, "the value of nf is not supported")
                }
            }
        }
    }
    if (_prof != NULL) {
        _prof->stop("reorder_buf2mem");
    }
}

void SwitchTopo::disp() {
    BEGIN_FUNC
    INFO("------------------------------------------\n");
    INFO("## Topo Swticher MPI\n");
    INFO("--- INPUT\n");
    INFO2("  - input axis = %d\n", _topo_in->axis());
    INFO4("  - input local = %d %d %d\n", _topo_in->nloc(0), _topo_in->nloc(1), _topo_in->nloc(2));
    INFO4("  - input global = %d %d %d\n", _topo_in->nglob(0), _topo_in->nglob(1), _topo_in->nglob(2));
    INFO4("  - istart = %d %d %d\n", _istart[0], _istart[1], _istart[2]);
    INFO4("  - iend = %d %d %d\n", _iend[0], _iend[1], _iend[2]);
    INFO("--- OUTPUT\n");
    INFO2("  - output axis = %d\n", _topo_out->axis());
    INFO4("  - output local = %d %d %d\n", _topo_out->nloc(0), _topo_out->nloc(1), _topo_out->nloc(2));
    INFO4("  - output global = %d %d %d\n", _topo_out->nglob(0), _topo_out->nglob(1), _topo_out->nglob(2));
    INFO4("  - ostart = %d %d %d\n", _ostart[0], _ostart[1], _ostart[2]);
    INFO4("  - oend = %d %d %d\n", _oend[0], _oend[1], _oend[2]);
    INFO("------------------------------------------\n");
}

void SwitchTopo_test() {
    BEGIN_FUNC

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {8, 8, 8};
    const int nproc[3] = {2, 2, 1};

    const int nglob_big[3] = {17, 8, 8};
    const int nproc_big[3] = {2, 2, 1};

    //===========================================================================
    // real numbers
    Topology* topo    = new Topology(0, nglob, nproc, false);
    Topology* topobig = new Topology(0, nglob_big, nproc_big, false);

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
    // printf("\n\n=============================\n");
    SwitchTopo* switchtopo = new SwitchTopo(topo, topobig, fieldstart, NULL);

    // printf("\n\n============ FORWARD =================\n");
    switchtopo->execute(data, FLUPS_FORWARD);

    hdf5_dump(topobig, "test_real_padd", data);

    // printf("\n\n============ BACKWARD =================\n");
    switchtopo->execute(data, FLUPS_BACKWARD);

    hdf5_dump(topo, "test_real_returned", data);

    fftw_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);

    //===========================================================================
    // complex numbers
    topo    = new Topology(0, nglob, nproc, true);
    topobig = new Topology(2, nglob_big, nproc_big, true);

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

    const int fieldstart2[3] = {4, 0, 0};
    // printf("\n=============================\n");
    switchtopo = new SwitchTopo(topo, topobig, fieldstart2, NULL);

    switchtopo->execute(data, FLUPS_FORWARD);

    hdf5_dump(topobig, "test_complex_padd", data);

    switchtopo->execute(data, FLUPS_BACKWARD);

    hdf5_dump(topo, "test_complex_returned", data);

    fftw_free(data);
    delete (switchtopo);
    delete (topo);
    delete (topobig);
}