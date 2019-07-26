/**
 * @file Reorder_MPI.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-25
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "Reorder_MPI.hpp"

inline static size_t localindex(const int i[3], const int n[3], const int axis)
{
    int axis1 = (axis + 1) % 3;
    int axis2 = (axis + 2) % 3;
    return i[axis] + n[axis] * (i[axis1] + n[axis1] * i[axis2]);
}
inline static int rankindex(const int rankd[3], const int nproc[3])
{
    return rankd[0] + nproc[0] * (rankd[1] + nproc[1] * rankd[2]);
}
inline static int rankindex(const int rankd_0, const int rankd_1, const int rankd_2, const int nproc_0, const int nproc_1, const int nproc_2)
{
    return rankd_0 + nproc_0 * (rankd_1 + nproc_1 * rankd_2);
}
inline static void ranksplit(const int rank, const int nproc[3], int rankd[3])
{
    rankd[0] = rank % nproc[0];
    rankd[1] = (rank % (nproc[0] * nproc[1])) / nproc[0];
    rankd[2] = rank / (nproc[0] * nproc[1]);
}

Reorder_MPI::Reorder_MPI(const int nglob[3], const int mynf, const int inproc[3], int myaxis0, int onproc[3], int myaxis1, double* mybufs[2]):
 _iaxis(myaxis0), _oaxis(myaxis1), _nf(mynf)
{
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    UP_CHECK2(nglob[0] % inproc[0] == 0, "unable to divide %d unknows into %d procs", nglob[0], inproc[0]);
    UP_CHECK2(nglob[1] % inproc[1] == 0, "unable to divide %d unknows into %d procs", nglob[1], inproc[1]);
    UP_CHECK2(nglob[2] % inproc[2] == 0, "unable to divide %d unknows into %d procs", nglob[2], inproc[2]);

    //-------------------------------------------------------------------------
    /** - get the rankd for input and output  */
    //-------------------------------------------------------------------------
    ranksplit(rank, inproc, _irankd);
    ranksplit(rank, onproc, _orankd);

    for (int id = 0; id < 3; id++)
    {
        _inproc[id] = inproc[id];
        _onproc[id] = onproc[id];
        _inloc[id] = nglob[id] / inproc[id];
        _onloc[id] = nglob[id] / onproc[id];
    }

    // sanity checks
    UP_CHECK2(_inloc[0] * _inloc[1] * _inloc[2] == _onloc[0] * _onloc[1] * _onloc[2], "Memory has changed from %d to %d", _inloc[0] * _inloc[1] * _inloc[2], _onloc[0] * _onloc[1] * _onloc[2]);

    //-------------------------------------------------------------------------
    /** - allocate the size arrays */
    //-------------------------------------------------------------------------
    _nsend = (int *)fftw_malloc(comm_size * sizeof(int));
    _nrecv = (int *)fftw_malloc(comm_size * sizeof(int));
    _ssend = (int *)fftw_malloc(comm_size * sizeof(int));
    _srecv = (int *)fftw_malloc(comm_size * sizeof(int));
    _count = (int *)fftw_malloc(comm_size * sizeof(int));

    std::memset(_nsend, 0, sizeof(int) * comm_size);
    std::memset(_nrecv, 0, sizeof(int) * comm_size);
    std::memset(_ssend, 0, sizeof(int) * comm_size);
    std::memset(_srecv, 0, sizeof(int) * comm_size);
    std::memset(_count, 0, sizeof(int) * comm_size);

    //-------------------------------------------------------------------------
    /** - allocate the buffer if needed */
    //-------------------------------------------------------------------------
    if (_bufsend == NULL)
    {
        INFOLOG("allocating memory for the buffers");
        _bufsend = (double *)fftw_malloc(sizeof(double) * _nf * _inloc[0] * _inloc[1] * _inloc[2]);
        _bufrecv = (double *)fftw_malloc(sizeof(double) * _nf * _onloc[0] * _onloc[1] * _onloc[2]);

        mybufs[0] = _bufsend;
        mybufs[1] = _bufrecv;

        do_deallocate = true;
    }
    else
    {
        UP_CHECK0(UP_ISALIGNED(mybufs[0]),"myfufs[0] not aligned");
        UP_CHECK0(UP_ISALIGNED(mybufs[1]),"myfufs[0] not aligned");

        _bufsend = mybufs[0];
        _bufrecv = mybufs[1];
    }

    //-------------------------------------------------------------------------
    /** - for each data get its destination rank */
    //-------------------------------------------------------------------------
    int dest_rankd[3];
    for (int i2 = 0; i2 < _inloc[2]; ++i2)
    {
        dest_rankd[2] = (_irankd[2] * _inloc[2] + i2) / _onloc[2];
        for (int i1 = 0; i1 < _inloc[1]; ++i1)
        {
            dest_rankd[1] = (_irankd[1] * _inloc[1] + i1) / _onloc[1];
            for (int i0 = 0; i0 < _inloc[0]; ++i0)
            {
                dest_rankd[0] = (_irankd[0] * _inloc[0] + i0) / _onloc[0];
                _nsend[rankindex(dest_rankd, onproc)] += _nf;
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - send the size to each proc and store it */
    //-------------------------------------------------------------------------
    MPI_Alltoall(_nsend, 1, MPI_INT, _nrecv, 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 1; i < comm_size; ++i)
    {
        _ssend[i] = _ssend[i - 1] + _nsend[i - 1];
        _srecv[i] = _srecv[i - 1] + _nrecv[i - 1];
    }
}

Reorder_MPI::~Reorder_MPI()
{
    if (_nsend != NULL)
        fftw_free(_nsend);
    if (_nrecv != NULL)
        fftw_free(_nrecv);
    if (_ssend != NULL)
        fftw_free(_ssend);
    if (_srecv != NULL)
        fftw_free(_srecv);
    if (_srecv != NULL)
        fftw_free(_count);
    if (_bufsend != NULL && do_deallocate)
        fftw_free(_bufsend);
    if (_bufrecv != NULL && do_deallocate)
        fftw_free(_bufrecv);
}

void Reorder_MPI::execute(opt_double_ptr v, const int sign)
{
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    opt_int_ptr ssend = _ssend;
    opt_int_ptr srecv = _srecv;
    opt_int_ptr nsend = _nsend;
    opt_int_ptr nrecv = _nrecv;
    opt_int_ptr count = _count;
    opt_double_ptr recvbuf = _bufrecv;
    opt_double_ptr sendbuf = _bufsend;

    //-------------------------------------------------------------------------
    /** - set the size counters*/
    //-------------------------------------------------------------------------
    int iaxis, oaxis;
    int inloc[3];
    int onloc[3];
    int inproc[3];
    int onproc[3];
    int irankd[3];
    int orankd[3];

    if (sign == FFTW_FORWARD)
    {
        iaxis     = _iaxis;
        oaxis     = _oaxis;

        for(int id=0;id<3; id++){
            // input parameters
            inloc[id]  = _inloc[id];
            irankd[id] = _irankd[id];
            inproc[id] = _inproc[id];
            // output params
            onloc[id]  = _onloc[id];
            orankd[id] = _orankd[id];
            onproc[id] = _onproc[id];
        }
    }
    else if (sign == FFTW_BACKWARD)
    {
        iaxis     = _oaxis;
        oaxis     = _iaxis;

        for(int id=0;id<3; id++){
            // input parameters
            inloc[id]  = _onloc[id];
            irankd[id] = _orankd[id];
            inproc[id] = _onproc[id];
            // output params
            onloc[id]  = _inloc[id];
            orankd[id] = _irankd[id];
            onproc[id] = _inproc[id];
        }
    }
    else
    {
        UP_CHECK0(false, 'the sign is not FFTW_FORWARD nor FFTW_BACKWARD');
    }

    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    std::memset(count, 0, sizeof(int) * comm_size);

    int dest_rankd[3];

    for (int i2 = 0; i2 < inloc[2]; i2++)
    {
        dest_rankd[2] = (irankd[2] * inloc[2] + i2) / onloc[2];
        for (int i1 = 0; i1 < inloc[1]; i1++)
        {
            dest_rankd[1] = (irankd[1] * inloc[1] + i1) / onloc[1];
            for (int i0 = 0; i0 < inloc[0]; i0++)
            {
                dest_rankd[0] = (irankd[0] * inloc[0] + i0) / onloc[0];

                const int dest_rank = rankindex(dest_rankd, onproc);
                const int buf_idx   = ssend[dest_rank] + count[dest_rank];
                const int id[3]     = {i0, i1, i2};
                const int my_idx    = localindex(id, inloc, iaxis);

                for (int i = 0; i < _nf; i++)
                    sendbuf[buf_idx + i] = v[my_idx + i];

                count[dest_rank] += _nf;
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - Send it */
    //-------------------------------------------------------------------------
    if(sign == FFTW_FORWARD)
        MPI_Alltoallv(sendbuf, nsend, ssend, MPI_DOUBLE, recvbuf, nrecv, srecv, MPI_DOUBLE, MPI_COMM_WORLD);
    else
        MPI_Alltoallv(recvbuf, nrecv, srecv, MPI_DOUBLE, sendbuf, nsend, ssend, MPI_DOUBLE, MPI_COMM_WORLD);
    

    //-------------------------------------------------------------------------
    /** - Fill the memory */
    //-------------------------------------------------------------------------
    std::memset(count, 0, sizeof(int) * comm_size);

    int orig_rankd[3];

    for (int i2 = 0; i2 < onloc[2]; i2++)
    {
        orig_rankd[2] = (orankd[2] * onloc[2] + i2) / inloc[2];
        for (int i1 = 0; i1 < onloc[1]; i1++)
        {
            orig_rankd[1] = (orankd[1] * onloc[1] + i1) / inloc[1];
            for (int i0 = 0; i0 < onloc[0]; i0++)
            {
                orig_rankd[0] = (orankd[0] * onloc[0] + i0) / inloc[0];

                const int orig_rank = rankindex(orig_rankd, inproc);
                const int buf_idx   = srecv[orig_rank] + count[orig_rank];
                const int id[3]     = {i0, i1, i2};
                const int my_idx    = localindex(id, onloc, oaxis);

                for (int i = 0; i < _nf; i++)
                    v[my_idx + i] = recvbuf[buf_idx + i];

                count[orig_rank] += _nf;
            }
        }
    }
}
