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
inline static void ranksplit(const int rank, const int nproc[3], int rankd[3])
{
    rankd[0] = rank % nproc[0];
    rankd[1] = (rank % (nproc[0] * nproc[1])) / nproc[0];
    rankd[2] = rank / (nproc[0] * nproc[1]);
}

Reorder_MPI::Reorder_MPI(const int nglob[3], const int nf, const int inproc[3], int axis0, int onproc[3], int axis1)
{
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //-------------------------------------------------------------------------
    /** - get the rankd for input and output  */
    //-------------------------------------------------------------------------
    ranksplit(rank, inproc, _irankd);
    ranksplit(rank, onproc, _orankd);

    _nf = nf;
    _axis0 = axis0;
    _axis1 = axis1;

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
    /** - allocate the buffer */
    //-------------------------------------------------------------------------
    _nsend = (int *)calloc(comm_size, sizeof(int));
    _nrecv = (int *)calloc(comm_size, sizeof(int));
    _ssend = (int *)calloc(comm_size, sizeof(int));
    _srecv = (int *)calloc(comm_size, sizeof(int));
    _count = (int *)calloc(comm_size, sizeof(int));

    _bufsend = (double *)fftw_malloc(sizeof(double) * nf * _inloc[0] * _inloc[1] * _inloc[2]);
    _bufrecv = (double *)fftw_malloc(sizeof(double) * nf * _onloc[0] * _onloc[1] * _onloc[2]);

    printf("allocating value of %d %d %d\n",_inloc[0],_inloc[1],_inloc[2]);
    printf("allocating value of %d %d %d\n",_onloc[0],_onloc[1],_onloc[2]);

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
                _nsend[rankindex(dest_rankd, _onproc)] += nf;
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
        free(_nsend);
    if (_nrecv != NULL)
        free(_nrecv);
    if (_ssend != NULL)
        free(_ssend);
    if (_srecv != NULL)
        free(_srecv);

    if (_bufsend != NULL)
        fftw_free(_bufsend);
    if (_bufrecv != NULL)
        fftw_free(_bufrecv);
}

void Reorder_MPI::execute(double *v)
{
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //-------------------------------------------------------------------------
    /** - fill the buffers */
    //-------------------------------------------------------------------------
    int dest_rankd[3];
    for(int ip=0; ip<comm_size; ip++) _count[ip] = 0;
    for (int i2 = 0; i2 < _inloc[2]; ++i2)
    {
        dest_rankd[2] = (_irankd[2] * _inloc[2] + i2) / _onloc[2];
        for (int i1 = 0; i1 < _inloc[1]; ++i1)
        {
            dest_rankd[1] = (_irankd[1] * _inloc[1] + i1) / _onloc[1];
            for (int i0 = 0; i0 < _inloc[0]; ++i0)
            {
                dest_rankd[0] = (_irankd[0] * _inloc[0] + i0) / _onloc[0];
                const int dest_rank = rankindex(dest_rankd, _onproc);
                const int id[3] = {i0, i1, i2};

                for (int i = 0; i < _nf; i++)
                {
                    _bufsend[_ssend[dest_rank] + (_count[dest_rank]++)] = v[localindex(id, _inloc, _axis0)];
                }
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - Send it */
    //-------------------------------------------------------------------------
    MPI_Alltoallv(_bufsend, _nsend, _ssend, MPI_DOUBLE, _bufrecv, _nrecv, _srecv, MPI_DOUBLE, MPI_COMM_WORLD);

    //-------------------------------------------------------------------------
    /** - Fill the memory */
    //-------------------------------------------------------------------------
    int orig_rankd[3];
    for(int ip=0; ip<comm_size; ip++) _count[ip] = 0;
    for (int i2 = 0; i2 < _onloc[2]; i2++)
    {
        orig_rankd[2] = (_orankd[2] * _onloc[2] + i2) / _inloc[2];
        for (int i1 = 0; i1 < _onloc[1]; i1++)
        {
            orig_rankd[1] = (_orankd[1] * _onloc[1] + i1) / _inloc[1];
            for (int i0 = 0; i0 < _onloc[0]; i0++)
            {
                orig_rankd[0] = (_orankd[0] * _onloc[0] + i0) / _inloc[0];
                int origrank = rankindex(orig_rankd, _inproc);
                const int id[3] = {i0, i1, i2};

                for (int i = 0; i < _nf; i++)
                {
                    v[localindex(id, _onloc, _axis1)] = _bufrecv[_srecv[origrank] + (_count[origrank]++)];   
                }
            }
        }
    }
}
