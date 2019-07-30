/**
 * @file topology.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-26
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "topology.hpp"

/**
 * @brief Construct a new Topology
 * 
 * @param nf the number of field (eg: real = 1, complex = 2)
 * @param axis the dimension of the fastest rotating index (eg: 0, 1 or 2)
 * @param nglob the global size per dim
 * @param nproc the number of proc per dim
 */
Topology::Topology(const int axis, const int nglob[3], const int nproc[3],const double h[3])
{

    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_comm_size);

    _axis = axis;

    //-------------------------------------------------------------------------
    /** - get the rankd for input and output  */
    //-------------------------------------------------------------------------
    ranksplit(_rank, nproc, _rankd);

    for (int id = 0; id < 3; id++)
    {
        // total dimension
        _nglob[id] = nglob[id];
        // lenght of the domain and h
        _h[id] = h[id];
        _L[id] = h[id]*nglob[id];
        // number of proc in each direction
        _nproc[id] = nproc[id];
        // number of unknows everywhere except the last one
        _nbyproc[id] = nglob[id] / nproc[id]; // integer division = floor
        // if we are the last rank in the direction, we take everything what is left
        _nloc[id] = (_rankd[id] != (nproc[id]-1)) ? _nbyproc[id] : std::max(_nbyproc[id],nglob[id]-_nbyproc[id]*_rankd[id]);
    }
}

/**
 * @brief Destroy the Topology
 * 
 */
Topology::~Topology(){}


/**
 * @brief compute the sarting and ending ids to use in order to fit inside the other topology's bounds
 * 
 * @param shift the shift between the 2 topos: current topo in (#_shift) = other topo in (0,0,0)
 * @param other the other topology
 * @param start the start index to use in the current topo
 * @param end the end index to use in the current topo
 */
void Topology::cmpt_intersect_id(const int shift[3], const Topology* other,int start[3],int end[3]) const
{
    for (int id = 0; id < 3; id++)
    {
        const int onglob = other->nglob(id);

        start[id] = 0;
        // for the input configuration
        for (int i = 0; i < _nloc[id]; ++i)
        {
            // get the global id in the other topology
            int oid_global = _rankd[id] * _nbyproc[id] + i + shift[id];
            if (oid_global <= 0) start[id] = i;
            if (oid_global < onglob) end[id] = i + 1;
        }
    }
}