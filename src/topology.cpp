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
 * @param isComplex indicate if the Topology uses complex indexing or not
 */
Topology::Topology(const int axis, const int nglob[3], const int nproc[3], const bool isComplex) {
    BEGIN_FUNC

    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_comm_size);

    //-------------------------------------------------------------------------
    /** - get memory axis and complex information  */
    //-------------------------------------------------------------------------
    _axis      = axis;
    _isComplex = isComplex;
    if (!_isComplex)
        _nf = 1;
    else
        _nf = 2;

    //-------------------------------------------------------------------------
    /** - get the rankd for input and output  */
    //-------------------------------------------------------------------------
    ranksplit(_rank, nproc, _rankd);

    for (int id = 0; id < 3; id++) {
        // total dimension
        _nglob[id] = nglob[id];
        // number of proc in each direction
        _nproc[id] = nproc[id];
        // number of unknows everywhere except the last one
        _nbyproc[id] = nglob[id] / nproc[id];  // integer division = floor
        // if we are the last rank in the direction, we take everything what is left
        _nloc[id] = (_rankd[id] != (nproc[id] - 1)) ? _nbyproc[id] : std::max(_nbyproc[id], nglob[id] - _nbyproc[id] * _rankd[id]);
    }
}

/**
 * @brief Destroy the Topology
 * 
 */
Topology::~Topology() {}

/**
 * @brief compute the sarting and ending ids to use in order to fit inside the other topology's bounds
 * 
 * @param shift the shift between the 2 topos: current topo in (#_shift) = other topo in (0,0,0)
 * @param other the other topology
 * @param start the start index to use in the current topo
 * @param end the end index to use in the current topo
 */
void Topology::cmpt_intersect_id(const int shift[3], const Topology* other, int start[3], int end[3]) const {
    BEGIN_FUNC
    UP_CHECK0(_isComplex == other->isComplex(), "The two topo have to be both complex or real");

    for (int id = 0; id < 3; id++) {
        const int onglob = other->nglob(id);

        start[id] = 0;
        end[id]   = 1;
        // for the input configuration
        for (int i = 0; i < _nloc[id]; ++i) {
            // get the global id in the other topology
            int oid_global = _rankd[id] * _nbyproc[id] + i + shift[id];
            if (oid_global <= 0) start[id] = i;
            if (oid_global < onglob) end[id] = i + 1;
        }
    }

    // printf("Intersection from %d %d %d to %d %d %d\n",start[0],start[1],start[2],end[0],end[1],end[2]);
}

void Topology::disp() const {
    BEGIN_FUNC

    INFO("------------------------------------------\n");
    INFO3("## Topology created on proc %d/%d", _rank, _comm_size);
    INFO2(" - axis = %d\n", _axis);
    INFO4(" - nglob = %d %d %d\n", _nglob[0], _nglob[1], _nglob[2]);
    INFO4(" - nloc = %d %d %d\n", _nloc[0], _nloc[1], _nloc[2]);
    INFO4(" - nproc = %d %d %d\n", _nproc[0], _nproc[1], _nproc[2]);
    INFO4(" - rankd = %d %d %d\n", _rankd[0], _rankd[1], _rankd[2]);
    INFO4(" - nbyproc = %d %d %d\n", _nbyproc[0], _nbyproc[1], _nbyproc[2]);
    INFO2(" - isComplex = %d\n", _isComplex);
    // INFO4(" - h = %f %f %f\n",_h[0],_h[1],_h[2]);
    // INFO4(" - L = %f %f %f\n",_L[0],_L[1],_L[2]);
    INFO("------------------------------------------\n");
}