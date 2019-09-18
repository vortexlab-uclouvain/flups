/**
 * @file Topology.cpp
 * @author Thomas Gillis & Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-07-26
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "Topology.hpp"

using namespace FLUPS;

/**
 * @brief Construct a new Topology
 * 
 * @param nf the number of field (eg: real = 1, complex = 2)
 * @param axis the dimension of the fastest rotating index (eg: 0, 1 or 2)
 * @param nglob the global size per dim
 * @param nproc the number of proc per dim
 * @param isComplex indicate if the Topology uses complex indexing or not
 * @param axproc gives the order of the rank decomposition (eg. (0,2,1) to start decomposing in X then Z then Y). If NULL is passed, use by default (0,1,2).
 */
Topology::Topology(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3]) {
    BEGIN_FUNC;

    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_comm_size);

    FLUPS_CHECK(nproc[0]*nproc[1]*nproc[2] == _comm_size,"the total number of procs (=%d) have to be = to the comm size (=%d)",nproc[0]*nproc[1]*nproc[2], _comm_size, LOCATION);

    //-------------------------------------------------------------------------
    /** - get proc information  */
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        // total dimension
        _nglob[id] = nglob[id];
        // number of proc in each direction
        _nproc[id] = nproc[id];
        // store the proc axis, used to split the rank
        _axproc[id] = (axproc == NULL) ? id : axproc[id];
    }
    // split the rank
    ranksplit(_rank, _axproc, _nproc, _rankd);

    //-------------------------------------------------------------------------
    /** - get memory axis and complex information  */
    //-------------------------------------------------------------------------
    _axis = axis;
    if (!isComplex) {
        _nf = 1;
    } else {
        _nf = 2;
    }

    //-------------------------------------------------------------------------
    /** - Get the sizes  */
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        // number of unknows everywhere except the last one
        _nbyproc[id] = nglob[id] / nproc[id];  // integer division = floor
        // if we are the last rank in the direction, we take everything what is left
        _nloc[id] = get_nloc(id, _rankd[id], this);
    }
}

/**
 * @brief Destroy the Topology
 * 
 */
Topology::~Topology() {}

/**
 * @brief compute the sarting and ending ids for the current topo in order to be inside the other topology's bounds
 * 
 * @param shift the shift between the 2 topos: current topo in (0,0,0) = other topo in (shift)
 * @param other the other topology
 * @param start the start index to use in the current topo
 * @param end the end index to use in the current topo
 */
void Topology::cmpt_intersect_id(const int shift[3], const Topology* other, int start[3], int end[3]) const {
    BEGIN_FUNC;
    FLUPS_CHECK(this->isComplex() == other->isComplex(), "The two topo have to be both complex or real", LOCATION);

    for (int id = 0; id < 3; id++) {
        const int onglob = other->nglob(id);

        start[id] = 0;
        end[id]   = 0;
        // for the input configuration
        for (int i = 0; i < _nloc[id]; ++i) {
            // get the global id in the other topology
            int oid_global = _rankd[id] * _nbyproc[id] + i + shift[id];
            if (oid_global <= 0) start[id] = i;
            if (oid_global < onglob) end[id] = i + 1;
        }
        FLUPS_CHECK(end[id]-start[id]>0,"iend has to be at least 1 bigger than istart: my nloc = %d %d %d vs other %d %d %d",_nloc[0],_nloc[1],_nloc[2],other->nloc(0),other->nloc(1),other->nloc(2), LOCATION);
    }
}

/**
 * @brief Compute the length (along _axis dimension) of an intersection between two topologies
 * 
 * For each proc of the current topology, we compute the lenght of the intersection with every other proc of the other topology.
 * The length is taken as the maximum continuous amount of data aligned in the current topo's axis
 * 
 * For example, if we want to compute the naxis between proc P0 on the current topo with the horizontal axis,
 * with proc P1 from other topology ( with a vertical axis), the length is 10
 * ```
 *                  P1
 *       ----+-----------------+----
 *           |   |             |
 *  |        |   |             |        |
 *  |        |   l=36          |        |
 *  +-----------------------------------+
 *  |   ---------- l=20 ---------->     |   P0
 *  |                                   |
 *  +--------+--------+--------+--------+
 *  |        |   |             |        |
 *  |        |   v             |        |
 *       ----+-----------------+----
 *              <--- n=10 --->
 * ```
 * 
 * @param other the other topology
 * @param istart the starting index for the current proc data (using #cmpt_intersect_id)
 * @param iend the end index of the current proc data (using #cmpt_intersect_id)
 * @param naxis array of the length from this proc to any other proc in the other topology
 */
void Topology::cmpt_intersect_naxis(const Topology* other, const int istart[3], const int iend[3], const int ishift[3], int* naxis) const {
    int in_idStart[3];
    get_istart_glob(in_idStart, this);
    in_idStart[0] = in_idStart[0] + istart[0];
    in_idStart[1] = in_idStart[1] + istart[1];
    in_idStart[2] = in_idStart[2] + istart[2];

    // go through each other proc
    for (int p2 = 0; p2 < other->nproc(2); p2++) {
        for (int p1 = 0; p1 < other->nproc(1); p1++) {
            for (int p0 = 0; p0 < other->nproc(0); p0++) {
                int p[3] = {p0, p1, p2};

                // compute the intersection volume
                int intersectVol[3];
                for (int id = 0; id < 3; id++) {
                    // left limit of the intersection
                    int in_idleft  = _rankd[id] * _nbyproc[id] + istart[id] + ishift[id];
                    int out_idleft = p[id] * other->nbyproc(id);
                    int limleft    = std::max(in_idleft, out_idleft);
                    // right limit of the intersection
                    int in_idright  = in_idleft + iend[id];
                    int out_idright = out_idleft + get_nloc(id, p[id], other);
                    int limright    = std::min(in_idright, out_idright);
                    // store inside the volume
                    intersectVol[id] = ((limright - limleft) > 0) ? limright - limleft : 0;
                }

                if (intersectVol[0] * intersectVol[1] * intersectVol[2] > 0) {
                    // if we have an intersection, store the length of it
                    naxis[rankindex(p, other)] = intersectVol[_axis];
                } else {
                    // if not, store 0
                    naxis[rankindex(p, other)] = 0;
                }
            }
        }
    }
}

/**
 * @brief display topology most important info
 * 
 */
void Topology::disp() const {
    BEGIN_FUNC;

    FLUPS_INFO("------------------------------------------");
    FLUPS_INFO("## Topology created on proc %d/%d", _rank, _comm_size);
    FLUPS_INFO(" - axis = %d", _axis);
    FLUPS_INFO(" - nglob = %d %d %d", _nglob[0], _nglob[1], _nglob[2]);
    FLUPS_INFO(" - nloc = %d %d %d", _nloc[0], _nloc[1], _nloc[2]);
    FLUPS_INFO(" - nproc = %d %d %d", _nproc[0], _nproc[1], _nproc[2]);
    FLUPS_INFO(" - rankd = %d %d %d", _rankd[0], _rankd[1], _rankd[2]);
    FLUPS_INFO(" - nbyproc = %d %d %d", _nbyproc[0], _nbyproc[1], _nbyproc[2]);
    FLUPS_INFO(" - axproc = %d %d %d", _axproc[0], _axproc[1], _axproc[2]);
    FLUPS_INFO(" - isComplex = %d", _nf == 2);
    // FLUPS_INFO(" - h = %f %f %f",_h[0],_h[1],_h[2]);
    // FLUPS_INFO(" - L = %f %f %f",_L[0],_L[1],_L[2]);
    FLUPS_INFO("------------------------------------------");
}

void Topology::disp_rank() const{
    double* rankdata = (double*) fftw_malloc(sizeof(double)*this->locmemsize());

    for(int i=0; i<this->locmemsize(); i++){
        rankdata[i] = this->rank();
    }

    std::string name = "rank_topo_axis" + std::to_string(this->axis());
    hdf5_dump(this, name, rankdata);
}