/**
 * @file topology.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-26
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "defines.hpp"
#include "mpi.h"

/**
 * @brief Class Topology
 * 
 */
class Topology
{

protected:
    int _rank; //**< @brief rank of the current process */
    int _comm_size; //**< @brief size of the communicator */
    int _axis;//**< @brief fastest rotating index in the topology  */
    int _nglob[3]; //**< @brief number of unknows per dim, global (always x,y,z)  */
    int _nloc[3];  //**< @brief real number of unknows perd dim, local  */
    int _nproc[3];//**< @brief number of procs per dim  */
    int _rankd[3];//**< @brief rank of the current process per dim  */
    int _nbyproc[3]; //**< @brief mean number of unkows per dim = nloc except for the last one  */

    double _h[3]; //**< @brief grid spacing */
    double _L[3];//**< @brief length of the domain  */

public:
    Topology(const int axis, const int nglob[3], const int nproc[3], const double h[3]);
    ~Topology();

    /**
     * @name getters
     * 
     * @{
     */
    inline int comm_size() const { return _comm_size; }
    inline int axis() const { return _axis; }
    inline int h(const int dim) const { return _h[dim]; }
    inline int L(const int dim) const { return _L[dim]; }
    inline int nglob(const int dim) const { return _nglob[dim]; }
    inline int nloc(const int dim) const { return _nloc[dim]; }
    inline int nproc(const int dim) const { return _nproc[dim]; }
    inline int rankd(const int dim) const { return _rankd[dim]; }
    inline int nbyproc(const int dim) const { return _nbyproc[dim]; }
    /**@} */


    /**
     * @name Usefull functions manipulating indexes
     * 
     * @{
     */
    inline size_t locsize() const{return _nloc[0] * _nloc[1] * _nloc[2];}
    void cmpt_intersect_id(const int shift[3], const Topology *other, int start[3], int end[3]) const;
    inline int  get_idstart_glob(const int dim) const { return _rankd[dim] * _nbyproc[dim]; }
    inline void get_idstart_glob(int istart[3]) const
    {
        istart[0] = _rankd[0] * _nbyproc[0];
        istart[1] = _rankd[1] * _nbyproc[1];
        istart[2] = _rankd[2] * _nbyproc[2];
    }
    /**
     * @brief compute the rank in the other topology of the data "i"
     * 
     * @param dim 
     * @param other 
     * @param i 
     * @return int 
     */
    inline int cmpt_matchrank(const int dim, const Topology *other, const int i) const
    {
        return std::min((_rankd[dim] * _nbyproc[dim] + i) / other->nbyproc(dim),other->nproc(dim)-1);
    }
    /**@} */
    
};

/**
 * @brief split the rank into rank per dimensions
 * 
 * @param rank 
 * @param nproc 
 * @param rankd 
 */
inline static void ranksplit(const int rank, const int nproc[3], int rankd[3])
{
    rankd[0] = rank % nproc[0];
    rankd[1] = (rank % (nproc[0] * nproc[1])) / nproc[0];
    rankd[2] = rank / (nproc[0] * nproc[1]);
}

/**
 * @brief get the rank from the rank per dimension
 * 
 * @param rankd 
 * @param topo 
 * @return int 
 */
inline static int rankindex(const int rankd[3], const Topology *topo)
{
    return rankd[0] + topo->nproc(0) * (rankd[1] + topo->nproc(1) * rankd[2]);
}

/**
 * @brief return the local index for the data "i", given the #_axis and #_nloc
 * 
 * @param i 
 * @param topo 
 * @return size_t 
 */
inline static size_t localindex(const int i0,const int i1,const int i2, const Topology *topo)
{
    const int i[3] = {i0,i1,i2};
    const int axis  = topo->axis();
    const int axis1 = (topo->axis() + 1) % 3;
    const int axis2 = (topo->axis() + 2) % 3;
    return i[axis] + topo->nloc(axis) * (i[axis1] + topo->nloc(axis1) * i[axis2]);
}

#endif