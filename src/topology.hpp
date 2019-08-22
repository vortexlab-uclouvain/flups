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

    bool _isComplex; //**< @brief indicate if the Topology uses complex indexing or not */
    int _nf; //**< @brief the number of doubles inside one unknows (if complex = 2, if real = 1) */
    int _rank; //**< @brief rank of the current process */
    int _comm_size; //**< @brief size of the communicator */
    int _axis;//**< @brief fastest rotating index in the topology  */
    int _nglob[3]; //**< @brief number of unknows per dim, global (always x,y,z)  */
    int _nloc[3];  //**< @brief real number of unknows perd dim, local  */
    int _nproc[3];//**< @brief number of procs per dim  */
    int _rankd[3];//**< @brief rank of the current process per dim  */
    int _nbyproc[3]; //**< @brief mean number of unkows per dim = nloc except for the last one  */

    // double _h[3]; //**< @brief grid spacing */
    // double _L[3];//**< @brief length of the domain  */

public:
    Topology(const int axis, const int nglob[3], const int nproc[3],const bool isComplex);
    ~Topology();

    inline void switch2complex(){
        if(!_isComplex){
            _nf = 2;
            _isComplex = true;
            _nglob[_axis] /= 2;
            _nloc[_axis] /= 2;
            _nbyproc[_axis] /= 2;
        }
    }
    inline void switch2real(){
        if(_isComplex){
            _nf = 1;
            _isComplex = false;
            _nglob[_axis] *= 2;
            _nloc[_axis]  *= 2;
            _nbyproc[_axis] *= 2;
        }
    }

    /**
     * @name getters
     * 
     * @{
     */
    inline int comm_size() const { return _comm_size; }
    inline int axis() const { return _axis; }
    inline int nf() const { return _nf; }
    inline int isComplex() const { return _isComplex; }
    // inline double h(const int dim) const { return _h[dim]; }
    // inline double L(const int dim) const { return _L[dim]; }
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
    inline size_t locmemsize() const{return _nloc[0] * _nloc[1] * _nloc[2] * _nf;}
    inline size_t globmemsize() const{return _nglob[0] * _nglob[1] * _nglob[2] * _nf;}
    void cmpt_intersect_id(const int shift[3], const Topology *other, int start[3], int end[3]) const;
    // inline int  get_idstart_glob(const int dim) const { return _rankd[dim] * _nbyproc[dim]; }
    
    /**
     * @brief compute the rank in the other topology of the data "i"
     * 
     * @param dim 
     * @param other 
     * @param i the current index in my topo that has to match another one
     * @return int 
     */
    inline int cmpt_matchrank(const int dim, const Topology *other, const int i) const
    {
        return std::min((_rankd[dim] * _nbyproc[dim] + i) / other->nbyproc(dim),other->nproc(dim)-1);
    }
    /**@} */
    

    void disp() const;
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
 * @brief return the starting local index for the data (ix,iy,iz) in the order of the dimensions
 * 
 * @param ix index in the X direction
 * @param iy index in the Y direction
 * @param iz index in the Z direction
 * @param topo 
 * @return size_t 
 */
inline static size_t localindex_xyz(const int ix, const int iy, const int iz, const Topology *topo)
{
    const int nf = topo->nf();

    const int i[3] = {ix, iy, iz};
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    return i[ax0] * nf + topo->nloc(ax0) * nf * (i[ax1] + topo->nloc(ax1) * i[ax2]);
}

/**
 * @brief return the starting local index for the data (i0,i1,i2) in the order of the axis
 * 
 * @param ix index in the X direction
 * @param iy index in the Y direction
 * @param iz index in the Z direction
 * @param topo 
 * @return size_t 
 */
inline static size_t localindex_ao(const int i0,const int i1, const int i2, const Topology *topo)
{
    const int nf = topo->nf();
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    
    return i0 * nf + topo->nloc(ax0) * nf * (i1 + topo->nloc(ax1) * i2);
}

inline static void get_idstart_glob(int istart[3], const Topology *topo)
{
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    istart[ax0] = topo->rankd(ax0) * topo->nbyproc(ax0);
    istart[ax1] = topo->rankd(ax1) * topo->nbyproc(ax1);
    istart[ax2] = topo->rankd(ax2) * topo->nbyproc(ax2);
}

#endif