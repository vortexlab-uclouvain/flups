/**
 * @file topology.hpp
 * @author Thomas Gillis & Denis-Gabriel Caprace
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
 * A topology describes the layout of the data on the current processor.
 * 
 */
class Topology {
   protected:
    int _nf;         /**<@brief the number of doubles inside one unknows (if complex = 2, if real = 1) */
    int _rank;       /**<@brief rank of the current process */
    int _comm_size;  /**<@brief size of the communicator */
    int _axis;       /**<@brief fastest rotating index in the topology  */
    int _nglob[3];   /**<@brief number of unknows per dim, global (always XYZ)  */
    int _nloc[3];    /**<@brief real number of unknows perd dim, local  */
    int _nproc[3];   /**<@brief number of procs per dim  */
    int _rankd[3];   /**<@brief rank of the current process per dim  */
    int _nbyproc[3]; /**<@brief mean number of unkows per dim = nloc except for the last one  */

    // double _h[3]; //**< @brief grid spacing */
    // double _L[3];//**< @brief length of the domain  */
    // -> We got rid of these, as L changes during a transform occuring in the associated topo, and the computation of h would also
    //      need to depend on the number of points (N, N+2 if we prepare a symmetric transform, etc.)

   public:
    Topology(const int axis, const int nglob[3], const int nproc[3], const bool isComplex);
    ~Topology();

    /**
     * @name getters
     * 
     * @{
     */
    inline int comm_size() const { return _comm_size; }
    inline int axis() const { return _axis; }
    inline int nf() const { return _nf; }
    inline int isComplex() const { return _nf == 2; }
    // inline double h(const int dim) const { return _h[dim]; }
    // inline double L(const int dim) const { return _L[dim]; }
    inline int nglob(const int dim) const { return _nglob[dim]; }
    inline int nloc(const int dim) const { return _nloc[dim]; }
    inline int nproc(const int dim) const { return _nproc[dim]; }
    inline int rankd(const int dim) const { return _rankd[dim]; }
    inline int nbyproc(const int dim) const { return _nbyproc[dim]; }
    /**@} */

    /**
     * @name Functions to compute intersection data with other Topologies
     * 
     * @{
     */
    void cmpt_intersect_id(const int shift[3], const Topology *other, int start[3], int end[3]) const;
    void cmpt_intersect_naxis(const Topology *other, const int istart[3], const int iend[3], const int ishift[3], int *naxis) const;
    /**@} */

    /**
     * @name Usefull functions manipulating indexes and memory
     * 
     * @{
     */
    /**
     * @brief returns the local size of the memory (in double!) on this proc
     * 
     * @return size_t 
     */
    inline size_t locmemsize() const { return _nloc[0] * _nloc[1] * _nloc[2] * _nf; }
    /**
     * @brief returns the global memory size (in double!)
     * 
     * @return size_t 
     */
    inline size_t globmemsize() const { return _nglob[0] * _nglob[1] * _nglob[2] * _nf; }

    /**
     * @brief compute the rank in the other topology of the processor containing the data i
     * 
     * @param dim the dimension along which i is measured in the (indexing of the other topology!)
     * @param other the other topo
     * @param i the global index in the current topo
     * @return int 
     */
    inline int cmpt_matchrank(const int dim, const Topology *other, const int i) const {
        return std::min((_rankd[dim] * _nbyproc[dim] + i) / other->nbyproc(dim), other->nproc(dim) - 1);
    }
    /**@} */

    /**
     * @brief switch the topology to a complex mode
     * 
     */
    inline void switch2complex() {
        if (_nf == 1) {
            _nf = 2;
            _nglob[_axis] /= 2;
            _nloc[_axis] /= 2;
            _nbyproc[_axis] /= 2;
        }
    }
    /**
     * @brief switch the topology to a real mode
     * 
     */
    inline void switch2real() {
        if (_nf == 2) {
            _nf = 1;
            _nglob[_axis] *= 2;
            _nloc[_axis] *= 2;
            _nbyproc[_axis] *= 2;
        }
    }

    void disp() const;
};

/**
 * @brief split the rank into rank per dimensions
 * 
 * @param rank the rank of the proc (from MPI)
 * @param nproc the number of procs along each direction
 * @param rankd the rank per dimension in XYZ format
 */
inline static void ranksplit(const int rank, const int nproc[3], int rankd[3]) {
    rankd[0] = rank % nproc[0];
    rankd[1] = (rank % (nproc[0] * nproc[1])) / nproc[0];
    rankd[2] = rank / (nproc[0] * nproc[1]);
}

/**
 * @brief get the rank from the rank per dimension
 * 
 * @param rankd the rank in XYZ format
 * @param topo the topology
 * @return int 
 */
inline static int rankindex(const int rankd[3], const Topology *topo) {
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
inline static size_t localindex_xyz(const int ix, const int iy, const int iz, const Topology *topo) {
    const int nf = topo->nf();

    const int i[3] = {ix, iy, iz};
    const int ax0  = topo->axis();
    const int ax1  = (ax0 + 1) % 3;
    const int ax2  = (ax0 + 2) % 3;

    return i[ax0] * nf + topo->nloc(ax0) * nf * (i[ax1] + topo->nloc(ax1) * i[ax2]);
}

/**
 * @brief return the local index in memory for the data (i0,i1,i2) in the order of the axis, and in double indexing
 * 
 * @param i0 index along the ax0 direction (the fast rotating index in the current topo)
 * @param i1 index along the ax1 direction
 * @param i2 index along the ax2 direction
 * @param topo 
 * @return size_t 
 */
inline static size_t localindex_ao(const int i0, const int i1, const int i2, const Topology *topo) {
    const int nf  = topo->nf();
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;

    return i0 * nf + topo->nloc(ax0) * nf * (i1 + topo->nloc(ax1) * i2);
}

/**
 * @brief return the starting local index for the data (i0,i1,i2) in the order of the axis given
 *
 * @param axis index of the axis corresponding to i0
 * @param ix index in the X direction
 * @param iy index in the Y direction
 * @param iz index in the Z direction
 * @param topo 
 * @return size_t 
 */
inline static size_t localindex(const int axis, const int i0, const int i1, const int i2, const Topology *topo) {
    const int nf   = topo->nf();
    const int i[3] = {i0, i1, i2};
    // compute the shift to perform from the axis reference to
    const int dax0 = (3 + topo->axis() - axis) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;

    // return localindex_xyz(i[0], i[1], i[2], topo);
    return i[dax0] * nf + topo->nloc(ax0) * nf * (i[dax1] + topo->nloc(ax1) * i[dax2]);
}

/**
 * @brief Get the istart in global indexing
 * 
 * @param istart start index along the ax0 direction (fast rotating index in current topo), ax1 and ax2
 * @param topo 
 */
inline static void get_istart_glob(int istart[3], const Topology *topo) {
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    istart[ax0] = topo->rankd(ax0) * topo->nbyproc(ax0);
    istart[ax1] = topo->rankd(ax1) * topo->nbyproc(ax1);
    istart[ax2] = topo->rankd(ax2) * topo->nbyproc(ax2);
}

/**
 * @brief return the number of local point for the proc index iproc in the dimension id 
 * 
 * @param id the dimension ID
 * @param iproc the id of the proc in the direction id
 * @param topo the topology
 */
inline static int get_nloc(const int id, const int iproc, const Topology *topo) {
    return (iproc != (topo->nproc(id) - 1)) ? topo->nbyproc(id) : std::max(topo->nbyproc(id), topo->nglob(id) - topo->nbyproc(id) * iproc);
}

inline static void cmpt_symID(const int axsrc, const int i0, const int i1, const int i2, const int istart[3], const double symstart[3], const int axtrg, int is[3]) {
    // get the global indexes in the axsrc configuration
    const int ie[3] = {(istart[axsrc] + i0), (istart[(axsrc + 1) % 3] + i1), (istart[(axsrc + 2) % 3] + i2)};
    // cmpt the shift in axis and the axis for the symstart
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;
    const int ax2  = (ax0 + 2) % 3;
    // fill the array in the axtrg configuration
    is[0] = (symstart[ax0] == 0.0 || ie[dax0] <= symstart[ax0]) ? ie[dax0] : std::max((int)fabs(2.0 * symstart[ax0] - ie[dax0]), 1);
    is[1] = (symstart[ax1] == 0.0 || ie[dax1] <= symstart[ax1]) ? ie[dax1] : std::max((int)fabs(2.0 * symstart[ax1] - ie[dax1]), 1);
    is[2] = (symstart[ax2] == 0.0 || ie[dax2] <= symstart[ax2]) ? ie[dax2] : std::max((int)fabs(2.0 * symstart[ax2] - ie[dax2]), 1);

}
#endif