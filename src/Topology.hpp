/**
 * @file Topology.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"

/**
 * @brief Class Topology
 * 
 * A topology describes the layout of the data on the current processor.
 * 
 */
class Topology {
   protected:
    int      _nproc[3];   /**<@brief number of procs per dim (012-indexing)  */
    int      _axproc[3];  /**<@brief axis of the procs for ranksplit  */
    int      _nf;         /**<@brief the number of doubles inside one unknows (if complex = 2, if real = 1) */
    int      _nloc[3];    /**<@brief real number of unknows perd dim, local (012-indexing)  */
    int      _nmem[3];    /**<@brief real number of unknows perd dim, local (012-indexing)  */
    int      _axis;       /**<@brief fastest rotating index in the topology  */
    int      _rankd[3];   /**<@brief rank of the current process per dim (012-indexing)  */
    int      _nglob[3];   /**<@brief number of unknows per dim, global (012-indexing)  */
    int      _nbyproc[3]; /**<@brief mean number of unkows per dim = nloc except for the last one (012-indexing)  */
    MPI_Comm _comm;       /**<@brief the comm associated with the topo, with ranks potentially optimized for switchtopos */

    // double _h[3]; //**< @brief grid spacing */
    // double _L[3];//**< @brief length of the domain  */
    // -> We got rid of these, as L changes during a transform occuring in the associated topo, and the computation of h would also
    //      need to depend on the number of points (N, N+2 if we prepare a symmetric transform, etc.)

   public:
    // Topology(const int axis, const int nglob[3], const int nproc[3], const bool isComplex);
    Topology(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment);
    ~Topology();

    /**
     * @name setters
     * 
     * @{
     */
    void set_comm(MPI_Comm comm);
    /**@} */

    /**
     * @name getters
     * 
     * @{
     */
    inline int axis() const { return _axis; }
    inline int nf() const { return _nf; }
    inline int isComplex() const { return _nf == 2; }

    inline int nglob(const int dim) const { return _nglob[dim]; }
    inline int nloc(const int dim) const { return _nloc[dim]; }
    inline int nmem(const int dim) const { return _nmem[dim]; }
    inline int nproc(const int dim) const { return _nproc[dim]; }
    inline int rankd(const int dim) const { return _rankd[dim]; }
    inline int nbyproc(const int dim) const { return _nbyproc[dim]; }
    inline int axproc(const int dim) const { return _axproc[dim]; }
    inline MPI_Comm get_comm() const {return _comm; }
    /**
     * @name Functions to compute intersection data with other Topologies
     * 
     * @{
     */
    void cmpt_intersect_id(const int shift[3], const Topology *other, int start[3], int end[3]) const;
    /**@} */

    /**
     * @name Usefull functions manipulating indexes and memory
     * 
     * @{
     */
    /**
     * @brief returns the local size of on this proc
     * 
     * @return size_t 
     */
    inline size_t locsize() const { return (size_t)(_nloc[0] * _nloc[1] * _nloc[2] * _nf); }
    /**
     * @brief returns the memory size of on this proc
     * 
     * @return size_t 
     */
    inline size_t memsize() const { return (size_t)(_nmem[0] * _nmem[1] * _nmem[2] * _nf); }

    /**
     * @brief returns the starting global index of the current proc
     * 
     */
    inline void get_istart_glob(int istart[3]) const {
        istart[0]   = _rankd[0] * _nbyproc[0];
        istart[1]   = _rankd[1] * _nbyproc[1];
        istart[2]   = _rankd[2] * _nbyproc[2];
    }

    /**
     * @brief switch the topology to a complex mode
     * 
     */
    inline void switch2complex() {
        if (_nf == 1) {
            _nf = 2;
            _nglob[_axis] /= 2;
            _nloc[_axis] /= 2;
            _nmem[_axis] /= 2;
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
            _nmem[_axis] *= 2;
            _nbyproc[_axis] *= 2;
        }
    }

    void disp() const;
    void disp_rank() ;
};

/**
 * @brief split the rank into rank per dimensions
 * 
 * @param rank the rank of the proc (from MPI, in the current communicator of the topo)
 * @param nproc the number of procs along each direction
 * @param rankd the rank per dimension in XYZ format
 */
inline static void ranksplit(const int rank, const int axproc[3], const int nproc[3], int rankd[3]) {
    const int ax0 = axproc[0];
    const int ax1 = axproc[1];
    const int ax2 = axproc[2];
    rankd[ax0]    = rank % nproc[ax0];
    rankd[ax1]    = (rank % (nproc[ax0] * nproc[ax1])) / nproc[ax0];
    rankd[ax2]    = rank / (nproc[ax0] * nproc[ax1]);
}

/**
 * @brief get the rank from the rank per dimension
 * 
 * @param rankd the rank in XYZ format
 * @param topo the topology
 * @return int 
 */
inline static int rankindex(const int rankd[3], const Topology *topo) {
    const int ax0 = topo->axproc(0);
    const int ax1 = topo->axproc(1);
    const int ax2 = topo->axproc(2);
    return rankd[ax0] + topo->nproc(ax0) * (rankd[ax1] + topo->nproc(ax1) * rankd[ax2]);
}

// /**
//  * @brief return the starting local index for the data (ix,iy,iz) in the order of the dimensions
//  * 
//  * @param ix index in the X direction
//  * @param iy index in the Y direction
//  * @param iz index in the Z direction
//  * @param topo 
//  * @return size_t 
//  */
// inline static size_t localindex_xyz(const int ix, const int iy, const int iz, const Topology *topo) {
//     const int nf = topo->nf();

//     const int i[3] = {ix, iy, iz};
//     const int ax0  = topo->axis();
//     const int ax1  = (ax0 + 1) % 3;
//     const int ax2  = (ax0 + 2) % 3;

//     return i[ax0] * nf + topo->nloc(ax0) * nf * (i[ax1] + topo->nloc(ax1) * i[ax2]);
// }

// /**
//  * @brief return the local index in memory for the data (i0,i1,i2) in the order of the axis, and in double indexing
//  * 
//  * @param i0 index along the ax0 direction (the fast rotating index in the current topo)
//  * @param i1 index along the ax1 direction
//  * @param i2 index along the ax2 direction
//  * @param topo 
//  * @return size_t 
//  */
// inline static size_t localindex_ao(const int i0, const int i1, const int i2, const Topology *topo) {
//     const int nf  = topo->nf();
//     const int ax0 = topo->axis();
//     const int ax1 = (ax0 + 1) % 3;

//     return i0 * nf + topo->nloc(ax0) * nf * (i1 + topo->nloc(ax1) * i2);
// }
// /**
//  * @brief return the starting local index for the data (i0,i1,i2) in the order of the axis given
//  *
//  * @param axis index of the axis corresponding to i0
//  * @param ix index in the X direction
//  * @param iy index in the Y direction
//  * @param iz index in the Z direction
//  * @param topo 
//  * @return size_t 
//  */
// inline static size_t localindex(const int axis, const int i0, const int i1, const int i2, const Topology *topo) {
//     const int nf   = topo->nf();
//     const int i[3] = {i0, i1, i2};
//     // compute the shift to perform from the axis reference to
//     const int dax0 = (3 + topo->axis() - axis) % 3;
//     const int dax1 = (dax0 + 1) % 3;
//     const int dax2 = (dax0 + 2) % 3;

//     const int ax0 = topo->axis();
//     const int ax1 = (ax0 + 1) % 3;

//     // return localindex_xyz(i[0], i[1], i[2], topo);
//     return i[dax0] * nf + topo->nloc(ax0) * nf * (i[dax1] + topo->nloc(ax1) * i[dax2]);
// }
/**
 * @brief compute the memory local index for a point (i0,i1,i2) in axsrc-indexing in a memory in the axtrg-indexing
 * 
 * @param axsrc the FRI for the point (i0,i1,i2)
 * @param i0
 * @param i1 
 * @param i2 
 * @param axtrg the target FRI
 * @param size the size of the memory (012-indexing)
 * @param nf the number of unknows in one element
 * @return size_t 
 */
static inline size_t localIndex(const int axsrc, const int i0, const int i1, const int i2,
                                const int axtrg, const int size[3], const int nf) {
    const int i[3] = {i0, i1, i2};
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;

    // return localindex_xyz(i[0], i[1], i[2], topo);
    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * i[dax2]);
}
/**
 * @brief compute the memory local index for a point (i0,id) in axsrc-indexing in a memory in the same indexing
 * 
 * The memory id is computed as the collapsed version of the 2 external loops
 * 
 * @param axsrc the FRI for the point (i0,i1,i2)
 * @param i0 the index aligned along the axsrc axis
 * @param id the collapsed id of the outer two loops
 * @param size the size of the memory (012-indexing)
 * @param nf the number of unknows in one element
 * @return size_t 
 */
static inline size_t collapsedIndex(const int axsrc, const int i0, const int id, const int size[3], const int nf) {
    const int ax0  = axsrc;
    return i0 * nf + size[ax0] * nf * id;
}

/**
 * @brief split a global index along the different direction using the FRI axtrg
 * 
 * @param id the global id
 * @param size the size in 012-indexing
 * @param axtrg the target axis
 * @param idv the indexes along each directions
 */
static inline void localSplit(const size_t id, const int size[3], const int axtrg, int idv[3], const int nf) {
    const int ax0   = axtrg;
    const int ax1   = (ax0 + 1) % 3;
    const int ax2   = (ax0 + 2) % 3;
    const int size0 = (size[ax0] * nf);

    idv[ax0] = id % size0;
    idv[ax1] = (id % (size0 * size[ax1])) / size0;
    idv[ax2] = id / (size0 * size[ax1]);
}
static inline void localSplit(const size_t id, const int size[3], const int axtrg, int *id0, int *id1, int *id2, const int nf) {
    const int ax0   = axtrg;
    const int ax1   = (ax0 + 1) % 3;
    const int size0 = (size[ax0] * nf);

    (*id0) = id % size0;
    (*id1) = (id % (size0 * size[ax1])) / size0;
    (*id2) = id / (size0 * size[ax1]);
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
 * @brief compute the global symmetrized index of a given point.
 * 
 * The 3 indexes are given along axsrc axis (i0,i1,i2) while the symmetrized output is given along the axtrg axis.
 * 
 * @param axsrc 
 * @param i0 
 * @param i1 
 * @param i2 
 * @param istart 
 * @param symstart 
 * @param axtrg 
 * @param is 
 */
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