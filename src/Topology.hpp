/**
 * @file Topology.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"
#include <cstring>

/**
 * @brief Class Topology
 * 
 * A topology describes the layout of the data on the current processor.
 * 
 * The number of unkowns in each direction owned by a rank divides them in two groups.
 * First, we compute the integer division, nbyproc, between nglob_ and nproc_.
 * 
 * The first group, named g0, owns nbyproc+1 unknowns. The group starts at rank 0 and ends in rank mod(nglob_,nproc_)-1, included.
 * The second group, named g1, owns nbyproc unknowns. The group starts at rank mod(nglob_,nproc_) to rank nproc_, included.
 * 
 */
class Topology {
   protected:
    int       nproc_[3];   /**<@brief number of procs per dim (012-indexing)  */
    int       axproc_[3];  /**<@brief axis of the procs for ranksplit  */
    int       nf_;         /**<@brief the number of doubles inside one unknowns (if complex = 2, if real = 1) */
    int       nloc_[3];    /**<@brief real number of unknowns per dim, for 1 component, local (012-indexing)  */
    int       nmem_[3];    /**<@brief real number of unknowns per dim, for 1 component, local (012-indexing)  */
    int       axis_;       /**<@brief fastest rotating index in the topology  */
    int       rankd_[3];   /**<@brief rank of the current process per dim (012-indexing)  */
    int       nglob_[3];   /**<@brief number of unknowns per dim, global (012-indexing)  */
    int       lda_;        /**<@brief leading dimension of array=the number of components (eg scalar=1, vector=3) */
    // int       nbyproc_[3]; /**<@brief mean number of unknowns per dim = nloc except for the last one (012-indexing)  */
    const int alignment_;
    MPI_Comm  comm_; /**<@brief the comm associated with the topo, with ranks potentially optimized for switchtopos */

    // double h_[3]; //**< @brief grid spacing */
    // double L_[3];//**< @brief length of the domain  */
    // -> We got rid of these, as L changes during a transform occuring in the associated topo, and the computation of h would also
    //      need to depend on the number of points (N, N+2 if we prepare a symmetric transform, etc.)

   public:
    Topology(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm);
    ~Topology();

    /**
     * @name setters
     * 
     * @{
     */
    void change_comm(MPI_Comm comm);
    /**@} */

    /**
     * @name getters
     * 
     * @{
     */
    inline int axis() const { return axis_; }
    inline int lda() const { return lda_; }
    inline int nf() const { return nf_; }
    inline int isComplex() const { return nf_ == 2; }

    inline int nglob(const int dim) const { return nglob_[dim]; }
    inline int nloc(const int dim) const { return nloc_[dim]; }
    inline int nmem(const int dim) const { return nmem_[dim]; }
    inline int nproc(const int dim) const { return nproc_[dim]; }
    inline int rankd(const int dim) const { return rankd_[dim]; }
    // inline int nbyproc(const int dim) const { return nbyproc_[dim]; }
    inline int      axproc(const int dim) const { return axproc_[dim]; }
    inline MPI_Comm get_comm() const { return comm_; }

    /**
     * @brief compute the scalar number of unknowns on each proc, i.e. the number of unkowns for one component
     * 
     * @param id 
     * @return int 
     */
    inline int cmpt_nbyproc(const int id) const {
        return (nglob_[id] / nproc_[id]) + 1 * ((nglob_[id] % nproc_[id]) > rankd_[id]);
        // const int start = (nglob_[id] * rankd_[id]) / nproc_[id];
        // const int end   = (nglob_[id] * (rankd_[id] + 1)) / nproc_[id];
        // return (end - start);
    }

    /**
     * @name Functions to compute the starting index of each topology
     * 
     * @param id the id for one component
     */
    inline int cmpt_start_id(const int id) const {
        return (rankd_[id]) * (nglob_[id] / nproc_[id]) + std::min(rankd_[id], nglob_[id] % nproc_[id]);
        // return (nglob_[id] * rankd_[id]) / nproc_[id];
    }

    /**
     * @brief compute the rank associated to a scalar global id
     * The domain decomposition of the points in ranks (@ref cmpt_start_id) give the following inequality 
     * global_id * nproc_ <= n_glob_*rank_id <= min(global_id +1 , nglob - 1)* nproc
     * Since we use integral division, we only take the right inequality to compute the 
     * rank associated with a global id
     *
     * @param global_id the scalar id of the point considered
     * @param id the direction of interest
     * @return int 
     */
    inline int cmpt_rank_fromid(const int global_id, const int id) const{
        const int nproc_g0 = nglob_[id]%nproc_[id]; // number of procs that have a +1 in their unkowns
        const int nbyproc = nglob_[id]/nproc_[id]; // the number of unknowns in the integer division
        const int global_g0 = nproc_g0*(nbyproc+1); // the number of unknowns in the first group of procs
        return (global_id < global_g0)? global_id/(nbyproc+1) : (global_id-global_g0)/nbyproc + nproc_g0;
        
        // return  (nproc_[id] * std::min(global_id + 1 , nglob_[id] - 1 ))/nglob_[id] ;
    }

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
    void cmpt_sizes();
    void memshift(const int sign,const int lia, double* data);

    /**
     * @brief returns the scalar local size on this proc, i.e. the number of unknowns for one component
     * 
     * @return size_t 
     */
    inline size_t locsize() const { return (size_t)nloc_[0] * (size_t)nloc_[1] * (size_t)nloc_[2] * (size_t)nf_; }

    /**
     * @brief returns the memory size of on this proc for one component
     * 
     * @return size_t 
     */
    inline size_t memdim() const { return (size_t)nmem_[0] * (size_t)nmem_[1] *(size_t) nmem_[2] * (size_t)nf_; }
    /**
     * @brief returns the memory size of on this proc, i.e. the number of dimension * the memory of one dimension
     * 
     * @return size_t 
     */
    inline size_t memsize() const { return (size_t)nmem_[0] * (size_t)nmem_[1] * (size_t)nmem_[2] * (size_t)nf_ * (size_t)lda_; }

    /**
     * @brief returns the starting global index on the current proc
     * 
     */
    inline void get_istart_glob(int istart[3]) const {
        istart[0]   = cmpt_start_id(0);
        istart[1]   = cmpt_start_id(1);
        istart[2]   = cmpt_start_id(2);
    }
    /**@} */

    /**
     * @brief switch the topology to a complex mode
     * 
     */
    inline void switch2complex() {
        if (nf_ == 1) {
            nf_ = 2;
            nglob_[axis_] /= 2;
            nloc_[axis_] /= 2;
            nmem_[axis_] /= 2;
            // nbyproc_[axis_] /= 2;
        }
    }
    /**
     * @brief switch the topology to a real mode
     * 
     */
    inline void switch2real() {
        if (nf_ == 2) {
            nf_ = 1;
            nglob_[axis_] *= 2;
            nloc_[axis_] *= 2;
            nmem_[axis_] *= 2;
            // nbyproc_[axis_] *= 2;
        }
    }

    /**
     * @brief swith the current topo to the scalar case
     * 
     */
    inline void switch2Scalar(){
        lda_=1;
    }
    /**
     * @brief switch the current topo to the vector state
     * 
     */
    inline void switch2Vector(){
        lda_=3;
    }

    void disp() const;
    void disp_rank() ;
};

/**
 * @brief split the rank into rank per dimensions
 * 
 * axproc is not used if comm is of type MPI_CART.
 * 
 * @param rank the rank of the proc (from MPI, in the current communicator of the topo)
 * @param nproc the number of procs along each direction
 * @param comm a communicator. If it is of type MPI_CART, we use that information in the returned rankd.
 * @param rankd the rank per dimension in XYZ format
 */
inline static void ranksplit(const int rank, const int axproc[3], const int nproc[3], MPI_Comm comm, int rankd[3]) {
    const int ax0 = axproc[0];
    const int ax1 = axproc[1];
    const int ax2 = axproc[2];

    int mpi_topo_type;
    MPI_Topo_test(comm, &mpi_topo_type);
    if (mpi_topo_type == MPI_CART) {
        MPI_Cart_coords(comm, rank, 3, rankd);
    } else {
        rankd[ax0] = rank % nproc[ax0];
        rankd[ax1] = (rank % (nproc[ax0] * nproc[ax1])) / nproc[ax0];
        rankd[ax2] = rank / (nproc[ax0] * nproc[ax1]);
    }
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

    int mpi_topo_type, rank;
    MPI_Topo_test(topo->get_comm(), &mpi_topo_type);
    if (mpi_topo_type == MPI_CART) {
        MPI_Cart_rank(topo->get_comm(), rankd, &rank);
    } else {
        rank = rankd[ax0] + topo->nproc(ax0) * (rankd[ax1] + topo->nproc(ax1) * rankd[ax2]);
    }
    return rank;
}

/**
 * @brief compute the memory local index for a point (i0,i1,i2) in axsrc-indexing in a memory in the axtrg-indexing
 * 
 * @param axsrc the FRI for the point (i0,i1,i2)
 * @param i0
 * @param i1 
 * @param i2 
 * @param lia leading index of array
 * @param axtrg the target FRI
 * @param size the size of the memory (012-indexing)
 * @param nf the number of unknows in one element
 * @param lda leading dimension of array, number of vector components
 * @return size_t 
 */
static inline size_t localIndex(const int axsrc, const int i0, const int i1, const int i2,
                                const int axtrg, const int size[3], const int nf, const int lda) {
    const int i[3] = {i0, i1, i2};
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3; 
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;
    const int ax2  = (ax0 + 2) % 3;

    // return localindex_xyz(i[0], i[1], i[2], topo);
    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * (i[dax2] + size[ax2] * lda) );
}
/**
 * @brief compute the memory local index for a point (i0,id) in axsrc-indexing in a memory in the same indexing
 * 
 * The memory id is computed as the collapsed version of the 3 external loops:
 * - the loop on the lda
 * - the loop on the ax2 direction
 * - the loop on the ax1 direction
 * 
 * @param axsrc the FRI for the point (i0,i1,i2)
 * @param i0 the index aligned along the axsrc axis
 * @param id the collapsed id of the outer three loops
 * @param size the size of the memory (012-indexing)
 * @param nf the number of unknows in one element
 * @return size_t 
 */
static inline size_t collapsedIndex(const int axsrc, const int i0, const int id, const int size[3], const int nf) {
    const int ax0  = axsrc;
    return i0 * nf + size[ax0] * nf * id;
}

/**
 * @brief split a global index along the different direction using the FRI axtrg, for one component
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
// static inline void localSplit(const size_t id, const int size[3], const int axtrg, int *id0, int *id1, int *id2, const int nf) {
//     const int ax0   = axtrg;
//     const int ax1   = (ax0 + 1) % 3;
//     const int size0 = (size[ax0] * nf);

//     (*id0) = id % size0;
//     (*id1) = (id % (size0 * size[ax1])) / size0;
//     (*id2) = id / (size0 * size[ax1]);
// }

/**
 * @brief compute the global symmetrized index of a given point, for one component
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
 * 
 * Symmetry computation:
 * We have to take the symmetry around symstart.
 * E.g. in X direction:
 *      `symstart[0] - (ix - symstart[0]) = 2 symstart[0] - ix`
 * In some cases when we have an R2C transform, it ask for 2 additional doubles.
 * The value is meaningless but we would like to avoid segfault and nan's.
 * To do so, we use 2 tricks:
 * - The `abs` is used to stay on the positivie side and hence avoid negative memory access
 * - The `max` is used to prevent the computation of the value in 0, which is never used in the symmetry.
 * 
 * The final formula is then ( in the X direction):
 *      `max( abs(2.0 symstart[0] - ix) , 1)`
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
    is[0] = (symstart[ax0] == 0.0 || ie[dax0] <= symstart[ax0]) ? ie[dax0] : -std::max((int)fabs(2.0 * symstart[ax0] - ie[dax0]), 1);
    is[1] = (symstart[ax1] == 0.0 || ie[dax1] <= symstart[ax1]) ? ie[dax1] : -std::max((int)fabs(2.0 * symstart[ax1] - ie[dax1]), 1);
    is[2] = (symstart[ax2] == 0.0 || ie[dax2] <= symstart[ax2]) ? ie[dax2] : -std::max((int)fabs(2.0 * symstart[ax2] - ie[dax2]), 1);
    
}
#endif