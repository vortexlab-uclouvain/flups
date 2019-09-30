/**
 * @file SwitchTopo.hpp
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

#ifndef SWITCHTOPO_HPP
#define SWITCHTOPO_HPP

#include <cstring>
#include "Topology.hpp"
#include "defines.hpp"
#include "mpi.h"
#include "omp.h"
#include "Profiler.hpp"

/**
 * @brief Defines the basic interface for the SwitchTopo objects.
 * 
 * A SwitchTopo reorganizes the memory between 2 different topologies, also accounting for a
 * "principal axis" which is aligned with the fast rotating index.
 * 
 * Communications are handled by packing data in blocks (i.e. chuncks). These 
 * smaller pieces are sent and recieved in a way which is left to the implementation choice.
 * 
 * A memory shift can be specified in the switch between the topologies, when
 * there is a need for skipping some data at the left or right side of a given direction.
 * 
 */
class FLUPS::SwitchTopo {
   protected:
    MPI_Comm _subcomm; /**<@brief the subcomm for this switchTopo */
    int _exSize[3]; /**<@brief exchanged size in each dimension (012-indexing) */

    int _nByBlock[3]; /**<@brief The number of data per blocks in each dim (!same on each process! and 012-indexing)  */
    int _istart[3]; /**<@brief the starting index for #_topo_in to be inside #_topo_out  */
    int _ostart[3]; /**<@brief the starting index for #_topo_out to be inside #_topo_in  */
    int _iend[3];   /**<@brief the ending index for #_topo_in to be inside #_topo_out  */
    int _oend[3];   /**<@brief the ending index for #_topo_out to be inside #_topo_in  */

    int _inBlock[3];  /**<@brief the local number of block in each dim in the input topology */
    int _onBlock[3];  /**<@brief the local number of block in each dim in the output topology  */

    int* _iBlockSize[3] = {NULL,NULL,NULL}; /**<@brief The number of data per blocks in each dim for each block (!same on each process! and 012-indexing)  */
    int* _oBlockSize[3] = {NULL,NULL,NULL}; /**<@brief The number of data per blocks in each dim for each block (!same on each process! and 012-indexing)  */

    opt_int_ptr _i2o_destRank = NULL; /**<@brief The destination rank in the output topo of each block */
    opt_int_ptr _o2i_destRank = NULL; /**<@brief The destination rank in the output topo of each block */

    const Topology *_topo_in  = NULL; /**<@brief input topology  */
    const Topology *_topo_out = NULL; /**<@brief  output topology */

    opt_double_ptr *_sendBuf = NULL; /**<@brief The send buffer for MPI send */
    opt_double_ptr *_recvBuf = NULL; /**<@brief The recv buffer for MPI recv */

    Profiler* _prof = NULL;

   public:
    virtual void          setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData) = 0;
    virtual void          execute(opt_double_ptr v, const int sign) const                 = 0;
    virtual void          disp() const                                                    = 0;

    /**
     * @brief return the memory size of a block (including the padding for odd numbers if needed)
     * 
     * @return size_t 
     */
    inline size_t get_blockMemSize() const {
        // get the max block size
        size_t total = 1;
        for (int id = 0; id < 3; id++) {
            total *= (size_t)(_nByBlock[id] + _exSize[id] % 2);
        }
        // the nf at the moment of the switchTopo is ALWAYS the one from the output topo!!
        total *= (size_t)_topo_out->nf();
        // return the total size
        return total;
    };
    /**
     * @brief return the buffer size for one proc = number of blocks * blocks memory size
     * 
     * @return size_t 
     */
    inline size_t get_bufMemSize() const {
        // nultiply by the number of blocks
        size_t total = (size_t) std::max(_inBlock[0] * _inBlock[1] * _inBlock[2], _onBlock[0] * _onBlock[1] * _onBlock[2]);
        total *= get_blockMemSize();
        // return the total size
        return total;
    };

   protected:
    void _cmpt_nByBlock();
    void _cmpt_blockDestRankAndTag(const int nBlock[3], const int blockIDStart[3], const FLUPS::Topology *topo, const int *nBlockEachProc, int *destRank, int *destTag);
    void _cmpt_blockSize(const int nBlock[3], const int blockIDStart[3], const int nByBlock[3], const int istart[3], const int iend[3], int *nBlockSize[3]);
    void _cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const FLUPS::Topology *topo, int nBlock[3], int blockIDStart[3], int *nBlockEachProc);

    void _cmpt_commSplit();
    void _setup_subComm(MPI_Comm newcomm, const int nBlock[3], int* destRank, int** count, int** start) ;
};

static inline int gcd(int a, int b) {
    return (a == 0) ? b : gcd(b % a, a);
}

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

#endif