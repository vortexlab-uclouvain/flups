/**
 * @file SwitchTopo_nb.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-28
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
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

#ifndef SWITCHTOPO_NB_HPP
#define SWITCHTOPO_NB_HPP

#include <cstring>
#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"
#include "Topology.hpp"
#include "Profiler.hpp"
#include "omp.h"
#include "SwitchTopo.hpp"

/**
 * @brief Takes care of the switch between to different topologies
 * 
 * Reorganize the memory between 2 different topologies, also accounting for a
 * "principal axis" which is aligned with the fast rotating index.
 * 
 * Communications are handled by packing data in blocks (i.e. chuncks). These 
 * smaller pieces are sent and recieved asynchronesouly, in a hope for reducing
 * the time required for communication.
 * 
 * The a memory shift can be specified in the switch between the topologies, when
 * there is a need for skipping some data at the left or right side of a given direction.
 * 
 */
class FLUPS::SwitchTopo_nb : public SwitchTopo {
   protected:
    int _selfBlockN=0;
    int* _iselfBlockID = NULL;
    int* _oselfBlockID = NULL;

    int _nByBlock[3]; /**<@brief The number of data per blocks in each dim (!same on each process! and 012-indexing)  */
    int _istart[3]; /**<@brief the starting index for #_topo_in to be inside #_topo_out  */
    int _ostart[3]; /**<@brief the starting index for #_topo_out to be inside #_topo_in  */
    int _iend[3];   /**<@brief the ending index for #_topo_in to be inside #_topo_out  */
    int _oend[3];   /**<@brief the ending index for #_topo_out to be inside #_topo_in  */

    int _inBlock[3];  /**<@brief the local number of block in each dim in the input topology */
    int _onBlock[3];  /**<@brief the local number of block in each dim in the output topology  */

    int _inBlockByProc[3]; /**<@brief The number of blocks in each dim in the input topo = send topo (!different on each process! and 012-indexing) */
    int _onBlockByProc[3]; /**<@brief The number of blocks in each dim in the output topo = recv topo (!different on each process! and 012-indexing) */

    int _iblockIDStart[3]; /**<@brief starting index of the block (0,0,0) in the input topo (012-indexing)    */
    int _oblockIDStart[3]; /**<@brief starting index of the block (0,0,0) in the output topo (012-indexing)    */

    int _ib2o_shift[3]; /**<@brief position in the output topology of the first block (0,0,0) matching the origin of the input topology  */
    int _ob2i_shift[3]; /**<@brief position in the input topology of the first block (0,0,0) matching the origin of the output topology  */

    int* _iBlockSize[3]; /**<@brief The number of data per blocks in each dim for each block (!same on each process! and 012-indexing)  */
    int* _oBlockSize[3]; /**<@brief The number of data per blocks in each dim for each block (!same on each process! and 012-indexing)  */

    opt_int_ptr _i2o_destRank = NULL; /**<@brief The destination rank in the output topo of each block */
    opt_int_ptr _o2i_destRank = NULL; /**<@brief The destination rank in the output topo of each block */

    opt_int_ptr _i2o_destTag = NULL; /**<@brief The destination rank in the output topo of each block */
    opt_int_ptr _o2i_destTag = NULL; /**<@brief The destination rank in the output topo of each block */

    const Topology *_topo_in  = NULL; /**<@brief input topology  */
    const Topology *_topo_out = NULL; /**<@brief  output topology */

    MPI_Request *_i2o_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_i2o_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */
    MPI_Request *_o2i_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_o2i_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */

    opt_double_ptr *_sendBuf = NULL; /**<@brief The send buffer for MPI send */
    opt_double_ptr *_recvBuf = NULL; /**<@brief The recv buffer for MPI recv */

    Profiler* _prof = NULL;

   public:
    SwitchTopo_nb(const Topology *topo_input, const Topology *topo_output, const int shift[3],Profiler* prof);
    ~SwitchTopo_nb();

    void setup_buffers(opt_double_ptr _sendBuf,opt_double_ptr _recvBuf);
    void execute(opt_double_ptr v, const int sign) const;

    /**
     * @brief return the memory size of a block (including the padding for odd numbers if needed)
     * 
     * @return size_t 
     */
    inline size_t get_blockMemSize() const {
        // // get the max block size
        // size_t total = 1;
        // for (int id = 0; id < 3; id++) {
        //     total *= (size_t)(_nByBlock[id] + _exSize[id] % 2);
        // }
        // // the nf at the moment of the switchTopo is ALWAYS the one from the output topo!!
        // total *= (size_t)_topo_out->nf();
        // // return the total size
        // return total;
        // the nf at the moment of the switchTopo is ALWAYS the one from the output topo!!
        const int nf = _topo_out->nf();
        // get the max block size
        size_t maxsize = 0;
        for (int ib = 0; ib < _inBlock[0] * _inBlock[1] * _inBlock[2]; ib++) {
            maxsize = std::max(maxsize,((size_t)_iBlockSize[0][ib])*((size_t)_iBlockSize[1][ib])*((size_t)_iBlockSize[2][ib])* ((size_t)nf));
        }
        // get the max block size
        for (int ib = 0; ib < _onBlock[0] * _onBlock[1] * _onBlock[2]; ib++) {
            maxsize = std::max(maxsize,((size_t)_oBlockSize[0][ib])*((size_t)_oBlockSize[1][ib])*((size_t)_oBlockSize[2][ib])*((size_t)nf));
        }
        // return
        return maxsize;
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

    void disp() const;
    void disp_rankgraph(const int id_in,const int id_out) const;
};

void SwitchTopo_test();

#endif