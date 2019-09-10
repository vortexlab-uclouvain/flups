/**
 * @file SwitchTopo.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-25
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef REODER_MPI_HPP
#define REODER_MPI_HPP

#include <cstring>
#include "defines.hpp"
#include "hdf5_io.hpp"
#include "mpi.h"
#include "Topology.hpp"
#include "Profiler.hpp"
#include "omp.h"

typedef int bcoord[3];

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
class FLUPS::SwitchTopo {
   protected:
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

    opt_int_ptr _i2o_destRank = NULL; /**<@brief The destination rank in the output topo of each block */
    opt_int_ptr _o2i_destRank = NULL; /**<@brief The destination rank in the output topo of each block */

    opt_int_ptr _i2o_destTag = NULL; /**<@brief The destination rank in the output topo of each block */
    opt_int_ptr _o2i_destTag = NULL; /**<@brief The destination rank in the output topo of each block */

    const Topology *_topo_in  = NULL; /**<@brief input topology  */
    const Topology *_topo_out = NULL; /**<@brief  output topology */

    MPI_Request *_sendRequest = NULL; /**<@brief The MPI Request generated on the send */
    MPI_Request *_recvRequest = NULL; /**<@brief The MPI Request generated on the recv */

    opt_double_ptr *_sendBuf = NULL; /**<@brief The send buffer for MPI send */
    opt_double_ptr *_recvBuf = NULL; /**<@brief The recv buffer for MPI recv */

    Profiler* _prof = NULL;

   public:
    SwitchTopo(const Topology *topo_input, const Topology *topo_output, const int shift[3],Profiler* prof);
    
    ~SwitchTopo();

    void execute(opt_double_ptr v, const int sign);

    void disp();
};

void SwitchTopo_test();

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
    const int ax0 = axtrg;
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    const int size0 = (size[ax0]*nf);

    idv[ax0] = id % size0;
    idv[ax1] = (id % (size0 * size[ax1])) / size0;
    idv[ax2] = id / (size0 * size[ax1]);
}
static inline void localSplit(const size_t id, const int size[3], const int axtrg, int *id0,int *id1,int *id2, const int nf) {
    const int ax0 = axtrg;
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    const int size0 = (size[ax0]*nf);

    (*id0) = id % size0;
    (*id1) = (id % (size0 * size[ax1])) / size0;
    (*id2) = id / (size0 * size[ax1]);
}

/**
 * @brief compute the destination rank for every block on the current processor
 * 
 * @param nBlockByProc the number of block on the current proc (012-indexing)
 * @param blockIDStart the global starting id of the block (0,0,0) in the current topo
 * @param topo the destination topology
 * @param nBlockOnProc the number of block on each proc in the destination topology
 * @param destRank the computed destination rank for each block
 */
static inline void cmpt_blockDestRankAndTag(const int nBlock[3], const int blockIDStart[3], const FLUPS::Topology *topo, const int *nBlockEachProc,
                                            int *destRank, int *destTag) {
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    // go through each block
    for (int ib2 = 0; ib2 < nBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < nBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < nBlock[0]; ib0++) {
                // get the global block index
                const int bidv[3] = {ib0 + blockIDStart[0],
                                     ib1 + blockIDStart[1],
                                     ib2 + blockIDStart[2]};
                // initialize the destrank
                int destrankd[3] = {0, 0, 0};
                int local_bid[3] = {0, 0, 0};

                // determine the dest rank for each dimension
                for (int id = 0; id < 3; id++) {
                    // we go trough every proc on the dim
                    int block_count = 0;
                    for (int ir = 0; ir < topo->nproc(id); ir++) {
                        // update the destination rank
                        destrankd[id] = ir;
                        // compute the local block id in this rank
                        local_bid[id] = bidv[id] - block_count;
                        // update the number of block already visited
                        block_count += nBlockEachProc[id * comm_size + rankindex(destrankd, topo)];
                        // if we have already visited more block than my block id then we have found the destination rank
                        if (bidv[id] < block_count) {
                            break;
                        }
                    }
                }

                // get the global destination rank
                int destrank = rankindex(destrankd, topo);
                // get the global destination rank
                int bid = localIndex(0,ib0,ib1,ib2,0,nBlock,1);
                destRank[bid] = destrank;
                // get the number of block in the destination rank
                int dest_nBlock[3] = {nBlockEachProc[0 * comm_size + destrank],
                                      nBlockEachProc[1 * comm_size + destrank],
                                      nBlockEachProc[2 * comm_size + destrank]};
                // store the destination tag = local block index in the destination rank
                destTag[bid] = localIndex(0,local_bid[0], local_bid[1], local_bid[2],0,dest_nBlock,1);
            }
        }
    }
}

/**
 * @brief compute the number of blocks, the starting indexes of the block (0,0,0) and the number of block in each proc
 * 
 * This function computes several usefull indexes for the block:
 * - the number of blocks on the current procs
 * - the starting index in the topo of the block (0,0,0)
 * - the number of block on each proc.
 * 
 * For a given proc, nBlockEachProc[comm_size * id + ip] is the number of proc in the dimension id on the proc ip
 * 
 * @param istart the starting indexes on this proc
 * @param iend the end indexes on this proc
 * @param nByBlock the number of unkowns in one block (012-indexing)
 * @param topo the current topology
 * @param nBlock the number of block in this proc
 * @param blockIDStart the starting point of the block (0,0,0)
 * @param nBlockEachProc the number of procs on each proc
 */
static inline void cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const FLUPS::Topology *topo,
                                     int nBlock[3], int blockIDStart[3], int *nBlockEachProc) {
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    for (int id = 0; id < 3; id++) {
        // send/recv number of block on my proc
        nBlock[id] = (iend[id] - istart[id]) / nByBlock[id];
        // get the list of number of procs
        MPI_Allgather(&(nBlock[id]), 1, MPI_INT, &(nBlockEachProc[comm_size * id]), 1, MPI_INT, MPI_COMM_WORLD);
        // set the starting indexes to 0
        blockIDStart[id] = 0;
        // compute the starting index
        const int myrankd  = topo->rankd(id);
        int       rankd[3] = {topo->rankd(0), topo->rankd(1), topo->rankd(2)};
        for (int ir = 0; ir < myrankd; ir++) {
            // update the rankd
            rankd[id] = ir;
            // increment the block counter
            blockIDStart[id] += nBlockEachProc[comm_size * id + rankindex(rankd, topo)];
        }
    }
}

#endif