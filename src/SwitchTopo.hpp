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
#include "topology.hpp"

typedef int bcoord[3];

/**
 * @brief Do the switch between to different topologies
 * 
 * The switch between two topologies is driven by the shift between the topologies.
 * 
 */
class SwitchTopo {
   protected:
    const Topology *_topo_in;  /**<@brief input topology  */
    const Topology *_topo_out; /**<@brief  output topology */

    int _inBlockByProc[3]; /**<@brief The number of blocks in each dim in the input topo = send topo (!different on each process! and 012-indexing) */
    int _onBlockByProc[3]; /**<@brief The number of blocks in each dim in the output topo = recv topo  */

    int _nByBlock[3]; /**<@brief The number of data per blocks in each dim (!same on each process! and 012-indexing)  */
    int _inBlock[3]; /**<@brief the local number of block in each dim in the input topology */
    int _onBlock[3]; /**<@brief the local number of block in each dim in the output topology  */

    int _iblockIDStart[3];/**<@brief starting index of the block (0,0,0) in the input topo (012-indexing)    */
    int _oblockIDStart[3];/**<@brief starting index of the block (0,0,0) in the output topo (012-indexing)    */

    MPI_Request* _sendRequest;
    MPI_Request* _recvRequest;

    int _ib2o_shift[3]; /**<@brief position in the output topology of the first block (0,0,0) matching the origin of the input topology  */
    int _ob2i_shift[3]; /**<@brief position in the input topology of the first block (0,0,0) matching the origin of the output topology  */

    int* _i2o_destRank; /**<@brief The destination rank in the output topo of each block */
    int* _o2i_destRank; /**<@brief The destination rank in the output topo of each block */

    int* _i2o_destTag; /**<@brief The destination rank in the output topo of each block */
    int* _o2i_destTag; /**<@brief The destination rank in the output topo of each block */

    opt_double_ptr* _sendBuf;
    opt_double_ptr* _recvBuf;


    //====================== OLD STUFFS =============================

    int _ishift[3]; /**<@brief the shift index for #_topo_in to match coordinate (0,0,0) in #__topo_out   */
    int _oshift[3]; /**<@brief the shift index for #_topo_out to match coordinate (0,0,0) in #_topo_in*/

    int _istart[3]; /**<@brief the starting index for #_topo_in to be inside #_topo_out  */
    int _iend[3];   /**<@brief the ending index for #_topo_in to be inside #_topo_out  */
    int _ostart[3]; /**<@brief the starting index for #_topo_out to be inside #_topo_in  */
    int _oend[3];   /**<@brief the ending index for #_topo_out to be inside #_topo_in  */

    int *_nsend          = NULL; /**<@brief number of unknowns to send to each proc  */
    int *_nrecv          = NULL; /**<@brief number of unknowns to receiv from each proc */
    int *_ssend          = NULL; /**<@brief start index in the buffer memory for data to send to each proc */
    int *_srecv          = NULL; /**<@brief start index in the buffer memory for data to receive from each proc */
    int *_count          = NULL; /**<@brief counter to know where to put the current readen data  */
    int *_naxis_i2o_send = NULL; /**<@brief naxis when sending from the _topo_in to _topo_out  */
    int *_naxis_o2i_send = NULL; /**<@brief naxis when sending from the _topo_out to _topo_in */
    int *_naxis_i2o_recv = NULL; /**<@brief naxis when receiving from the _topo_in to _topo_out  */
    int *_naxis_o2i_recv = NULL; /**<@brief naxis when receiving from the _topo_out to _topo_in */

    double *_bufsend = NULL; /**<@brief buffer to send data */
    double *_bufrecv = NULL; /**<@brief buffer to receive data */

   public:
    SwitchTopo(const Topology *topo_input, const Topology *topo_output, const int shift[3]);
    ~SwitchTopo();

    void execute(opt_double_ptr v, const int sign);

    void disp();
};

void SwitchTopo_test();

static inline int gcd(int a, int b) {
    return (a == 0) ? b : gcd(b % a, a);
}

static inline void cmpt_exchangedSize(const Topology* topo_in, const Topology* topo_out, int exchangeSize[3]){
    exchangeSize[0] = std::min(topo_in->nglob(0),topo_out->nglob(0));
    exchangeSize[1] = std::min(topo_in->nglob(1),topo_out->nglob(1));
    exchangeSize[2] = std::min(topo_in->nglob(2),topo_out->nglob(2));
}

// /**
//  * @brief Compute the mean number by block given the current and the other topology
//  * 
//  * @param current the current topology
//  * @param other the other topology 
//  * @param nByBlock the mean number of unkowns by proc
//  */
// static inline void cmpt_nByBlock(const Topology *current, const Topology *other, int nByBlock[3]) {
//     // get the exchanged size
//     int exSize[3];
//     cmpt_exchangedSize(current, other, exSize);
//     // compute the number of unknowns by proc
//     for (int id = 0; id < 3; id++) {
//         // get the least common multiplier between the two topologies
//         int lcm = (current->nproc(id) * other->nproc(id)) / gcd(current->nproc(id), other->nproc(id));
//         // the size is the min of the two topos
//         nByBlock[id] = exSize[id] / lcm;
//     }
// }

// /**
//  * @brief get the block size on the current toplogy
//  * 
//  * @param current the current topology
//  * @param exchangeSize the total size to exchange
//  * @param nBlock the number of blocks (012-indexing)
//  * @param nByBlock the mean number of points per block (012-indexing)
//  * @param blockID the local block id (012-indexing)
//  * @param size the returned size of the block (012-indexing)
//  */
// static inline void cmpt_blockSize(const Topology* current,const Topology* other, const int nBlock[3],const int nByBlock[3], const int blockID[3], int size[3]){
//     // get the exchanged size between topo's
//     int exSize[3];
//     cmpt_exchangedSize(current,other,exSize);
//     // get the actual size
//     for(int id=0; id<3; id++){
//         // is the block the last one in the direction and get the starting index of the block
//         const bool islast   = blockID[id] == (nBlock[id] - 1) && current->rankd(id) == (current->nproc(id) - 1);
//         const int  id_start = current->rankd(id) * current->nbyproc(id) + nByBlock[id] * blockID[id];
//         // get the actual size of the block
//         size[id] = (islast) ? std::max(nByBlock[id], exSize[id] - id_start) : nByBlock[id];
//     }
// }

/**
 * @brief return the local block index from the splitted indexing (012-indexing)
 * 
 * @param ib0 block index in the 0 direction
 * @param ib1 block index in the 1 direction
 * @param ib2 block index in the 2 direction
 * @param nBlock the local number of blocks
 * @return int 
 */
static inline int blockID(const int ib0, const int ib1, const int ib2, const int nBlock[3]) {
    return ib0 + nBlock[0] * (ib1 + nBlock[1] * ib2);
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
 * @param nf 
 * @return size_t 
 */
static inline size_t localIndex(const int axsrc, const int i0, const int i1, const int i2,
                                const int axtrg, const int size[3], const int nf) {
    // const int ax0 = axis;
    // const int ax1 = (ax0 + 1) % 3;
    // const int ax2 = (ax0 + 2) % 3;
    // return i0 * nf + size[ax0] * nf * (i1 + size[ax1] * i2);

    const int i[3] = {i0, i1, i2};
    const int dax0 = (3 + axtrg - axsrc) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;

    // return localindex_xyz(i[0], i[1], i[2], topo);
    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * i[dax2]);
}

static inline void localSplit(const int id, const int size[3], const int axtrg, int idv[3]) {
    const int ax0  = axtrg;
    const int ax1  = (ax0 + 1) % 3;
    const int ax2  = (ax0 + 2) % 3;

    idv[ax0] = id % size[ax0];
    idv[ax1] = (id % (size[0] * size[ax1])) / size[ax0];
    idv[ax2] = id / (size[ax0] * size[ax1]);
}

static inline void blocksplit(const int bid, const int nBlock[3], int blockid[3]) {
    blockid[0] = bid % nBlock[0];
    blockid[1] = (bid % (nBlock[0] * nBlock[1])) / nBlock[0];
    blockid[2] = bid / (nBlock[0] * nBlock[1]);
}

inline static size_t localIndex_withBlock_ao(const int i0, const int i1, const int i2, const Topology *topo, const int blockID[3], const int nByBlock[3]) {
    const int nf  = topo->nf();
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int loc_i0 = nByBlock[ax0] * blockID[ax0] + i0;
    const int loc_i1 = nByBlock[ax1] * blockID[ax1] + i1;
    const int loc_i2 = nByBlock[ax2] * blockID[ax2] + i2;

    return loc_i0 * nf + topo->nloc(ax0) * nf * (loc_i1 + topo->nloc(ax1) * loc_i2);
}

inline static size_t localindex_withBlock(const int axis, const int i0, const int i1, const int i2, const Topology *topo, const int blockID[3], const int nByBlock[3]) {
    const int nf   = topo->nf();
    
    // compute the shift to perform from the axis reference to
    const int dax0 = (3 + topo->axis() - axis) % 3;
    const int dax1 = (dax0 + 1) % 3;
    const int dax2 = (dax0 + 2) % 3;

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    const int loc_i[3] = {nByBlock[ax0] * blockID[ax0]+ i0, nByBlock[ax1] * blockID[ax1]+i1, nByBlock[ax2] * blockID[ax2]+i2};

    // return localindex_xyz(i[0], i[1], i[2], topo);
    return loc_i[dax0] * nf + topo->nloc(ax0) * nf * (loc_i[dax1] + topo->nloc(ax1) * loc_i[dax2]);
}

inline static size_t blockLocalIndex_ao(const int i0, const int i1, const int i2, const Topology *topo, const int blockSize[3]) {
    const int nf  = topo->nf();
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;

    return i0 * nf + blockSize[ax0] * nf * (i1 + blockSize[ax1] * i2);
}


/**
 * @brief For a given block in the current topology compute the corresponding rank in the other topology
 * 
 * @param current the current topology
 * @param other the other topology
 * @param b_shift the position of the current topology (0,0,0) in the other topology (012-indexing)
 * @param nByBlock the mean number of points per block (012-indexing)
 * @param bid the local block id (012-indexing)
 * @return int 
 */
static inline int cmpt_blockRank(const Topology* current, const Topology* other,const int b_shift[3],const int nBlockByProc[3],const int nByBlock[3], const int bid[3])
{
    int rankd[3];
    rankd[0] = ((current->rankd(0) * nBlockByProc[0] + bid[0])* nByBlock[0] + b_shift[0]) / other->nbyproc(0);
    rankd[1] = ((current->rankd(1) * nBlockByProc[1] + bid[1])* nByBlock[1] + b_shift[1]) / other->nbyproc(1);
    rankd[2] = ((current->rankd(2) * nBlockByProc[2] + bid[2])* nByBlock[2] + b_shift[2]) / other->nbyproc(2);

    return rankindex(rankd,other);
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
static inline void cmpt_blockDestRankAndTag(const int nBlock[3], const int blockIDStart[3], const Topology *topo, const int *nBlockEachProc,
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
                destRank[blockID(ib0, ib1, ib2, nBlock)] = destrank;
                // get the number of block in the destination rank
                int dest_nBlock[3] = {nBlockEachProc[0 * comm_size + destrank],
                                      nBlockEachProc[1 * comm_size + destrank],
                                      nBlockEachProc[2 * comm_size + destrank]};
                // store the destination tag = local block index in the destination rank
                destTag[blockID(ib0, ib1, ib2, nBlock)] = blockID(local_bid[0], local_bid[1], local_bid[2], dest_nBlock);
            }
        }
    }
}

static inline void cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const Topology *topo,
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
        const int myrankd        = topo->rankd(id);
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