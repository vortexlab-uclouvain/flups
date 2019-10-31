/**
 * @file SwitchTopo.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-30
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


#include "SwitchTopo.hpp"
#include "Topology.hpp"
#include <limits>

void SwitchTopo::_cmpt_nByBlock(int istart[3], int iend[3], int ostart[3], int oend[3],int nByBlock[3]){
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(_inComm,&comm_size);

    int* onProc = (int*)flups_malloc(comm_size * sizeof(int));

    for (int id = 0; id < 3; id++) {
        // get the gcd between send and receive
        int isend = (iend[id] - istart[id]);
        int osend = (oend[id] - ostart[id]);
        // // compute the exchanged size same if from the input or output
        // MPI_Allreduce(&isend, &_exSize[id], 1, MPI_INT, MPI_SUM, _inComm);
        // // we have summed the size nproc(id+1)*size nproc(id+2) * size, so we divide
        // _exSize[id] /= _topo_in->nproc((id+1)%3) * _topo_in->nproc((id+2)%3);

        // // if I am the last one, I decrease the blocksize by one if needed
        // if (_topo_in->rankd(id) == (_topo_in->nproc(id) - 1)) {
        //     isend = isend - _exSize[id] % 2;
        // }
        // if (_topo_out->rankd(id) == (_topo_out->nproc(id) - 1)) {
        //     osend = osend - _exSize[id] % 2;
        // }
        int npoints = gcd(isend,osend);
        // gather on each proc the gcd
        MPI_Allgather(&npoints, 1, MPI_INT, onProc, 1, MPI_INT, _inComm);
        // get the Greatest Common Divider among every process
        int my_gcd = onProc[0];
        for (int ip = 1; ip < comm_size; ip++) {
            my_gcd = gcd(my_gcd, onProc[ip]);
        }
        // store it as the block dimension
        nByBlock[id] = my_gcd;
    }
    flups_free(onProc);
    END_FUNC;
}

/**
 * @brief compute the destination rank for every block on the current processor
 * 
 * @param nBlock the number of block on the current proc (012-indexing)
 * @param blockIDStart the global starting id of the block (0,0,0) in the current topo
 * @param topo the destination topology
 * @param nBlockOnProc the number of block on each proc in the destination topology
 * @param destRank the computed destination rank for each block
 */
void SwitchTopo::_cmpt_blockDestRankAndTag(const int nBlock[3], const int blockIDStart[3], const Topology *topo, const int *startBlockEachProc, const int *nBlockEachProc, int *destRank, int *destTag) {
    BEGIN_FUNC;
    int comm_size;
    MPI_Comm_size(_inComm, &comm_size);
    // go through each block
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        // get the split index
        int bidv[3];
        localSplit(ib, nBlock, 0, bidv, 1);
        // initialize the destrank
        int global_bid[3] = {0, 0, 0};
        int destrankd[3] = {0, 0, 0};
        // determine the dest rank for each dimension
        for (int id = 0; id < 3; id++) {
            // we go trough every rank on the given dim
            global_bid[id] = bidv[id] + blockIDStart[id];
            for (int ir = 0; ir < topo->nproc(id); ir++) {
                // update the destination rank
                destrankd[id] = ir;

                // update the number of block already visited
                int minBlockLocal = startBlockEachProc[id * comm_size + rankindex(destrankd, topo)];
                int maxBlockLocal = minBlockLocal + nBlockEachProc[id * comm_size + rankindex(destrankd, topo)];

                // if we have already visited more block than my block id then we have found the destination rank
                if (global_bid[id] >= minBlockLocal && global_bid[id] < maxBlockLocal) {
                    break;
                }
            }
        }

        // get the global destination rank
        const int destrank = rankindex(destrankd, topo);
        // get the global destination tag
        destRank[ib] = destrank;

        FLUPS_CHECK(destrank < comm_size, "the destination rank is > than the commsize: %d = %d %d %d vs %d", destrank, destrankd[0], destrankd[1], destrankd[2], comm_size, LOCATION);
        if (destTag != NULL) {
            // get the number of block in the destination rank
            int dest_nBlock[3] = {nBlockEachProc[0 * comm_size + destrank],
                                  nBlockEachProc[1 * comm_size + destrank],
                                  nBlockEachProc[2 * comm_size + destrank]};
            // store the destination tag = local block index in the destination rank
            // get the number of block in the destination rank
            int dest_iBlock[3] = {global_bid[0]-startBlockEachProc[0 * comm_size + destrank],
                                  global_bid[1]-startBlockEachProc[1 * comm_size + destrank],
                                  global_bid[2]-startBlockEachProc[2 * comm_size + destrank]};
            // create the tag 
            destTag[ib] = localIndex(0, dest_iBlock[0], dest_iBlock[1], dest_iBlock[2], 0, dest_nBlock, 1);
        }
    }

    //if the communicator of topo is not the same as the reference communicator, we need to adapt the destrank
    //for now, it has been computed in the comm of topo. We thus change for the reference _inComm.
    translate_ranks(nBlock[0] * nBlock[1] * nBlock[2], destRank, topo->get_comm(), _inComm);

    END_FUNC;
}

/**
 * @brief compute the destination rank for every block on the current processor
 * 
 * @param nBlock the number of block on the current proc (012-indexing)
 * @param blockIDStart the global starting id of the block (0,0,0) in the current topo
 * @param topo the destination topology
 * @param nBlockOnProc the number of block on each proc in the destination topology
 * @param destRank the computed destination rank for each block
 */
void SwitchTopo::_cmpt_blockDestRank(const int nBlock[3],const int nByBlock[3],const int shift[3],const int istart[3],const Topology* topo_in, const Topology *topo_out, int *destRank) {
    BEGIN_FUNC;
    int comm_size;
    MPI_Comm_size(_inComm, &comm_size);
    // go through each block
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        // get the split index
        int bidv[3];
        localSplit(ib, nBlock, 0, bidv, 1);
        // initialize the destrank
        int global_id[3] = {0, 0, 0};
        int destrankd[3] = {0, 0, 0};

        // determine the dest rank for each dimension
        for (int id = 0; id < 3; id++) {
            // get the global starting index in my current topo = topo_in
            global_id[id] = bidv[id] * nByBlock[id] + topo_in->nbyproc(id) * topo_in->rankd(id) + istart[id];
            // the (0,0,0) in topo in is located in shift in topo_out
            FLUPS_INFO("block %d starts at %d / %d ",ib,(global_id[id] + shift[id]),topo_out->nbyproc(id));
            destrankd[id] = (global_id[id] + shift[id]) / topo_out->nbyproc(id);
            // if the last proc has more data than the other ones, we need to max the destrank
            destrankd[id] = std::min(destrankd[id],topo_out->nproc(id)-1);
        }
        destRank[ib] = rankindex(destrankd, topo_out);
        
        FLUPS_INFO("block %d will go to proc %d",ib,destRank[ib]);
    }
    //if the communicator of topo is not the same as the reference communicator, we need to adapt the destrank
    //for now, it has been computed in the comm of topo. We thus change for the reference _inComm.
    translate_ranks(nBlock[0] * nBlock[1] * nBlock[2], destRank, topo_out->get_comm(), _inComm);

    END_FUNC;
}

/**
 * @brief compute the size of the blocks inside the given topology
 * 
 * @param nBlock 
 * @param blockIDStart 
 * @param nByBlock 
 * @param topo 
 * @param nBlockSize 
 */
void SwitchTopo::_cmpt_blockSize(const int nBlock[3], const int blockIDStart[3], const int nByBlock[3], const int istart[3], const int iend[3], int *nBlockSize[3]) {
    BEGIN_FUNC;
    // go through each block
    for (int ib2 = 0; ib2 < nBlock[2]; ib2++) {
        for (int ib1 = 0; ib1 < nBlock[1]; ib1++) {
            for (int ib0 = 0; ib0 < nBlock[0]; ib0++) {
                // get the global block index
                const int bidv[3] = {ib0, ib1, ib2};
                const int bid     = localIndex(0, ib0, ib1, ib2, 0, nBlock, 1);
                // determine the size in each direction
                for (int id = 0; id < 3; id++) {
                    //if I am the last block, I forgive a small difference between the blocksizes
                    if (bidv[id] == (nBlock[id] - 1)) {
                        nBlockSize[id][bid] = (iend[id] - istart[id]) - bidv[id] * nByBlock[id];
                    } else {
                        nBlockSize[id][bid] = nByBlock[id];
                    }
                }
            }
        }
    }
    END_FUNC;
}

/**
 * @brief given a topology, try to merge the blocks that go to the same destination
 * 
 * @param [in] topo the topology
 * @param [in] nByBlock the number of unknowns by blocks
 * @param [in] istart the start id on the current proc
 * @param [in] nBlockv the number of blocks in 3D
 * @param [in/out] blockSize the size of each block
 * @param [in/out] blockiStart the starting index of each proc in the local mem (maybe {NULL,NULL,NULL} if no estimation is available)
 * @param [out] nBlock the total number of block 
 * @param [in/out] destRank the destination rank for each block
 * @param [in/out] destTag the destianation tag for each block (may be NULL)
 */
void SwitchTopo::_gather_blocks(const Topology* topo, int nByBlock[3], int istart[3],int iend[3], int nBlockv[3], int* blockSize[3], int* blockiStart[3], int* nBlock, int** destRank) {
    BEGIN_FUNC;
    // get the communicator
    MPI_Comm comm = topo->get_comm();
    int      commsize;
    MPI_Comm_size(comm, &commsize);

    int* nblockToEachProc = (int*)flups_malloc(sizeof(int) * commsize);
    std::memset(nblockToEachProc, 0, sizeof(int) * commsize);

    const int old_nBlock = nBlockv[0] * nBlockv[1] * nBlockv[2];
    
    for (int ib = 0; ib < old_nBlock; ib++) {
        nblockToEachProc[(*destRank)[ib]] += 1;
    }

    int newNBlock = 0;
    for (int ir = 0; ir < commsize; ir++) {
        if (nblockToEachProc[ir] > 0) {
            newNBlock += 1;
        }
    }

    // allocate the new arrays: rank, tag, blocksize, block istart
    int* newBlockSize[3]   = {NULL, NULL, NULL};
    int* newblockiStart[3] = {NULL, NULL, NULL};
    int* newDestRank       = (int*)flups_malloc(sizeof(int) * newNBlock);

    newBlockSize[0]   = (int*)flups_malloc(sizeof(int) * newNBlock);
    newBlockSize[1]   = (int*)flups_malloc(sizeof(int) * newNBlock);
    newBlockSize[2]   = (int*)flups_malloc(sizeof(int) * newNBlock);
    newblockiStart[0] = (int*)flups_malloc(sizeof(int) * newNBlock);
    newblockiStart[1] = (int*)flups_malloc(sizeof(int) * newNBlock);
    newblockiStart[2] = (int*)flups_malloc(sizeof(int) * newNBlock);

    // compute the default value of each array
    int count = 0;
    for (int ir = 0; ir < commsize; ir++) {
        if (nblockToEachProc[ir] > 0) {
            // store the destination rank
            newDestRank[count] = ir;
            // set a blocksize of 0
            newBlockSize[0][count] = -INT_MAX;
            newBlockSize[1][count] = -INT_MAX;
            newBlockSize[2][count] = -INT_MAX;
            // set a istart to max int
            newblockiStart[0][count] = INT_MAX/2;
            newblockiStart[1][count] = INT_MAX/2;
            newblockiStart[2][count] = INT_MAX/2;
            // increment the counter
            count += 1;
        }
    }
    // free the temp array
    flups_free(nblockToEachProc);

    // recompute the number of blocks on each proc and the starting index

    // loop over the blocks and store the information
    for (int nib = 0; nib < newNBlock; nib++) {
        // FLUPS_INFO(">>> looking for new block %d", nib);
        for (int ib = 0; ib < old_nBlock; ib++) {
            // if we have a matching block
            if ((*destRank)[ib] == newDestRank[nib]) {
                // FLUPS_INFO("old block %d go to the same proc: %d vs %d", ib, (*destRank)[ib], newDestRank[nib]);
                // get the block id in 3D
                int bidv[3];
                localSplit(ib, nBlockv, 0, bidv, 1);

                // compute the block size
                int myblockSize[3] = {0,0,0};
                for (int id = 0; id < 3; id++) {
                    myblockSize[id] = (bidv[id] < (nBlockv[id] - 1))? nByBlock[id] : (iend[id] - istart[id]) - bidv[id] * nByBlock[id];
                }
                // int myblockSize[3] = {blockSize[0][ib],blockSize[1][ib],blockSize[2][ib]};
                
                // get the current indexes of the block
                int ib_start[3] = {istart[0] + bidv[0] * nByBlock[0], istart[1] + bidv[1] * nByBlock[1], istart[2] + bidv[2] * nByBlock[2]};
                int ib_end[3]   = {ib_start[0] + myblockSize[0], ib_start[1] + myblockSize[1], ib_start[2] + myblockSize[2]};

                // get the last index of the block
                int nib_start[3] = {newblockiStart[0][nib], newblockiStart[1][nib], newblockiStart[2][nib]};
                int nib_end[3]   = {nib_start[0] + newBlockSize[0][nib], nib_start[1] + newBlockSize[1][nib], nib_start[2] + newBlockSize[2][nib]};
                
                // FLUPS_INFO(">>> old block lim = %d %d %d -> %d %d %d",ib_start[0],ib_start[1],ib_start[2],ib_end[0],ib_end[1],ib_end[2]);
                // FLUPS_INFO(">>> new block lim = %d %d %d -> %d %d %d",nib_start[0],nib_start[1],nib_start[2],nib_end[0],nib_end[1],nib_end[2]);

                // get the new starting index (and overwrittes the INT_MAX if any!!)
                newblockiStart[0][nib] = std::min(nib_start[0], ib_start[0]);
                newblockiStart[1][nib] = std::min(nib_start[1], ib_start[1]);
                newblockiStart[2][nib] = std::min(nib_start[2], ib_start[2]);
                // FLUPS_INFO(">>> new block start = %d %d %d",newblockiStart[0][nib],newblockiStart[1][nib],newblockiStart[2][nib]);
                
                // 2. the size
                newBlockSize[0][nib] = std::max(nib_end[0], ib_end[0]) - newblockiStart[0][nib];
                newBlockSize[1][nib] = std::max(nib_end[1], ib_end[1]) - newblockiStart[1][nib];
                newBlockSize[2][nib] = std::max(nib_end[2], ib_end[2]) - newblockiStart[2][nib];
                // FLUPS_INFO(">>> new block end = %d %d %d",std::max(nib_end[0], ib_end[0]),std::max(nib_end[1], ib_end[1]),std::max(nib_end[2], ib_end[2]));
                // FLUPS_INFO(">>> new block size = %d %d %d",newBlockSize[0][nib],newBlockSize[1][nib],newBlockSize[2][nib]);
            }
        }
        
    }
    // store the new block number
    (*nBlock) = newNBlock;

    // delete the old arrays
    flups_free((*destRank));
    if (blockSize[0] != NULL) {
        flups_free(blockSize[0]);
    }
    if (blockSize[1] != NULL) {
        flups_free(blockSize[1]);
    }
    if (blockSize[2] != NULL) {
        flups_free(blockSize[2]);
    }
    if (blockiStart[0] != NULL) {
        flups_free(blockiStart[0]);
    }
    if (blockiStart[1] != NULL) {
        flups_free(blockiStart[1]);
    }
    if (blockiStart[2] != NULL) {
        flups_free(blockiStart[2]);
    }
    // store the new arrays
    (*destRank) = newDestRank;
    blockSize[0]   = newBlockSize[0];
    blockSize[1]   = newBlockSize[1];
    blockSize[2]   = newBlockSize[2];
    blockiStart[0] = newblockiStart[0];
    blockiStart[1] = newblockiStart[1];
    blockiStart[2] = newblockiStart[2];

    END_FUNC;
}

/**
 * @brief get the tag for a block after a gather_blocks function
 * 
 * @param comm the communicator to use
 * @param inBlock the number of block in the input topo
 * @param onBlock the number of blocks in the output topo
 * @param i2o_destRank the destination rank array fo input to output
 * @param o2i_destRank the destination rank array for output to input
 * @param i2o_destTag the destination tag for input to output
 * @param o2i_destTag the destination tag for output to input
 */
void SwitchTopo::_gather_tags(MPI_Comm comm, const int inBlock, const int onBlock, const int* i2o_destRank, const int* o2i_destRank, int** i2o_destTag, int** o2i_destTag) {
    BEGIN_FUNC;
    // free the memory if it has been allocated
    if((*i2o_destTag) != NULL){
        flups_free(*i2o_destTag);
    }
    if((*o2i_destTag) != NULL){
        flups_free(*o2i_destTag);
    }
    // reallocate it at the correct size
    (*i2o_destTag) = (int*)flups_malloc(sizeof(int) * inBlock);
    (*o2i_destTag) = (int*)flups_malloc(sizeof(int) * onBlock);


    // allocate the requests
    MPI_Request* irequest = (MPI_Request*)flups_malloc(inBlock * sizeof(MPI_Request));
    MPI_Request* orequest = (MPI_Request*)flups_malloc(onBlock * sizeof(MPI_Request));

    //----------------------------
    // for each block in the output configuration, I give its ID to the sender
    for (int ib = 0; ib < onBlock; ib++) {
        MPI_Isend(&ib, 1, MPI_INT, o2i_destRank[ib], 0, comm,orequest+ib);
    }
    // for each block in the input configuration, I receive from the destinator of this block the tag to put in the comm
    for (int ib = 0; ib < inBlock; ib++) {
        MPI_Irecv((*i2o_destTag) + ib, 1, MPI_INT, i2o_destRank[ib], 0, comm,irequest+ib);
    }
    // for for everything to be done
    MPI_Waitall(inBlock,irequest,MPI_STATUSES_IGNORE);
    MPI_Waitall(onBlock,orequest,MPI_STATUSES_IGNORE);
    //----------------------------
    // same but in the backward direction
    for (int ib = 0; ib < inBlock; ib++) {
        MPI_Isend(&ib, 1, MPI_INT, i2o_destRank[ib], 1, comm,irequest+ib);
    }
    for (int ib = 0; ib < onBlock; ib++) {
        MPI_Irecv(&((*o2i_destTag)[ib]), 1, MPI_INT, o2i_destRank[ib], 1, comm, orequest+ib);
    }
    // for for everything to be done
    MPI_Waitall(inBlock,irequest,MPI_STATUSES_IGNORE);
    MPI_Waitall(onBlock,orequest,MPI_STATUSES_IGNORE);

    // free the requests
    flups_free(irequest);
    flups_free(orequest);

    END_FUNC;
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
 * @param blockIDStart the starting id of the block (0,0,0)
 * @param nBlockEachProc the number of blocks on each proc
 */
// void SwitchTopo::_cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const Topology *topo,
//                                      int nBlock[3], int blockIDStart[3], int *startBlockEachProc, int *nBlockEachProc) {
void SwitchTopo::_cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const Topology *topo,int nBlock[3]) {
    BEGIN_FUNC;
    int comm_size;
    MPI_Comm_size(_inComm, &comm_size);
    for (int id = 0; id < 3; id++) {
        // send/recv number of block on my proc
        nBlock[id] = (iend[id] - istart[id]) / nByBlock[id];
        // // get the list of number of procs
        // MPI_Allgather(&(nBlock[id]), 1, MPI_INT, &(nBlockEachProc[comm_size * id]), 1, MPI_INT, topo->get_comm());
        // // set the starting indexes to 0
        // blockIDStart[id] = 0;
        // // compute the starting index
        // const int myrankd  = topo->rankd(id);
        // int       rankd[3] = {topo->rankd(0), topo->rankd(1), topo->rankd(2)};
        // for (int ir = 0; ir < myrankd; ir++) {
        //     // update the rankd
        //     rankd[id] = ir;
        //     // increment the block counter
        //     blockIDStart[id] += nBlockEachProc[comm_size * id + rankindex(rankd, topo)];
        // }
        // do some checks
        FLUPS_CHECK(nBlock[id] > 0, "The number of proc in one direction cannot be 0: istart = %d %d %d to iend = %d %d %d ", istart[0], istart[1], istart[2], iend[0], iend[1], iend[2], LOCATION);

        //everybody needs to know the startID of the first block in each proc
        // MPI_Allgather(&(blockIDStart[id]), 1, MPI_INT, &(startBlockEachProc[comm_size * id]), 1, MPI_INT, topo->get_comm());
    }
    END_FUNC;
}

/**
 * @brief split the _inComm communicator into subcomms
 * 
 * We here find the colors of the call graph, i.e. ranks communicating together have the same color.
 * Once the color are known, we divide the graph into subcomms.
 * 
 * 
 */
void SwitchTopo::_cmpt_commSplit(){
    BEGIN_FUNC;
    // get my rank and use-it as the initial color
    int comm_size, rank;
    MPI_Comm_rank(_inComm,&rank);
    MPI_Comm_size(_inComm,&comm_size);

    
    //-------------------------------------------------------------------------
    /** - Set the starting color and determine who I wish to get in my group */
    //-------------------------------------------------------------------------
    // allocate colors and inMyGroup array
    int*  colors    = (int*)flups_malloc(comm_size * sizeof(int));
    bool* inMyGroup = (bool*)flups_malloc(comm_size * sizeof(bool));


    for (int ir = 0; ir < comm_size; ir++) {
        inMyGroup[ir] = false;
    }
    inMyGroup[rank] = true;

    // do a first pass and give a color + who is in my group
    int mycolor = rank;
    for (int ib = 0; ib < _inBlock; ib++) {
        mycolor                     = std::min(mycolor, _i2o_destRank[ib]);
        inMyGroup[_i2o_destRank[ib]] = true;
    }
    for (int ib = 0; ib < _onBlock; ib++) {
        mycolor                      = std::min(mycolor, _o2i_destRank[ib]);
        inMyGroup[_o2i_destRank[ib]] = true;
    }

    //-------------------------------------------------------------------------
    /** - count how much ranks are in my group and assumes they don't have the same color as I do */
    //-------------------------------------------------------------------------
    int n_left       = 0;  // the global counter of wrong colors
    int n_wrongColor = 0;  // the local counter of wrong colors

    for (int ir = 0; ir < comm_size; ir++) {
        if (inMyGroup[ir]) {
            n_wrongColor += 1;
        }
    }
    // compute among everybody, if we need to continue
    MPI_Allreduce(&n_wrongColor, &n_left, 1, MPI_INT, MPI_SUM, _inComm);

    //-------------------------------------------------------------------------
    /** - Among everybody in group, get the minimum color.
     * The loop is used to force information to travel: 
     * If 0 is in the group of 1 and 1 in the group of 2 and 2 in the group of 3, etc.
     * we need to spread the color 0 among the group
     */
    //-------------------------------------------------------------------------
    int iter = 0;
    while (n_left > 0 && iter < comm_size) {
        // gather the color info from everyone
        MPI_Allgather(&mycolor, 1, MPI_INT, colors, 1, MPI_INT, _inComm);
        // iterate on the proc
        n_wrongColor = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            // if the rank is in my group and it has not the same color
            if (inMyGroup[ir] && (colors[ir] != mycolor)) {
                // we first increment the counter flagging that one is missing
                n_wrongColor += 1;
                // then we change the color if we are able to do so....
                // remove 1 if we are able to solve the color issue <=> my color > colors[ir]
                n_wrongColor = n_wrongColor - (colors[ir] < mycolor);
                // changing the color if possible
                mycolor = std::min(mycolor, colors[ir]);
            }
        }
        // compute among everybody, if we need to continue
        MPI_Allreduce(&n_wrongColor, &n_left, 1, MPI_INT, MPI_SUM, _inComm);
        iter++;
    }
    // if we failed to create the subcom, uses the default one
    if (n_left > 0) {
        _subcomm = _inComm;
        FLUPS_WARNING("I failed to create the subcomm: max iter reached, every group is not complete", LOCATION);
    } else {
        //If there is only 1 color left on all procs, it is 0, and I can still use COMM_WORLD
        int sumColor = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            sumColor += colors[ir];
        }
        // if nleft = 0 -> everybody is inside the same color = the rank = 0
        // we do not create a new comm if it is not necessary
        if (sumColor == 0) {
            // avoids the creation of a communicator
            _subcomm = _inComm;
            FLUPS_INFO("I did not create a new comm since I did not find a way to subdivise the initial comm");
        } else {
            // create the communicator and give a name
            MPI_Comm_split(_inComm, mycolor, rank, &_subcomm);
            std::string commname = "comm-" + std::to_string(mycolor);
            MPI_Comm_set_name(_subcomm, commname.c_str());
        }
    }
    // free the vectors
    flups_free(colors);
    flups_free(inMyGroup);

    END_FUNC;
}

/**
 * @brief setup the lists according to the master and sub communicators
 * 
 * We setup the following lists:
 * - destRank: transformed from the values in the world comm to the values in the new comm.
 * - count: the number of elements send to each proc form this proc
 * - start: the starting position of the data to send to each proc in the buffer
 * 
 * @param nBlock the number of blocks
 * @param BlockSize the size of the blocks
 * @param destRank the destination rank on the world comm. returns the new destination rank in the newcomm
 * @param count the number of information to send to each proc
 * @param start the id in the buffer where the information starts for each proc
 */
void SwitchTopo::_setup_subComm(const int nBlock,int* blockSize[3], int* destRank, int** count, int** start) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the new source & destination ranks    */
    //-------------------------------------------------------------------------
    int inrank, subrank, worldsize;
    MPI_Comm_size(_inComm, &worldsize);
    MPI_Comm_rank(_subcomm, &subrank);
    MPI_Comm_rank(_inComm, &inrank);
    
    // get the ranks of everybody in all communicators
    int* subRanks = (int*)flups_malloc(worldsize * sizeof(int));
    MPI_Allgather(&subrank, 1, MPI_INT, subRanks, 1, MPI_INT, _inComm);


    // int* destRank_cpy = (int*) flups_malloc(nBlock[0] * nBlock[1] * nBlock[2] * sizeof(int));
    // memcpy(destRank,destRank_cpy,nBlock[0] * nBlock[1] * nBlock[2] * sizeof(int));    

    // replace the old ranks by the newest ones
    for (int ib = 0; ib < nBlock; ib++) {
        destRank[ib] = subRanks[destRank[ib]];
    }
    // flups_free(destRank_cpy);
    flups_free(subRanks);
    
    //-------------------------------------------------------------------------
    /** - build the size vector of block to each procs    */
    //-------------------------------------------------------------------------
    if (count != NULL) {
        _cmpt_start_and_count(_subcomm, nBlock,blockSize, destRank, count, start);
    }

    END_FUNC;
}

void SwitchTopo::_cmpt_start_and_count(MPI_Comm comm, const int nBlock,int* blockSize[3], int* destRank, int** count, int** start) {
    BEGIN_FUNC;
    const int nf = std::max(_topo_in->nf(),_topo_out->nf());
    int size;
    MPI_Comm_size(comm, &size);
    // count the number of blocks to each process
    (*count) = (int*)flups_malloc(size * sizeof(int));
    (*start) = (int*)flups_malloc(size * sizeof(int));
    std::memset((*count), 0, size * sizeof(int));
    std::memset((*start), 0, size * sizeof(int));
    
    // count the number of blocks per rank
    for (int ib = 0; ib < nBlock; ib++) {
        // get the size per block
        const int blockMem = get_blockMemSize(ib,nf,blockSize);
        (*count)[destRank[ib]] += blockMem;
    }
    // compute the start indexes
    if (start != NULL) {
        (*start)[0] = 0;
        for (int ir = 1; ir < size; ir++) {
            (*start)[ir] = (*start)[ir - 1] + (*count)[ir - 1];
        }
    }

    END_FUNC;
}

/**
 * @brief setup the suffle plan to do the reordering of the indexes inside a block array form the axis of topo_in to topo_out
 * 
 * The 3D array is split into a rectangular 2D array:
 * - the dimension of the current FRI
 * - the dimension of the targeted FRI
 * 
 * The last dimension is aggregated with eather the current FRI or the targeted, whichever comes on its left
 * e.g.:
 * - if the shuffle is between 2 and 1, the default order is then 2 0 1, hence the dimensions will be (2 * 0) x (1)
 * - if the suffle is between 1 and 2, the  default order is then 1 2 0, hence the dimensions will be (1) x (2 * 0)
 * - if the suffle is between 0 and 1, the  default order is then 1 2 0, hence the dimensions will be (0) x (1 * 2)
 * 
 * @param bSize the block size
 * @param topo_in the topo_in with the current axis
 * @param topo_out the topo_out with the desired axis
 * @param data the data on which to apply the transformation
 * @param shuffle the suffle plan
 */
void SwitchTopo::_setup_shuffle(const int bSize[3], const Topology* topo_in, const Topology* topo_out, double* data, fftw_plan* shuffle) {
    BEGIN_FUNC;

    // the nf will always be the max of both topologies !!
    const int nf = std::max(topo_in->nf(),topo_out->nf());

    // enable the multithreading for this plan
    fftw_plan_with_nthreads(omp_get_max_threads());

    fftw_iodim dims[2];
    // dim[0] = dimension of the targeted FRI (FFTW-convention)
    dims[0].n  = 1;
    dims[0].is = 1;
    dims[0].os = 1;
    // dim[1] = dimension of the current FRI (FFTW-convention)
    dims[1].n  = 1;
    dims[1].is = 1;
    dims[1].os = 1;

    int iaxis[3] = {topo_in->axis(), (topo_in->axis() + 1) % 3, (topo_in->axis() + 2) % 3};
    int oaxis[3] = {topo_out->axis(), (topo_out->axis() + 1) % 3, (topo_out->axis() + 2) % 3};

    // compute the size and the stride of the array
    for (int id = 0; id < 3; id++) {
        if (iaxis[id] != topo_out->axis()) {
            dims[1].n  = dims[1].n * bSize[iaxis[id]];
            dims[0].is = dims[0].is * bSize[iaxis[id]];
        } else {
            break;
        }
    }
    for (int id = 0; id < 3; id++) {
        if (oaxis[id] != topo_in->axis()) {
            dims[0].n  = dims[0].n * bSize[oaxis[id]];
            dims[1].os = dims[1].os * bSize[oaxis[id]];
        } else {
            break;
        }
    }
    // display some info
    FLUPS_INFO("shuffle: setting up the shuffle form %d to %d",topo_in->axis(),topo_out->axis());
    FLUPS_INFO("shuffle: nf = %d, blocksize = %d %d %d",nf,bSize[0],bSize[1],bSize[2]);
    FLUPS_INFO("shuffle: DIM 0: n = %d, is=%d, os=%d",dims[0].n,dims[0].is,dims[0].os);
    FLUPS_INFO("shuffle: DIM 1: n = %d, is=%d, os=%d",dims[1].n,dims[1].is,dims[1].os);

    // plan the real or complex plan
    // the nf is driven by the OUT topology ALWAYS
    if (nf == 1) {
        *shuffle = fftw_plan_guru_r2r(0, NULL, 2, dims, data, data, NULL, FFTW_FLAG);
        FLUPS_CHECK(*shuffle != NULL, "Plan has not been setup", LOCATION);
    } else if (nf == 2) {
        *shuffle = fftw_plan_guru_dft(0, NULL, 2, dims, (fftw_complex*)data, (fftw_complex*)data, FLUPS_FORWARD, FFTW_FLAG);
        FLUPS_CHECK(*shuffle != NULL, "Plan has not been setup", LOCATION);
    }

    END_FUNC;
}

/**
 * @brief Determine and add the weights of the edges in the communication graph
 * 
 * @param sourcesW the weights associated to the edge between other processors communicating to me
 * @param destsW the weights associated to the edge betwenn me communicating to other processors
 */
void SwitchTopo::add_toGraph(int* sourcesW, int* destsW) const{
    BEGIN_FUNC;

    FLUPS_INFO("we call the graph on %d and %d blocks",_inBlock,_onBlock);

    // count the number of out edges
    for (int ib = 0; ib < _inBlock; ib++) {
        destsW[_i2o_destRank[ib]] += _iBlockSize[0][ib]*_iBlockSize[1][ib]*_iBlockSize[2][ib];
    }

    // count the number of in edges
    for (int ib = 0; ib < _onBlock; ib++) {
        sourcesW[_o2i_destRank[ib]] += _oBlockSize[0][ib]*_oBlockSize[1][ib]*_oBlockSize[2][ib];
    }

    // Note: by counting the edges like this on every process, we actually obtain
    // twice the number of edges in the total final graph, as the in and out edges 
    // between 2 processes have been accounted by both procs. However, the weight
    // is relative so it doesnt matter.
    END_FUNC;
}