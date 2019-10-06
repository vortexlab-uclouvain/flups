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

using namespace FLUPS;


void SwitchTopo::_cmpt_nByBlock(){
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);

    int* onProc = (int*)fftw_malloc(comm_size * sizeof(int));

    for (int id = 0; id < 3; id++) {
        // get the gcd between send and receive
        int isend = (_iend[id] - _istart[id]);
        int osend = (_oend[id] - _ostart[id]);
        // compute the exchanged size same if from the input or output
        MPI_Allreduce(&isend, &_exSize[id], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        // we have summed the size nproc(id+1)*size nproc(id+2) * size, so we divide
        _exSize[id] /= _topo_in->nproc((id+1)%3) * _topo_in->nproc((id+2)%3);

        // if I am the last one, I decrease the blocksize by one if needed
        if (_topo_in->rankd(id) == (_topo_in->nproc(id) - 1)) {
            isend = isend - _exSize[id] % 2;
        }
        if (_topo_out->rankd(id) == (_topo_out->nproc(id) - 1)) {
            osend = osend - _exSize[id] % 2;
        }
        int npoints = gcd(isend,osend);
        // gather on each proc the gcd
        MPI_Allgather(&npoints, 1, MPI_INT, onProc, 1, MPI_INT, MPI_COMM_WORLD);
        // get the Greatest Common Divider among every process
        int my_gcd = onProc[0];
        for (int ip = 1; ip < comm_size; ip++) {
            my_gcd = gcd(my_gcd, onProc[ip]);
        }
        // store it as the block dimension
        _nByBlock[id] = my_gcd;
    }
    fftw_free(onProc);
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
void SwitchTopo::_cmpt_blockDestRankAndTag(const int nBlock[3], const int blockIDStart[3], const FLUPS::Topology *topo, const int *nBlockEachProc, int *destRank, int *destTag) {
    BEGIN_FUNC;
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    // go through each block
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        // get the split index
        int bidv[3];
        localSplit(ib, nBlock, 0, bidv, 1);
        // initialize the destrank
        int destrankd[3] = {0, 0, 0};
        int local_bid[3] = {0, 0, 0};

        // determine the dest rank for each dimension
        for (int id = 0; id < 3; id++) {
            // we go trough every rank on the given dim
            int block_count = 0;
            for (int ir = 0; ir < topo->nproc(id); ir++) {
                // update the destination rank
                destrankd[id] = ir;
                // compute the local block id in this rank
                local_bid[id] = bidv[id] + blockIDStart[id] - block_count;
                // update the number of block already visited
                block_count += nBlockEachProc[id * comm_size + rankindex(destrankd, topo)];
                // if we have already visited more block than my block id then we have found the destination rank
                if ((bidv[id] + blockIDStart[id]) < block_count) {
                    break;
                }
            }
        }

        // get the global destination rank
        const int destrank = rankindex(destrankd, topo);
        // get the global destination rank
        destRank[ib] = destrank;
        FLUPS_CHECK(destrank < comm_size, "the destination rank is > than the commsize: %d = %d %d %d vs %d", destrank, destrankd[0], destrankd[1], destrankd[2], comm_size, LOCATION);
        if (destTag != NULL) {
            // get the number of block in the destination rank
            int dest_nBlock[3] = {nBlockEachProc[0 * comm_size + destrank],
                                  nBlockEachProc[1 * comm_size + destrank],
                                  nBlockEachProc[2 * comm_size + destrank]};
            // store the destination tag = local block index in the destination rank
            destTag[ib] = localIndex(0, local_bid[0], local_bid[1], local_bid[2], 0, dest_nBlock, 1);
        }
    }
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
void SwitchTopo::_cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const FLUPS::Topology *topo,
                                     int nBlock[3], int blockIDStart[3], int *nBlockEachProc) {
    BEGIN_FUNC;
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
        // do some checks
        FLUPS_CHECK(nBlock[id] > 0, "The number of proc in one direction cannot be 0: istart = %d %d %d to iend = %d %d %d ", istart[0], istart[1], istart[2], iend[0], iend[1], iend[2], LOCATION);
    }
    END_FUNC;
}

void SwitchTopo::_cmpt_commSplit(){
    BEGIN_FUNC;
    // get my rank and use-it as the initial color
    int comm_size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);

    int mycolor = rank;
    
    // allocate colors and inMyGroup array
    int*  colors    = (int*)fftw_malloc(comm_size * sizeof(int));
    bool* inMyGroup = (bool*)fftw_malloc(comm_size * sizeof(bool));

    for (int ir = 0; ir < comm_size; ir++) {
        inMyGroup[ir] = false;
    }

    // do a first pass and give a color + who is in my group
    for (int ib = 0; ib < _inBlock[0] * _inBlock[1] * _inBlock[2]; ib++) {
        mycolor                     = std::min(mycolor, _i2o_destRank[ib]);
        inMyGroup[_i2o_destRank[ib]] = true;
    }
    for (int ib = 0; ib < _onBlock[0] * _onBlock[1] * _onBlock[2]; ib++) {
        mycolor                      = std::min(mycolor, _o2i_destRank[ib]);
        inMyGroup[_o2i_destRank[ib]] = true;
    }

    // count how much is should be in my group
    // by default we assume that nobody is in the same group
    int nleft = 0;
    for (int ir = 0; ir < comm_size; ir++) {
        if (inMyGroup[ir]) {
            nleft += 1;
        }
    }
    // continue while we haven't found a solution
    while (nleft > 0) {
        // gather the color info from everyone
        MPI_Allgather(&mycolor, 1, MPI_INT, colors, 1, MPI_INT, MPI_COMM_WORLD);
        // iterate on the proc
        int n_notInMyGroup = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            // if it is reachable and the color is not already the same
            if (inMyGroup[ir] && (colors[ir] != mycolor)) {
                // we first increment the counter flagging that one is missing
                n_notInMyGroup += 1;
                // then we solve the problem if we are able to do so....
                // remove 1 if we are able to solve the issue <=> my color > colors[ir]
                n_notInMyGroup = n_notInMyGroup - (colors[ir] < mycolor);
                // changing if possible
                mycolor = std::min(mycolor, colors[ir]);
            }
        }
        // compute among everybody, if we need to continue
        MPI_Allreduce(&n_notInMyGroup, &nleft, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        FLUPS_INFO("stil %d to find (@ my proc: %d)", nleft, n_notInMyGroup);
    }
    // create the communicator and give a name
    MPI_Comm_split(MPI_COMM_WORLD, mycolor, rank, &_subcomm);

    std::string commname = "comm-" + std::to_string(mycolor);
    MPI_Comm_set_name(_subcomm, commname.c_str());

    // free the vectors
    fftw_free(colors);
    fftw_free(inMyGroup);
    END_FUNC;
}

/**
 * @brief setup the comm lists
 * 
 * We setup the following lists:
 * - destRank: transformed from the values in the world comm to the values in the new comm.
 * - count: the number of elements send to each proc form this proc
 * - start: the starting position of the data to send to each proc in the buffer
 * 
 * @param newcomm the new comm
 * @param nBlock the number of blocks
 * @param destRank the destination rank on the world comm. returns the new destination rank in the newcomm
 * @param count the number of information to send to each proc
 * @param start the id in the buffer where the information starts for each proc
 */
void SwitchTopo::_setup_subComm(MPI_Comm newcomm, const int nBlock[3], int* destRank, int** count, int** start) {
    BEGIN_FUNC;
    //-------------------------------------------------------------------------
    /** - get the new destination ranks    */
    //-------------------------------------------------------------------------
    int newrank, worldsize;
    MPI_Comm_rank(newcomm, &newrank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);

    // get the new ranks from my old friends in the old communicator
    int* newRanks = (int*)fftw_malloc(worldsize * sizeof(int));
    MPI_Allgather(&newrank, 1, MPI_INT, newRanks, 1, MPI_INT, MPI_COMM_WORLD);
    // replace the old ranks by the newest ones
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        destRank[ib] = newRanks[destRank[ib]];
    }
    fftw_free(newRanks);

    //-------------------------------------------------------------------------
    /** - build the size vector of block to each procs    */
    //-------------------------------------------------------------------------
    if (count != NULL) {
        int subsize;
        MPI_Comm_size(_subcomm, &subsize);
        // count the number of blocks to each process
        (*count) = (int*)fftw_malloc(subsize * sizeof(int));
        (*start) = (int*)fftw_malloc(subsize * sizeof(int));
        std::memset((*count), 0, subsize * sizeof(int));
        std::memset((*start), 0, subsize * sizeof(int));
        // get the size per block
        const int blockMem = get_blockMemSize();
        // count the number of blocks per rank
        for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
            (*count)[destRank[ib]] += blockMem;
        }
        // compute the start indexes
        if (start != NULL) {
            (*start)[0] = 0;
            for (int ir = 1; ir < subsize; ir++) {
                (*start)[ir] = (*start)[ir - 1] + (*count)[ir - 1];
            }
        }
    }
    END_FUNC;
}