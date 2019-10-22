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


void SwitchTopo::_cmpt_nByBlock(){
    BEGIN_FUNC;

    int comm_size;
    MPI_Comm_size(_inComm,&comm_size);

    int* onProc = (int*)flups_malloc(comm_size * sizeof(int));

    for (int id = 0; id < 3; id++) {
        // get the gcd between send and receive
        int isend = (_iend[id] - _istart[id]);
        int osend = (_oend[id] - _ostart[id]);
        // compute the exchanged size same if from the input or output
        MPI_Allreduce(&isend, &_exSize[id], 1, MPI_INT, MPI_SUM, _inComm);
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
        MPI_Allgather(&npoints, 1, MPI_INT, onProc, 1, MPI_INT, _inComm);
        // get the Greatest Common Divider among every process
        int my_gcd = onProc[0];
        for (int ip = 1; ip < comm_size; ip++) {
            my_gcd = gcd(my_gcd, onProc[ip]);
        }
        // store it as the block dimension
        _nByBlock[id] = my_gcd;
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
void SwitchTopo::_cmpt_blockDestRankAndTag(const int nBlock[3], const int blockIDStart[3], const Topology *topo, const int *nBlockEachProc, int *destRank, int *destTag) {
    BEGIN_FUNC;
    int comm_size;
    MPI_Comm_size(_inComm, &comm_size);
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
        // get the global destination tag
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

    //if the communicator of topo is not the same as the reference communicator, we need to adapt the destrank
    //for now, it has been computed in the comm of topo. We thus change for the reference _inComm.
    translate_ranks(nBlock[0] * nBlock[1] * nBlock[2], destRank, topo->get_comm(), _inComm);

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
void SwitchTopo::_cmpt_blockIndexes(const int istart[3], const int iend[3], const int nByBlock[3], const Topology *topo,
                                     int nBlock[3], int blockIDStart[3], int *nBlockEachProc) {
    BEGIN_FUNC;
    int comm_size;
    MPI_Comm_size(_inComm, &comm_size);
    for (int id = 0; id < 3; id++) {
        // send/recv number of block on my proc
        nBlock[id] = (iend[id] - istart[id]) / nByBlock[id];
        // get the list of number of procs
        MPI_Allgather(&(nBlock[id]), 1, MPI_INT, &(nBlockEachProc[comm_size * id]), 1, MPI_INT, _inComm);
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
    MPI_Comm_rank(_inComm,&rank);
    MPI_Comm_size(_inComm,&comm_size);

    int mycolor = rank;
    
    // allocate colors and inMyGroup array
    int*  colors    = (int*)flups_malloc(comm_size * sizeof(int));
    bool* inMyGroup = (bool*)flups_malloc(comm_size * sizeof(bool));

    for (int ir = 0; ir < comm_size; ir++) {
        inMyGroup[ir] = false;
    }
    inMyGroup[rank] = true;

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
    int iter = 0;
    while (nleft > 0 && iter < comm_size) {
        // gather the color info from everyone
        MPI_Allgather(&mycolor, 1, MPI_INT, colors, 1, MPI_INT, _inComm);
        // iterate on the proc
        int n_notInMyGroup = 0;
        for (int ir = 0; ir < comm_size; ir++) {
            // if(rank==3){
                printf("[ITER%d] color[%d] %d,  inMyGroup:%d\n",iter,ir,colors[ir],inMyGroup[ir]);
            // }
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
        MPI_Allreduce(&n_notInMyGroup, &nleft, 1, MPI_INT, MPI_SUM, _inComm);
        iter++;
    }
    // FLUPS_CHECK(iter<comm_size,"Could not divide the current comm in subcom.",LOCATION);
    if(iter>=comm_size){
        _subcomm = _inComm;
    }

    //If there is only 1 color left on all procs, it is 0, and I can still use COMM_WORLD
    nleft=0;
    for(int ir = 0; ir < comm_size; ir++){
        nleft+=colors[ir];
    }
    if(nleft==0){
        // avoids the creation of a communicator
        _subcomm = _inComm;
        FLUPS_INFO("I did not create a new comm since I did not find a way to subdivise master",LOCATION);
        
    } else {
        // create the communicator and give a name
        MPI_Comm_split(_inComm, mycolor, rank, &_subcomm);

        std::string commname = "comm-" + std::to_string(mycolor);
        MPI_Comm_set_name(_subcomm, commname.c_str());
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
 * @param destRank the destination rank on the world comm. returns the new destination rank in the newcomm
 * @param count the number of information to send to each proc
 * @param start the id in the buffer where the information starts for each proc
 */
void SwitchTopo::_setup_subComm(const int nBlock[3], int* destRank, int** count, int** start) {
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
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        destRank[ib] = subRanks[destRank[ib]];
    }
    // flups_free(destRank_cpy);
    flups_free(subRanks);
    
    //-------------------------------------------------------------------------
    /** - build the size vector of block to each procs    */
    //-------------------------------------------------------------------------
    if (count != NULL) {
        _cmpt_start_and_count(_subcomm, nBlock, destRank, count, start);
    }

    END_FUNC;
}

void SwitchTopo::_cmpt_start_and_count(MPI_Comm comm, const int nBlock[3], int* destRank, int** count, int** start) {
    int size;
    MPI_Comm_size(comm, &size);
    // count the number of blocks to each process
    (*count) = (int*)flups_malloc(size * sizeof(int));
    (*start) = (int*)flups_malloc(size * sizeof(int));
    std::memset((*count), 0, size * sizeof(int));
    std::memset((*start), 0, size * sizeof(int));
    // get the size per block
    const int blockMem = get_blockMemSize();
    // count the number of blocks per rank
    for (int ib = 0; ib < nBlock[0] * nBlock[1] * nBlock[2]; ib++) {
        (*count)[destRank[ib]] += blockMem;
    }
    // compute the start indexes
    if (start != NULL) {
        (*start)[0] = 0;
        for (int ir = 1; ir < size; ir++) {
            (*start)[ir] = (*start)[ir - 1] + (*count)[ir - 1];
        }
    }
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

    // count the number of out edges
    for (int ib = 0; ib < _inBlock[0] * _inBlock[1] * _inBlock[2]; ib++) {
        destsW[_i2o_destRank[ib]] += _iBlockSize[0][ib]*_iBlockSize[1][ib]*_iBlockSize[2][ib];
    }

    // count the number of in edges
    for (int ib = 0; ib < _onBlock[0] * _onBlock[1] * _onBlock[2]; ib++) {
        sourcesW[_o2i_destRank[ib]] += _oBlockSize[0][ib]*_oBlockSize[1][ib]*_oBlockSize[2][ib];
    }

    // Note: by counting the edges like this on every process, we actually obtain
    // twice the number of edges in the total final graph, as the in and out edges 
    // between 2 processes have been accounted by both procs. However, the weight
    // is relative so it doesnt matter.
    END_FUNC;
}