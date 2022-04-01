/**
 * @file SwitchTopo.hpp
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

#ifndef SWITCHTOPO_HPP
#define SWITCHTOPO_HPP

#include <cstring>
#include "Topology.hpp"
#include "defines.hpp"
#include "mpi.h"
#include "omp.h"
// #include "h3lpr/profiler.hpp"

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
class SwitchTopo {
   protected:
    MPI_Comm inComm_ = NULL; /**<@brief the reference input communicator */
    MPI_Comm outComm_ = NULL; /**<@brief the reference output communicator */
    MPI_Comm subcomm_ = NULL; /**<@brief the subcomm for this switchTopo */
    int shift_[3]; /**<@brief the shift in memory */

    int inBlock_; /**<@brief the local number of block in each dim in the input topology */
    int onBlock_; /**<@brief the local number of block in each dim in the output topology  */

    int* iBlockiStart_[3] = {NULL, NULL, NULL}; /**<@brief the local starting index for a block in the input topo  */
    int* oBlockiStart_[3] = {NULL, NULL, NULL}; /**<@brief the local starting index for a block in the output topo  */

    int* iBlockSize_[3] = {NULL, NULL, NULL}; /**<@brief The number of data per blocks in each dim for each block (!same on each process! and 012-indexing)  */
    int* oBlockSize_[3] = {NULL, NULL, NULL}; /**<@brief The number of data per blocks in each dim for each block (!same on each process! and 012-indexing)  */

    int* i2o_destRank_ = NULL; /**<@brief The destination rank in the output topo of each block */
    int* o2i_destRank_ = NULL; /**<@brief The destination rank in the output topo of each block */

    const Topology *topo_in_  = NULL; /**<@brief input topology  */
    const Topology *topo_out_ = NULL; /**<@brief  output topology */

    opt_double_ptr *sendBuf_ = NULL; /**<@brief The send buffer for MPI send */
    opt_double_ptr *recvBuf_ = NULL; /**<@brief The recv buffer for MPI recv */

    fftw_plan* i2o_shuffle_ = NULL;
    fftw_plan* o2i_shuffle_ = NULL;

    H3LPR::Profiler* prof_    = NULL;
    int       iswitch_ = -1;

   public:
    virtual ~SwitchTopo() {};
    virtual void setup()                                                                    = 0;
    virtual void setup_buffers(opt_double_ptr sendData, opt_double_ptr recvData)            = 0;
    virtual void execute(opt_double_ptr v, const int sign) const                            = 0;
    virtual void disp() const                                                               = 0;

    /**
     * @brief Get the memory size of a block padded to ensure alignment
     *
     * The size the size of one block, for 1 component of lda
     * @warning
     * Since we use gathered blocks, it is NOT STRAIGHTFORWARD to impose a common size for every block on every proc.
     * Therefore, we chose not to do it!!
     * 
     * @param ib the block id
     * @param nf the number of fields inside an element
     * @param blockSize the number of element in 3d for each block
     * @return size_t the padded size in memory of one block
     */
    inline size_t get_blockMemSize(const int ib, const int nf,const int * const blockSize[3]) const {
        // get the in and out sizes
        size_t total = (size_t)(blockSize[0][ib]) * (size_t)(blockSize[1][ib]) * (size_t)(blockSize[2][ib]) * (size_t)(nf);
        // add the difference with the alignement to be always aligned
        size_t alignDelta = ((total * sizeof(double)) % FLUPS_ALIGNMENT == 0) ? 0 : (FLUPS_ALIGNMENT - (total * sizeof(double)) % FLUPS_ALIGNMENT) / sizeof(double);
        total             = total + alignDelta;
        FLUPS_CHECK((total * sizeof(double)) % FLUPS_ALIGNMENT == 0, "The total size of one block HAS to match the alignement size");
        // return the total size
        return total;
    };
    /**
     * @brief return the buffer size for one proc = number of blocks * blocks memory size * lda component
     * 
     * @return size_t 
     */
    inline size_t get_bufMemSize() const {
        // the nf is the maximum between in and out
        const int nf = std::max(topo_in_->nf(),topo_out_->nf());
        // nultiply by the number of blocks
        size_t total = 0;
        for(int ib=0; ib<inBlock_; ib++){
            total += get_blockMemSize(ib,nf,iBlockSize_) * ((size_t)topo_in_->lda());
        }
        for(int ib=0; ib<onBlock_; ib++){
            total += get_blockMemSize(ib,nf,oBlockSize_) * ((size_t)topo_in_->lda());
        }
        // return the total size
        return total;
    };

    void add_toGraph(int* sourcesW, int* destsW) const;

   protected:
    void cmpt_nByBlock_(int istart[3], int iend[3], int ostart[3], int oend[3],int nByBlock[3]);
    void cmpt_blockDestRank_(const int nBlock[3], const int nByBlock[3], const int shift[3], const int istart[3], const Topology* topo_in, const Topology* topo_out, int* destRank);
    void cmpt_blockIndexes_(const int istart[3], const int iend[3], const int nByBlock[3], const Topology *topo,int nBlock[3]);

    void cmpt_commSplit_();
    void setup_subComm_(const int nBlock, const int lda, int* blockSize[3], int* destRank, int** count, int** start);
    void cmpt_start_and_count_(MPI_Comm comm, const int nBlock, const int lda, int* blockSize[3], int* destRank, int** count, int** start);
    void setup_shuffle_(const int bSize[3], const Topology* topo_in, const Topology* topo_out, double* data, fftw_plan* shuffle);
    void gather_blocks_(const Topology* topo, int nByBlock[3], int istart[3],int iend[3], int nBlockv[3], int* blockSize[3], int* blockiStart[3], int* nBlock, int** destRank);
    void gather_tags_(MPI_Comm comm, const int inBlock, const int onBlock, const int* i2o_destRank, const int* o2i_destRank, int** i2o_destTag, int** o2i_destTag);
};

static inline int gcd(int a, int b) {
    return (a == 0) ? b : gcd(b % a, a);
}

/**
 * @brief translate a list of ranks of size size from inComm to outComm
 * 
 * ranks are replaced with their new values.
 * 
 * @param size 
 * @param ranks a list of ranks expressed in the inComm
 * @param inComm input communicator
 * @param outComm output communicator
 */
inline static void translate_ranks(int size, int* ranks, MPI_Comm inComm, MPI_Comm outComm) {
    BEGIN_FUNC;

    int comp;
    MPI_Comm_compare(inComm, outComm, &comp);
    FLUPS_CHECK(size!=0,"size cant be 0.");

    int* tmprnks = (int*) m_calloc(size*sizeof(int));
    std::memcpy(tmprnks,ranks,size*sizeof(int));

    int err;
    if (comp != MPI_IDENT) {
        MPI_Group group_in, group_out;
        err = MPI_Comm_group(inComm, &group_in);
        FLUPS_CHECK(err==MPI_SUCCESS,"wrong group in");
        err = MPI_Comm_group(outComm, &group_out);
        FLUPS_CHECK(err==MPI_SUCCESS,"wrong group out");

        err = MPI_Group_translate_ranks(group_in, size, tmprnks, group_out, ranks);
        FLUPS_CHECK(err == MPI_SUCCESS, "Could not find a correspondance between incomm and outcomm.");
    }

    m_free(tmprnks);
    END_FUNC;
}

#endif