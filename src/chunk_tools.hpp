#ifndef CHUNKTOOLS_HPP_
#define CHUNKTOOLS_HPP_

#include "Topology.hpp"
#include "defines.hpp"

typedef struct
{
    int istart[3];  //!< local start index for the memory chunk (012 indexing)
    int isize[3];   // !< local size for the memory chunks (012 indexing)
    int dest_rank;  //!< destination ranks
} MemChunk;

void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk* chunks);

void void PlanShuffleChunk(const bool iscomplex, MemChunk* chunk, double* data);

/**
 * @brief returns the memory size (padded to a multiple of alignment) of a MemChunk
 *
 * We want to ensure that the next block will begin in an aligned position!
 * Therefore we return the padded size of the current chunk
 *
 * @param nf
 * @param chunk
 * @return size_t
 */
inline size_t get_ChunkMemSize(const int nf, const MemChunk* chunk) {
    //----------------------------------------------------------------------
    const size_t total     = (size_t)(chunk->isize[0]) * (size_t)(chunk->isize[1]) * (size_t)(chunk->isize[2]) * (size_t)(nf);
    const size_t total_ext = total + (M_ALIGNMENT - 1);
    return total_ext - (total_ext % M_ALIGNMENT);
    //----------------------------------------------------------------------
}

#endif