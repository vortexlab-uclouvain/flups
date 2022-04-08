#ifndef CHUNKTOOLS_HPP_
#define CHUNKTOOLS_HPP_

#include "Topology.hpp"
#include "defines.hpp"

/**
 * @brief A "chunk" is a memory block belonging to an input topology. The block is sent over to the output topology and shuffled
 *
 */
typedef struct
{
    int iaxis;  //!< the principal axis in the input topology
    int oaxis;  //!< the principal axis in the output topology

    int istart[3];  //!< local start index for the memory chunk (012 indexing)
    int isize[3];   // !< local size for the memory chunks (012 indexing)

    int dest_rank;  //!< destination ranks in the output topology

    fftw_plan shuffle;  //!< the shuffle plan used by FFTW to reorder data

    size_t         nda;          //!< the number of data array (1 if scalar, 3 if vector)
    size_t         size_padded;  //!< padded size for the data ptr
    opt_double_ptr data;         //!< the pointer to the data in the buffer

} MemChunk;

void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk* chunks);
void PlanShuffleChunk(const bool iscomplex, MemChunk* chunk);

/**
 * @brief returns the memory size (padded to a multiple of alignment) of a MemChunk
 *
 * We want to ensure that the next block will begin in an aligned position!
 * Therefore we return the padded size of the current chunk
 *
 * @param nf the number of field, = 1 if real, = 2 if complex
 * @param chunk the chunk of which the padded size is requested
 * @return size_t
 */
inline size_t get_ChunkPaddedSize(const size_t nf, const MemChunk* chunk) {
    //----------------------------------------------------------------------
    const size_t total     = (size_t)(chunk->isize[0]) * (size_t)(chunk->isize[1]) * (size_t)(chunk->isize[2]) * nf;
    const size_t total_ext = total + (M_ALIGNMENT - 1);
    return total_ext - (total_ext % M_ALIGNMENT);
    //----------------------------------------------------------------------
}

#endif