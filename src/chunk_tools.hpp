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
    int axis;       //!< the principal axis in the "input" topology, where the chunk is defined
    int istart[3];  //!< local start index for the memory chunk (012 indexing), in the "input" topology
    int isize[3];   // !< local size for the memory chunks (012 indexing), in the "input" topology

    int      dest_rank;  //!< destination ranks in the destination topology
    int      dest_axis;  //!< the principal axis in the destination topology
    MPI_Comm comm;       //!< the communicator to be used for the communication, also dictates the dest_rank id

    fftw_plan shuffle;  //!< the shuffle plan used by FFTW to reorder data

    size_t       offset;  //!< offset in memory in the "input" topology
    MPI_Datatype dtype;   //!< datatype in the "input" topology

    int            nda;          //!< the number of data array (1 if scalar, 3 if vector)
    int            nf;           //!< the number of double per data (1 if real, 2 if complex)
    size_t         size_padded;  //!< padded size for the data ptr
    opt_double_ptr data;         //!< the pointer to the data in the buffer

} MemChunk;

void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk** chunks, int* self_comm);
void PlanShuffleChunk(const bool iscomplex, MemChunk* chunk);
void DoShuffleChunk(MemChunk* chunk);

void CopyChunk2Data(const MemChunk* chunk, const int nmem[3], opt_double_ptr data);
void CopyData2Chunk(const int nmem[3], const opt_double_ptr data, MemChunk* chunk);

void ChunkToMPIDataType(const int nmem[3], MemChunk* chunk);//, size_t* offset, MPI_Datatype* type_xyzd);

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
    FLUPS_CHECK(M_ALIGNMENT % sizeof(double) == 0, "The alignement %d must be a multiple of %zu", M_ALIGNMENT, sizeof(double));
    //----------------------------------------------------------------------
    const size_t total     = (size_t)(chunk->isize[0]) * (size_t)(chunk->isize[1]) * (size_t)(chunk->isize[2]) * nf;
    const size_t align     = M_ALIGNMENT / sizeof(double);
    const size_t total_ext = total + (align - 1);

    return total_ext - (total_ext % align);
    //----------------------------------------------------------------------
}

#endif