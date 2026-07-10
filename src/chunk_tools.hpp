#ifndef CHUNKTOOLS_HPP_
#define CHUNKTOOLS_HPP_

#include "Topology.hpp"
#include "backend.hpp"
#include "defines.hpp"

/**
 * @brief A "chunk" is a memory block belonging to an input topology. The block is sent over to the output topology and shuffled
 *
 */
struct MemChunk {
    int axis;       //!< the principal axis in the "input" topology, where the chunk is defined
    int istart[3];  //!< local start index for the memory chunk (012 indexing), in the "input" topology
    int isize[3];   // !< local size for the memory chunks (012 indexing), in the "input" topology

    int      dest_rank;  //!< destination ranks in the destination topology
    int      dest_axis;  //!< the principal axis in the destination topology
    MPI_Comm comm;       //!< the communicator to be used for the communication, also dictates the dest_rank id

    flups_shuffle_plan shuffle;       //!< the shuffle plan used by FFTW to reorder data
    opt_double_ptr     shuffle_buff;  //!< the array for the out-of-place shuffle

    size_t       offset;  //!< offset in memory in the "input" topology
    MPI_Datatype dtype;   //!< datatype in the "input" topology

    MPI_Datatype dest_dtype;  //!< datatype in the "output" topology

    int            nda;          //!< the number of data array (1 if scalar, 3 if vector)
    int            nf;           //!< the number of double per data (1 if real, 2 if complex)
    size_t         size_padded;  //!< padded size for the data ptr
    opt_double_ptr data;         //!< the pointer to the data in the buffer

    ~MemChunk() {
        MPI_Type_free(&dtype);
        MPI_Type_free(&dest_dtype);
    }
};

void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk** chunks);

void ChunkToNewComm(const MPI_Comm new_comm, const MPI_Group new_group, MemChunk* chunk, bool* is_in_comm);
void PlanShuffleChunk(const bool iscomplex, MemChunk* chunk);
void DoShuffleChunk(MemChunk* chunk, stream_t stream);

void CopyChunk2Data(const MemChunk* chunk, const int nmem[3], opt_double_ptr data, bool is_buff, stream_t stream);
void CopyData2Chunk(const int nmem[3], const opt_double_ptr data, MemChunk* chunk, stream_t stream);

/**
 * @brief flat descriptor of one (chunk, lia) copy, precomputed on the host for the batched kernels
 *
 * the chunk-side memory is contiguous per (chunk, lia): element (i0,i1,i2) in the axis frame
 * sits at the local linear index i0 + nmax*(i1 + n1*i2). The field side is NOT precomputed:
 * the topologies switch between real and complex states over the solve, so nmem is only
 * known at execute time (like the base copy functions) and the kernels get it as an argument
 */
struct ChunkCopyDesc {
    size_t  start;       //!< first global element index of this (chunk, lia) in the batch
    size_t  tile_start;  //!< first global tile index of this (chunk, lia), for the tiled transpose kernel
    double* chunk_data;  //!< chunk->data + lia * size_padded
    int     listart[3];  //!< chunk->istart in the axis frame (ax0, ax1, ax2)
    int     lia;         //!< the component of this descriptor
    int     nmax;        //!< isize[ax0] * nf = length of a contiguous line (in doubles)
    int     n1;          //!< isize[ax1]
    int     shuffle_n0;  //!< transpose dims in nf-units (= dims[0].n of PlanShuffleChunk)
    int     shuffle_n1;  //!< transpose dims in nf-units (= dims[1].n of PlanShuffleChunk)
    int     ntile_j;     //!< number of transpose tiles along the shuffle_n1 dimension
    int     do_shuffle;  //!< true if the chunk must be shuffled upon reception
};

//! tile size of the shared-memory transpose used by the shuffled batched unpack
#define CHUNK_TRANSPOSE_TILE 32

/**
 * @brief a batch of (chunk, lia) copy descriptors, built once at setup and device-visible
 */
struct ChunkBatch {
    int            ndesc      = 0;     //!< number of (chunk, lia) descriptors
    size_t         ntotal     = 0;     //!< total number of doubles over all descriptors
    size_t         ntiles     = 0;     //!< total number of transpose tiles over all descriptors
    int            nf         = 1;     //!< nf of the chunks, uniform over the batch
    int            axis       = 0;     //!< home axis of the chunks, uniform over the batch
    int            do_shuffle = 0;     //!< shuffle flag, uniform over the batch
    ChunkCopyDesc* descs      = NULL;  //!< device-visible descriptor array
};

ChunkBatch ChunkBatchSetup(const int nchunks, const MemChunk* chunks);
void       ChunkBatchFree(ChunkBatch* batch);

void CopyData2Chunk_batch(const ChunkBatch* batch, const int nmem[3], const opt_double_ptr data, stream_t stream);
void CopyChunk2Data_batch(const ChunkBatch* batch, const int nmem[3], opt_double_ptr data, stream_t stream);

void ChunkToMPIDataType(const int nmem[3], MemChunk* chunk);  //, size_t* offset, MPI_Datatype* type_xyzd);
void ChunkToDestMPIDataType(MemChunk* chunk);

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
    FLUPS_CHECK(FLUPS_ALIGNMENT % sizeof(double) == 0, "The alignement %d must be a multiple of %zu", FLUPS_ALIGNMENT, sizeof(double));
    //----------------------------------------------------------------------
    const size_t total     = (size_t)(chunk->isize[0]) * (size_t)(chunk->isize[1]) * (size_t)(chunk->isize[2]) * nf;
    const size_t align     = FLUPS_ALIGNMENT / sizeof(double);
    const size_t total_ext = total + (align - 1);

    return total_ext - (total_ext % align);
    //----------------------------------------------------------------------
}

#endif
