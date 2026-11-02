#include "chunk_tools.hpp"

#include <limits>

/**
 * @brief Decomposes topo_in into chunks, each of them belonging to a different rank in topo_out
 *
 * A chunk is a memory region associated to topo_in that will be existing on a given rank in topo_out.
 * The chunk is given by the starting indexes and the length
 *
 * @param shift location of the (0,0,0) of topo_in in the topo_out
 * @param topo_in the input topology = the home topology of the chunks
 * @param topo_out the output topology, the destination of the chunks (comm and rank depend on topo_out)
 * @param n_chunks the number of chunks created to intersect the topologies
 * @param chunks array of created chunks
 */
void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk** chunks) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // NOT NEEDED as we imposed identical communicators
    FLUPS_INFO("topo-in: nf = %d, axis = %d, local sizes = %d %d %d vs mem size = %d %d %d -- global sizes = %d %d %d",
               topo_in->nf(),
               topo_in->axis(),
               topo_in->nloc(0),
               topo_in->nloc(1),
               topo_in->nloc(2),
               topo_in->nmem(0),
               topo_in->nmem(1),
               topo_in->nmem(2),
               topo_in->nglob(0),
               topo_in->nglob(1),
               topo_in->nglob(2));

    FLUPS_INFO("topo-out: nf = %d, axis = %d, local sizes = %d %d %d vs mem size = %d %d %d -- global sizes = %d %d %d",
               topo_out->nf(),
               topo_out->axis(),
               topo_out->nloc(0),
               topo_out->nloc(1),
               topo_out->nloc(2),
               topo_out->nmem(0),
               topo_out->nmem(1),
               topo_out->nmem(2),
               topo_out->nglob(0),
               topo_out->nglob(1),
               topo_out->nglob(2));

    //--------------------------------------------------------------------------
    /** - get the number of chunks to build */
    //--------------------------------------------------------------------------
    // get the real start and end index that will exist in both topologies
    int topoi_gstart[3], topoi_gend[3];
    topo_in->cmpt_intersect_id(shift, topo_out, topoi_gstart, topoi_gend);
    // Put them in the global reference frame to compute the number of rank spanned by the intersection of topos
    for (int id = 0; id < 3; ++id) {
        topoi_gstart[id] += topo_in->cmpt_start_id(id);
        topoi_gend[id] += topo_in->cmpt_start_id(id);
    }

    // if the current rank has nothing to send (no unknowns in the input topo, or no intersection with the
    // output topo), we have no chunk at all and the rank computations below are not defined
    bool is_empty = false;
    for (int id = 0; id < 3; ++id) {
        is_empty = is_empty || (topoi_gend[id] <= topoi_gstart[id]);
    }
    if (is_empty) {
        FLUPS_INFO("no intersection between the topologies on this rank -> nchunks = 0");
        n_chunks[0] = 0;
        *chunks     = reinterpret_cast<MemChunk*>(m_calloc(0));
        return;
    }

    // upper bound on the number of chunks: the number of ranks spanned by the intersection.
    // some of them might own nothing in the output topology (more ranks than unknowns in a direction),
    // they are skipped below and the actual number of chunks is the number of non-empty ones
    int n_chunks_max = 1;
    int srank[3], erank[3];
    for (int id = 0; id < 3; ++id) {
        // find the rank that will own the starting/ending point of the input topology
        srank[id] = topo_out->cmpt_rank_fromid(topoi_gstart[id] + shift[id], id);
        // the ending point is taken with -1 to be sure to include the end rank
        erank[id] = topo_out->cmpt_rank_fromid(topoi_gend[id] - 1 + shift[id], id);
        // accumulate the number of chunks over the dimensions
        // both ranks are included!
        const int n_ranks = (erank[id] - srank[id] + 1);
        n_chunks_max *= m_max(n_ranks, 0);
    }
    FLUPS_INFO("sranks = %d %d %d - eranks = %d %d %d -> at most nchunks = %d", srank[0], srank[1], srank[2], erank[0], erank[1], erank[2], n_chunks_max);
    *chunks = reinterpret_cast<MemChunk*>(m_calloc(n_chunks_max * sizeof(MemChunk)));

    //--------------------------------------------------------------------------
    /** - Retrieve the id of the current process in the output topology to check
     * if there is a self communication */
    //--------------------------------------------------------------------------
    int out_rank;
    MPI_Comm_rank(topo_out->get_comm(), &out_rank);

    //--------------------------------------------------------------------------
    /** - fill the chunks only if we have some chunks */
    //--------------------------------------------------------------------------
    int chunk_counter = 0;
    for (int ir2 = srank[2]; ir2 <= erank[2]; ++ir2) {
        for (int ir1 = srank[1]; ir1 <= erank[1]; ++ir1) {
            for (int ir0 = srank[0]; ir0 <= erank[0]; ++ir0) {
                // store the current rank in a XYZ format
                const int irank[3]  = {ir0, ir1, ir2};
                const int dest_rank = rankindex(irank, topo_out);
                MemChunk* cchunk    = *chunks + chunk_counter;

                // the destination rank might own nothing in one direction, then there is no chunk to send
                bool is_empty_chunk = false;
                for (int id = 0; id < 3; ++id) {
                    is_empty_chunk = is_empty_chunk || (topo_out->cmpt_start_id_from_rank(irank[id], id) == topo_out->cmpt_start_id_from_rank(irank[id] + 1, id));
                }
                if (is_empty_chunk) {
                    FLUPS_INFO("rank %d owns nothing in the output topology -> no chunk", dest_rank);
                    continue;
                }

                for (int id = 0; id < 3; ++id) {
                    // get the start and end index of the topo OUT and the topo IN in the global reference
                    const int topoo_start = topo_out->cmpt_start_id_from_rank(irank[id], id);
                    const int topoo_end   = topo_out->cmpt_start_id_from_rank(irank[id] + 1, id);

                    FLUPS_INFO("chunk in dim %d topoo_start = %d, topoo_end=%d", id, topoo_start, topoo_end);

                    // Compute the start and end of the chunk in the global frame of reference
                    // make sure that the chunks belongs to the topo_in
                    const int cchunk_gstart = m_max(topoo_start - shift[id], topoi_gstart[id]);
                    const int cchunk_gend   = m_min(topoo_end - shift[id], topoi_gend[id]);

                    // The chunk start in the current topo must be in the local frame of the input topo
                    cchunk->istart[id] = cchunk_gstart - topo_in->cmpt_start_id(id);
                    cchunk->isize[id]  = cchunk_gend - cchunk_gstart;

                    // // we have to now switch the start to the local indexing:
                    // cchunk->istart[id] = cchunk->istart[id];
                    FLUPS_CHECK(cchunk->istart[id] >= 0, "the start id = %d should be >= 0", cchunk->istart[id]);
                    FLUPS_CHECK(cchunk->isize[id] > 0, "The size of the chunk is %d = %d - %d which is not expecte", cchunk->isize[id], cchunk_gend, cchunk_gstart);
                }

                // get the destination rank, no translation is required here as we imposed identical communicators
                cchunk->dest_rank   = dest_rank;
                cchunk->comm        = topo_out->get_comm();
                cchunk->nda         = (size_t)topo_in->lda();
                cchunk->nf          = (size_t)topo_in->nf();
                cchunk->size_padded = get_ChunkPaddedSize(topo_in->nf(), cchunk);

                // fill the axis
                cchunk->axis      = topo_in->axis();
                cchunk->dest_axis = topo_out->axis();

                // setup the offset and the MPI datatype
                const int nmem_in[3] = {topo_in->nmem(0), topo_in->nmem(1), topo_in->nmem(2)};
                ChunkToMPIDataType(nmem_in, cchunk);
                ChunkToDestMPIDataType(cchunk);

                FLUPS_CHECK(topo_in->nf() == topo_out->nf(), "the 2 topo must have matching nfs: %d vs %d", topo_in->nf(), topo_out->nf());
                FLUPS_INFO("chunks going from %d %d %d with size %d %d %d and destination rank %d", cchunk->istart[0], cchunk->istart[1], cchunk->istart[2], cchunk->isize[0], cchunk->isize[1], cchunk->isize[2], cchunk->dest_rank);

                // increasee the chunk counter
                chunk_counter++;
            }
        }
    }

    FLUPS_CHECK(chunk_counter <= n_chunks_max, "the chunk counter = %d must be <= %d", chunk_counter, n_chunks_max);
    FLUPS_CHECK(chunk_counter > 0, "a non-empty intersection must lead to at least one chunk");
    // the actual number of chunks is the number of non-empty ones
    n_chunks[0] = chunk_counter;
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief sets the dtype argument of a chunk, the datatype in the home topology of the chunk
 *
 * @param nmem the memory strides associated to the chunk topo_in
 * @param chunk the memory chunk
 */
void ChunkToMPIDataType(const int nmem[3], MemChunk* chunk) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int nf         = chunk->nf;
    const int ax0        = chunk->axis;  // its home axis
    const int ax[3]      = {ax0, (ax0 + 1) % 3, (ax0 + 2) % 3};
    const int listart[3] = {chunk->istart[ax[0]], chunk->istart[ax[1]], chunk->istart[ax[2]]};

    //..........................................................................
    // the offset is the location in memory for the first dimension
    chunk->offset = localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, 0);
    // get the sride for the different nda in the data
    const size_t offset_dim = localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, 1) - chunk->offset;
    FLUPS_INFO("Offset between each vector component = %zu", offset_dim);

#ifndef NDEBUG
    size_t* offset_check = reinterpret_cast<size_t*>(m_calloc(sizeof(size_t) * chunk->nda));
    for (int ida = 0; ida < chunk->nda; ++ida) {
        offset_check[ida] = localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, ida);
    }
    for (int ida = 1; ida < chunk->nda; ++ida) {
        FLUPS_CHECK(offset_dim == (offset_check[ida] - offset_check[ida - 1]), "the offsets must be =");
    }
    m_free(offset_check);
#endif

    //..........................................................................
    FLUPS_INFO("nmem = %d %d %d", nmem[ax[0]], nmem[ax[1]], nmem[ax[2]]);
    const int    size[4]        = {chunk->isize[ax[0]] * nf, chunk->isize[ax[1]], chunk->isize[ax[2]], chunk->nda};
    const size_t stride_byte[4] = {sizeof(double),
                                   sizeof(double) * nf * nmem[ax[0]],
                                   sizeof(double) * nf * nmem[ax[0]] * nmem[ax[1]],
                                   sizeof(double) * offset_dim};
    // create the 3D datatype
    // MPI_Datatype type_x, type_xy, type_xyz;
    MPI_Datatype type_xy, type_xyz;
    // stride in x = 1, count = chunk size
    // MPI_Type_create_hvector(size[0], 1, stride_byte[0], MPI_DOUBLE, &type_x);
    FLUPS_INFO("puting %d %d-doubles together with strides = %zu", size[0], nf, stride_byte[0]);
    // stride in y = nmem[0]a, count = chunk sie
    // MPI_Type_create_hvector(size[1], 1, stride_byte[1], type_x, &type_xy);
    MPI_Type_create_hvector(size[1], size[0], stride_byte[1], MPI_DOUBLE, &type_xy);
    FLUPS_INFO("puting %d type_x together with strides = %zu", size[1], stride_byte[1]);
    // stride in z = nmem[0]*nmem[1], count = chunk size
    MPI_Type_create_hvector(size[2], 1, stride_byte[2], type_xy, &type_xyz);
    FLUPS_INFO("puting %d type_xy together with strides = %zu", size[2], stride_byte[2]);
    // finally get the different dimensions together
    if (size[3] > 1) {
        MPI_Type_create_hvector(size[3], 1, stride_byte[3], type_xyz, &(chunk->dtype));
    } else {
        MPI_Type_dup(type_xyz, &(chunk->dtype));
    }
    FLUPS_INFO("puting %d type_xyz together with strides = %zu", size[3], stride_byte[3]);

    //..........................................................................
    // commit the new type
    MPI_Type_commit(&(chunk->dtype));

    // free the now useless types
    // MPI_Type_free(&type_x);
    MPI_Type_free(&type_xy);
    MPI_Type_free(&type_xyz);
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief constructs the destination datatype for the given chunk, i.e. its datatype in the recv buffer
 *
 * due to padding in the chunk memory, this datatype is not "just continuous doubles"
 *
 * @param chunk
 */
void ChunkToDestMPIDataType(MemChunk* chunk) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    int      count        = chunk->nda;                                                       // number of blocks
    int      block_length = chunk->isize[0] * chunk->isize[1] * chunk->isize[2] * chunk->nf;  // Number of element per block
    MPI_Aint stride       = chunk->size_padded;                                               // number of bytes between start of each block
    MPI_Type_vector(count, block_length, stride, MPI_DOUBLE, &(chunk->dest_dtype));

    // commit the new type
    MPI_Type_commit(&(chunk->dest_dtype));
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief update the dest_rank and the comm member variable fo the chunk if it belongs to the given comm
 *
 * @param new_comm the new communicator
 * @param new_group the new group associated ot new_comm
 * @param chunk the chunk to update
 * @return true if the chunk belongs to the new group and has been updated
 * @return false if the chunk does not belong to the new comm
 */
void ChunkToNewComm(const MPI_Comm new_comm, const MPI_Group new_group, MemChunk* chunk, bool* is_in_comm) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // get the group associated with the chunk
    MPI_Group chunk_group;
    MPI_Comm_group(chunk->comm, &chunk_group);

    // get the dest_rank of the chunk in the new group
    int new_dest_rank;
    MPI_Group_translate_ranks(chunk_group, 1, &(chunk->dest_rank), new_group, &new_dest_rank);
    is_in_comm[0] = (MPI_UNDEFINED != new_dest_rank);
    if (is_in_comm[0]) {
        chunk->comm      = new_comm;
        chunk->dest_rank = new_dest_rank;
    }
    MPI_Group_free(&chunk_group);
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Prepare the plan for the shuffle for each chunk
 *
 * initilizing mutliple plans will rely on the wisdom of FFTW as soon as no fftw_cleanup is called
 *
 * @param iscomplex indicate if the input topo or the output topo is complex
 * @param chunk the chunk that will store the plan
 */
void PlanShuffleChunk(const bool iscomplex, MemChunk* chunk) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    flups_iodim dims[2];
    // dim[0] = dimension of the targeted FRI (FFTW-convention)
    dims[0].n  = 1;
    dims[0].is = 1;
    dims[0].os = 1;
    // dim[1] = dimension of the current FRI (FFTW-convention)
    dims[1].n  = 1;
    dims[1].is = 1;
    dims[1].os = 1;

    // the chunk arrives in its destination topology and must be shuffled into its correct topology
    // (after reception we are not shuffled yet!)
    int iaxis[3] = {chunk->dest_axis, (chunk->dest_axis + 1) % 3, (chunk->dest_axis + 2) % 3};
    int oaxis[3] = {chunk->axis, (chunk->axis + 1) % 3, (chunk->axis + 2) % 3};

    // note: this is a bit magical and I forgot the exact why we do that
    // loop over the input axis and update the dim[0] and dim[1] data structures accordingly
    for (int id = 0; id < 3; id++) {
        // if the axis we are currently viewing is not the output axis, update the size and the input stride
        if (iaxis[id] != oaxis[0]) {
            // input dimension, change the stride
            dims[0].is = dims[0].is * chunk->isize[iaxis[id]];
            // output dimension change the size
            dims[1].n = dims[1].n * chunk->isize[iaxis[id]];
        } else {
            break;
        }
    }
    // loop over the output axis and update the
    for (int id = 0; id < 3; id++) {
        if (oaxis[id] != iaxis[0]) {
            dims[0].n  = dims[0].n * chunk->isize[oaxis[id]];
            dims[1].os = dims[1].os * chunk->isize[oaxis[id]];
        } else {
            break;
        }
    }
    // display some info
    FLUPS_INFO("shuffle: setting up the shuffle form %d to %d", chunk->axis, chunk->dest_axis);
    FLUPS_INFO("shuffle: iscomplex = %d, blocksize = %d %d %d", iscomplex, chunk->isize[0], chunk->isize[1], chunk->isize[2]);
    FLUPS_INFO("shuffle: DIM 0: n = %d, is=%d, os=%d", dims[0].n, dims[0].is, dims[0].os);
    FLUPS_INFO("shuffle: DIM 1: n = %d, is=%d, os=%d", dims[1].n, dims[1].is, dims[1].os);

    // Allocate temp array to do out-of-place transpose
    chunk->shuffle_buff = (opt_double_ptr)backend_malloc(chunk->nda * chunk->size_padded * sizeof(double));

    // plan the real or complex plan
    // the nf is driven by the OUT topology ALWAYS
    chunk->shuffle = flups_plan_shuffle(dims, chunk->data, chunk->shuffle_buff, iscomplex);
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief executes the shuffle planed by PlanShuffleChunk()
 *
 * @param chunk
 */
void DoShuffleChunk(MemChunk* chunk, stream_t stream) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    FLUPS_CHECK(chunk->has_shuffle, "this chunk has no out-of-place shuffle plan/buffer (axis %d -> %d), see SwitchTopoX::need_shuffle_plan", chunk->axis, chunk->dest_axis);
    for (int ida = 0; ida < chunk->nda; ++ida) {
        opt_double_ptr data_in  = chunk->data + ida * chunk->size_padded;
        opt_double_ptr data_out = chunk->shuffle_buff + ida * chunk->size_padded;

        flups_execute_shuffle(chunk->shuffle, data_in, data_out, 2 == chunk->nf, stream);
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief Copy the memory from the chunk to the data pointer
 *
 * This function uses the memcpy algorithm, which should be the fastest memory copy possible
 * ex of the libc implementation (v2.31)
 *  https://sourceware.org/git/?p=glibc.git;a=blob;f=string/memcpy.c;h=2cb4c76515f476f36a9a8d5dd258ea98e36792b2;hb=9ea3686266dca3f004ba874745a4087a89682617
 *
 * the alignement is automatically performed and exploited, there is not need to do it by hand
 *
 * @param topo the topology in which the chunk and the data are located, must be the input topo of the chunk
 * @param chunk the chunk of memory to copy
 * @param data the vector of data corresponding to the current memory
 */
void CopyChunk2Data(const MemChunk* chunk, const int nmem[3], opt_double_ptr data, bool is_buff, stream_t stream) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // get the current ax as the topo_in one (otherwise the copy doesn't make sense)
    const int nf         = chunk->nf;
    const int ax0        = chunk->axis;
    const int ax[3]      = {ax0, (ax0 + 1) % 3, (ax0 + 2) % 3};
    const int listart[3] = {chunk->istart[ax[0]], chunk->istart[ax[1]], chunk->istart[ax[2]]};

    // get the indexes to copy: n_loop lines of nmax contiguous doubles
    const int n1     = chunk->isize[ax[1]];
    const int n_loop = chunk->isize[ax[1]] * chunk->isize[ax[2]];
    const int nmax   = chunk->isize[ax[0]] * nf;
    // strides (in doubles) of the lines and of the planes, in the chunk and in the data (= localIndex(ax0, 0, 1, 0, ...))
    const size_t s1_chunk = (size_t)nmax;
    const size_t s2_chunk = (size_t)nmax * (size_t)n1;
    const size_t s1_data  = (size_t)nmem[ax[0]] * (size_t)nf;
    const size_t s2_data  = s1_data * (size_t)nmem[ax[1]];

    FLUPS_INFO("copying data at %d %d %d", chunk->istart[0], chunk->istart[1], chunk->istart[2]);

    for (int lia = 0; lia < chunk->nda; ++lia) {
        // get the starting address for the chunk, taking into account the padding
        opt_double_ptr src_data;
        if (is_buff)
            src_data = chunk->shuffle_buff + chunk->size_padded * lia;
        else
            src_data = chunk->data + chunk->size_padded * lia;
        opt_double_ptr trg_data = data + localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, lia);

        // the chunk must be aligned all the time
        FLUPS_CHECK(m_isaligned(src_data, FLUPS_ALIGNMENT), "The chunk memory should be aligned");

        // define the operator: one line per outer index, the point in the line is the inner (contiguous) index
        const auto op = m_backend_lambda(const int io, const int ii) {
            // get the line indexes
            const int i1 = io % n1;
            const int i2 = io / n1;
            // get the starting address of the line
            const double* vsrc = src_data + i1 * s1_chunk + i2 * s2_chunk;
            double*       vtrg = trg_data + i1 * s1_data + i2 * s2_data;
            // copy the data
            vtrg[ii] = vsrc[ii];
        };
        parallel_for_2d(op, n_loop, nmax, "CopyChunk2Data", stream);
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

void CopyData2Chunk(const int nmem[3], const opt_double_ptr data, MemChunk* chunk, stream_t stream) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // get the current ax as the topo_in one (otherwise the copy doesn't make sense)
    const int nf         = chunk->nf;
    const int ax0        = chunk->axis;
    const int ax[3]      = {ax0, (ax0 + 1) % 3, (ax0 + 2) % 3};
    const int listart[3] = {chunk->istart[ax[0]], chunk->istart[ax[1]], chunk->istart[ax[2]]};

    // get the indexes to copy: n_loop lines of nmax contiguous doubles
    const int n1     = chunk->isize[ax[1]];
    const int n_loop = chunk->isize[ax[1]] * chunk->isize[ax[2]];
    const int nmax   = chunk->isize[ax[0]] * nf;
    // strides (in doubles) of the lines and of the planes, in the chunk and in the data (= localIndex(ax0, 0, 1, 0, ...))
    const size_t s1_chunk = (size_t)nmax;
    const size_t s2_chunk = (size_t)nmax * (size_t)n1;
    const size_t s1_data  = (size_t)nmem[ax[0]] * (size_t)nf;
    const size_t s2_data  = s1_data * (size_t)nmem[ax[1]];

    FLUPS_CHECK((chunk->istart[0] + chunk->isize[0]) <= nmem[0], "istart = %d + size = %d must be smaller than the local size %d", chunk->istart[0], chunk->isize[0], nmem[0]);
    FLUPS_CHECK((chunk->istart[1] + chunk->isize[1]) <= nmem[1], "istart = %d + size = %d must be smaller than the local size %d", chunk->istart[1], chunk->isize[1], nmem[1]);
    FLUPS_CHECK((chunk->istart[2] + chunk->isize[2]) <= nmem[2], "istart = %d + size = %d must be smaller than the local size %d", chunk->istart[2], chunk->isize[2], nmem[2]);

    for (int lia = 0; lia < chunk->nda; ++lia) {
        // get the starting address for the chunk, taking into account the padding
        opt_double_ptr trg_data = chunk->data + chunk->size_padded * lia;
        opt_double_ptr src_data = data + localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, lia);

        // we alwas know that the chunk memory is aligned
        FLUPS_CHECK(m_isaligned(trg_data, FLUPS_ALIGNMENT), "The chunk memory should be aligned, size_padded = %ld", chunk->size_padded);
        FLUPS_INFO("pointers are %p and %p", trg_data, src_data);
        FLUPS_INFO("copy %d %d %d from data to chunk", chunk->isize[0], chunk->isize[1], chunk->isize[2]);
        FLUPS_INFO("copy %d bytes in %d loops", nmax * (int)sizeof(double), n_loop);
        FLUPS_INFO("local memory is %d %d %d, listart is %d %d %d", nmem[0], nmem[1], nmem[2], chunk->istart[0], chunk->istart[1], chunk->istart[2]);

        // define the operator: one line per outer index, the point in the line is the inner (contiguous) index
        const auto op = m_backend_lambda(const int io, const int ii) {
            // get the line indexes
            const int i1 = io % n1;
            const int i2 = io / n1;
            // get the starting address of the line
            const double* vsrc = src_data + i1 * s1_data + i2 * s2_data;
            double*       vtrg = trg_data + i1 * s1_chunk + i2 * s2_chunk;
            // copy the data
            vtrg[ii] = vsrc[ii];
        };
        parallel_for_2d(op, n_loop, nmax, "CopyData2chunk", stream);
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief computes the transpose dimensions of the shuffle of a chunk, in nf-units
 *
 * replicates the dims[0].n / dims[1].n computation of PlanShuffleChunk() without touching
 * the flups_shuffle_plan (opaque fftw_plan on the CPU backend); the shuffle then reads
 * out[u] = in[(u % n0) * n1 + u / n0] for u in [0, n0*n1)
 */
static void ChunkShuffleDims(const MemChunk* chunk, int* n0, int* n1) {
    const int iaxis[3] = {chunk->dest_axis, (chunk->dest_axis + 1) % 3, (chunk->dest_axis + 2) % 3};
    const int oaxis[3] = {chunk->axis, (chunk->axis + 1) % 3, (chunk->axis + 2) % 3};
    n0[0]              = 1;
    n1[0]              = 1;
    for (int id = 0; id < 3; id++) {
        if (iaxis[id] != oaxis[0]) {
            n1[0] *= chunk->isize[iaxis[id]];
        } else {
            break;
        }
    }
    for (int id = 0; id < 3; id++) {
        if (oaxis[id] != iaxis[0]) {
            n0[0] *= chunk->isize[oaxis[id]];
        } else {
            break;
        }
    }
}

/**
 * @brief builds the device-visible descriptor array used by the batched copy kernels
 *
 * geometry, chunk pointers, and shuffle dims are all fixed once SwitchTopoX::setup_buffers
 * has run, so the descriptors are built and uploaded once; the field pointer and the nmem
 * (which depends on the real/complex state of the topo) change from one execute to the next
 * and are passed to the kernels instead
 *
 * @param nchunks number of chunks in the array
 * @param chunks the chunk array, with data pointers already assigned
 */
ChunkBatch ChunkBatchSetup(const int nchunks, const MemChunk* chunks) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    ChunkBatch batch;
    if (nchunks <= 0) {
        return batch;
    }

    int ndesc = 0;
    for (int ic = 0; ic < nchunks; ++ic) {
        ndesc += chunks[ic].nda;
    }
    ChunkCopyDesc* h_descs = reinterpret_cast<ChunkCopyDesc*>(m_calloc(ndesc * sizeof(ChunkCopyDesc)));

    size_t ntotal = 0;
    int    id     = 0;
    for (int ic = 0; ic < nchunks; ++ic) {
        const MemChunk* chunk = chunks + ic;
        const int       ax0   = chunk->axis;
        const int       ax[3] = {ax0, (ax0 + 1) % 3, (ax0 + 2) % 3};

        // the kernels assume nf, the axis, and the shuffle flag are uniform over the batch
        FLUPS_CHECK(chunk->nf == chunks[0].nf, "the nf must be uniform over the batch: %d vs %d", chunk->nf, chunks[0].nf);
        FLUPS_CHECK(chunk->axis == chunks[0].axis, "the axis must be uniform over the batch: %d vs %d", chunk->axis, chunks[0].axis);
        FLUPS_CHECK(chunk->dest_axis == chunks[0].dest_axis, "the dest_axis must be uniform over the batch: %d vs %d", chunk->dest_axis, chunks[0].dest_axis);

        int shuffle_n0, shuffle_n1;
        ChunkShuffleDims(chunk, &shuffle_n0, &shuffle_n1);

        // number of transpose tiles of this chunk, for the tiled shuffled unpack
        const int tile    = CHUNK_TRANSPOSE_TILE;
        const int ntile_i = (shuffle_n0 + tile - 1) / tile;
        const int ntile_j = (shuffle_n1 + tile - 1) / tile;

        for (int lia = 0; lia < chunk->nda; ++lia) {
            ChunkCopyDesc* d = h_descs + id;
            d->start         = ntotal;
            d->tile_start    = batch.ntiles;
            d->chunk_data    = chunk->data + chunk->size_padded * lia;
            d->listart[0]    = chunk->istart[ax[0]];
            d->listart[1]    = chunk->istart[ax[1]];
            d->listart[2]    = chunk->istart[ax[2]];
            d->lia           = lia;
            d->nmax          = chunk->isize[ax[0]] * chunk->nf;
            d->n1            = chunk->isize[ax[1]];
            d->shuffle_n0    = shuffle_n0;
            d->shuffle_n1    = shuffle_n1;
            d->ntile_j       = ntile_j;
            d->do_shuffle    = (chunk->axis != chunk->dest_axis);

            FLUPS_CHECK(m_isaligned(d->chunk_data, FLUPS_ALIGNMENT), "The chunk memory should be aligned, size_padded = %ld", chunk->size_padded);

            ntotal += (size_t)(d->nmax) * (size_t)(d->n1) * (size_t)(chunk->isize[ax[2]]);
            batch.ntiles += (size_t)ntile_i * (size_t)ntile_j;
            ++id;
        }
    }
    // the kernel indexing is int-based (parallel_for limitation, same as the r2r plans)
    FLUPS_CHECK(ntotal <= (size_t)std::numeric_limits<int>::max(), "the batched copy is too big for the int-based kernel indexing: %zu", ntotal);

    batch.ndesc      = ndesc;
    batch.ntotal     = ntotal;
    batch.nf         = chunks[0].nf;
    batch.axis       = chunks[0].axis;
    batch.do_shuffle = (chunks[0].axis != chunks[0].dest_axis);

    // upload the descriptors, one-time cost on the default stream
    batch.descs = reinterpret_cast<ChunkCopyDesc*>(backend_malloc(ndesc * sizeof(ChunkCopyDesc)));
    backend_memcpy(reinterpret_cast<double*>(batch.descs), reinterpret_cast<double*>(h_descs), ndesc * sizeof(ChunkCopyDesc), memcpyH2D, (stream_t)(0));
    backend_sync();
    m_free(h_descs);
    //--------------------------------------------------------------------------
    END_FUNC;
    return batch;
}

void ChunkBatchFree(ChunkBatch* batch) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    if (batch->descs != NULL) {
        backend_free(batch->descs);
    }
    batch->descs  = NULL;
    batch->ndesc  = 0;
    batch->ntotal = 0;
    //--------------------------------------------------------------------------
    END_FUNC;
}

/**
 * @brief batched equivalent of CopyData2Chunk over a whole chunk array: one single kernel launch
 *
 * @param batch the descriptors built by ChunkBatchSetup for these chunks
 * @param nmem the memory strides of the home topology of the chunks, at execute time
 * @param data the field array in the home topology of the chunks
 */
void CopyData2Chunk_batch(const ChunkBatch* batch, const int nmem[3], const opt_double_ptr data, stream_t stream) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int            n     = (int)batch->ntotal;
    const int            ndesc = batch->ndesc;
    const int            nf    = batch->nf;
    const int            ax0   = batch->axis;
    const ChunkCopyDesc* descs = batch->descs;
    if (n <= 0) {
        return;
    }
    // field strides in the axis frame, in doubles (captured by value)
    const size_t s1 = (size_t)nf * nmem[ax0];
    const size_t s2 = s1 * nmem[(ax0 + 1) % 3];
    const size_t sa = s2 * nmem[(ax0 + 2) % 3];  // stride between the vector components

    auto op = m_backend_lambda(const int idx) {
        // binary search for the descriptor owning this element
        int lo = 0;
        int hi = ndesc - 1;
        while (lo < hi) {
            const int mid = (lo + hi + 1) / 2;
            if (descs[mid].start <= (size_t)idx)
                lo = mid;
            else
                hi = mid - 1;
        }
        const ChunkCopyDesc d     = descs[lo];
        const int           local = idx - (int)d.start;
        // get the local indexes
        const int i0 = local % d.nmax;
        const int i1 = (local / d.nmax) % d.n1;
        const int i2 = local / (d.nmax * d.n1);
        // the chunk side is contiguous, the field side is strided
        const size_t data_idx = (size_t)(d.listart[0] * nf + i0) + s1 * (d.listart[1] + i1) + s2 * (d.listart[2] + i2) + sa * d.lia;
        d.chunk_data[local]   = data[data_idx];
    };
    parallel_for(op, n, "CopyData2Chunk_batch", stream);
    //--------------------------------------------------------------------------
    END_FUNC;
}

#if M_HAVE_GPU
/**
 * @brief tiled shared-memory version of the shuffled batched unpack (one block = one transpose tile)
 *
 * the naive fused shuffle reads the chunk with a stride of shuffle_n1 between consecutive
 * threads; here a tile is staged through shared memory so that both the chunk reads and the
 * field writes are coalesced. The +1 padding avoids shared memory bank conflicts.
 *
 * templated on NF so a complex (nf=2) unit moves as two adjacent doubles of the same tile row
 */
template <int NF>
__global__ void ChunkTransposeKernel(const ChunkCopyDesc* __restrict__ descs, const int ndesc, double* __restrict__ data,
                                     const size_t s1, const size_t s2, const size_t sa) {
    const int         tile = CHUNK_TRANSPOSE_TILE;
    __shared__ double smem[CHUNK_TRANSPOSE_TILE][CHUNK_TRANSPOSE_TILE * NF + 1];

    // binary search for the descriptor owning this tile
    const size_t bid = blockIdx.x;
    int          lo  = 0;
    int          hi  = ndesc - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (descs[mid].tile_start <= bid)
            lo = mid;
        else
            hi = mid - 1;
    }
    const ChunkCopyDesc d   = descs[lo];
    const int           tid = (int)(bid - d.tile_start);
    const int           sn0 = d.shuffle_n0;
    const int           sn1 = d.shuffle_n1;
    // tile coordinates: bj spans the input-fast dimension (sn1), bi the output-fast one (sn0)
    const int bj = tid % d.ntile_j;
    const int bi = tid / d.ntile_j;

    // load phase: coalesced read of the chunk, in[i * sn1 + j] with j the fast index
    for (int ty = threadIdx.y; ty < tile; ty += blockDim.y) {
        const int i = bi * tile + ty;
        const int j = bj * tile + threadIdx.x;
        if (i < sn0 && j < sn1) {
#pragma unroll
            for (int c = 0; c < NF; ++c) {
                smem[ty][threadIdx.x * NF + c] = d.chunk_data[(size_t)(i * sn1 + j) * NF + c];
            }
        }
    }
    __syncthreads();

    // write phase: coalesced write to the field, out flat index u = i + j * sn0 with i fast
    for (int ty = threadIdx.y; ty < tile; ty += blockDim.y) {
        const int i = bi * tile + threadIdx.x;
        const int j = bj * tile + ty;
        if (i < sn0 && j < sn1) {
            const int u = i + j * sn0;
#pragma unroll
            for (int c = 0; c < NF; ++c) {
                const int    local    = u * NF + c;
                const int    i0       = local % d.nmax;
                const int    i1       = (local / d.nmax) % d.n1;
                const int    i2       = local / (d.nmax * d.n1);
                const size_t data_idx = (size_t)(d.listart[0] * NF + i0) + s1 * (d.listart[1] + i1) + s2 * (d.listart[2] + i2) + sa * d.lia;
                data[data_idx]        = smem[threadIdx.x][ty * NF + c];
            }
        }
    }
}

/**
 * @brief launches the tiled transpose kernel, one block of (tile x 8) threads per tile
 */
static void LaunchChunkTranspose(const ChunkBatch* batch, const size_t s1, const size_t s2, const size_t sa,
                                 opt_double_ptr data, stream_t stream) {
    const dim3 grid((unsigned int)batch->ntiles);
    const dim3 block(CHUNK_TRANSPOSE_TILE, 8);

    const ChunkCopyDesc* descs  = batch->descs;
    const int            ndesc  = batch->ndesc;
    const void*          args[] = {&descs, &ndesc, &data, &s1, &s2, &sa};
#if defined(BACKEND_CUDA)
    void* kernel = (batch->nf == 2) ? (void*)ChunkTransposeKernel<2> : (void*)ChunkTransposeKernel<1>;
    m_cudart_check_errors(cudaLaunchKernel(kernel, grid, block, const_cast<void**>(args), 0, stream));
#elif defined(BACKEND_HIP)
    void* kernel = (batch->nf == 2) ? (void*)ChunkTransposeKernel<2> : (void*)ChunkTransposeKernel<1>;
    m_hiprt_check_errors(hipLaunchKernel(kernel, grid, block, const_cast<void**>(args), 0, stream));
#endif
}
#endif  // M_HAVE_GPU

/**
 * @brief batched equivalent of DoShuffleChunk + CopyChunk2Data over a whole chunk array
 *
 * the shuffle (a per-chunk 2D transpose) is fused into the read indexing, so the
 * shuffle_buff is never touched and a single kernel launch unpacks all the chunks
 *
 * @param batch the descriptors built by ChunkBatchSetup for these chunks
 * @param nmem the memory strides of the home topology of the chunks, at execute time
 * @param data the field array in the home topology of the chunks
 */
void CopyChunk2Data_batch(const ChunkBatch* batch, const int nmem[3], opt_double_ptr data, stream_t stream) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    const int            n     = (int)batch->ntotal;
    const int            ndesc = batch->ndesc;
    const int            nf    = batch->nf;
    const int            ax0   = batch->axis;
    const ChunkCopyDesc* descs = batch->descs;
    if (n <= 0) {
        return;
    }
    // field strides in the axis frame, in doubles (captured by value)
    const size_t s1 = (size_t)nf * nmem[ax0];
    const size_t s2 = s1 * nmem[(ax0 + 1) % 3];
    const size_t sa = s2 * nmem[(ax0 + 2) % 3];  // stride between the vector components

#if M_HAVE_GPU
    // shuffled batches go through the tiled transpose so both sides stay coalesced
    if (batch->do_shuffle) {
        LaunchChunkTranspose(batch, s1, s2, sa, data, stream);
        return;
    }
#endif

    auto op = m_backend_lambda(const int idx) {
        // binary search for the descriptor owning this element
        int lo = 0;
        int hi = ndesc - 1;
        while (lo < hi) {
            const int mid = (lo + hi + 1) / 2;
            if (descs[mid].start <= (size_t)idx)
                lo = mid;
            else
                hi = mid - 1;
        }
        const ChunkCopyDesc d     = descs[lo];
        const int           local = idx - (int)d.start;
        // get the local indexes
        const int i0 = local % d.nmax;
        const int i1 = (local / d.nmax) % d.n1;
        const int i2 = local / (d.nmax * d.n1);
        // fused shuffle = transposed read of the chunk, in nf-units
        int src = local;
        if (d.do_shuffle) {
            const int u = local / nf;
            const int c = local - u * nf;
            src         = ((u % d.shuffle_n0) * d.shuffle_n1 + u / d.shuffle_n0) * nf + c;
        }
        const size_t data_idx = (size_t)(d.listart[0] * nf + i0) + s1 * (d.listart[1] + i1) + s2 * (d.listart[2] + i2) + sa * d.lia;
        data[data_idx]        = d.chunk_data[src];
    };
    parallel_for(op, n, "CopyChunk2Data_batch", stream);
    //--------------------------------------------------------------------------
    END_FUNC;
}
