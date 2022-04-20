#include "chunk_tools.hpp"

void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk **chunks) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // NOT NEEDED as we imposed identical communicators
    /*
    // get the communicators and the associated groups
    MPI_Comm in_comm = topo_in->get_comm();
    MPI_Comm out_comm = topo_out->get_comm();
    MPI_Group group_in, group_out;
    err = MPI_Comm_group(inComm, &group_in);
    FLUPS_CHECK(err == MPI_SUCCESS, "wrong group in");
    err = MPI_Comm_group(outComm, &group_out);
    FLUPS_CHECK(err == MPI_SUCCESS, "wrong group out");
    */

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
    for(int id = 0; id < 3; ++id){
        topoi_gstart[id] += topo_in->cmpt_start_id(id);
        topoi_gend[id]   += topo_in->cmpt_start_id(id);
    }

    n_chunks[0] = 1;
    int srank[3], erank[3];
    for (int id = 0; id < 3; ++id) {
        // find the rank that will own the starting/ending point of the input topology
        // the starting point is always properly defined
        srank[id] = topo_out->cmpt_rank_fromid(topoi_gstart[id] + shift[id], id);
        // the ending point is taken with -1 to be sure to include the end rank
        // unless start = end, then we just use the start id
        if (topoi_gend[id] > topoi_gstart[id]) {
            erank[id] = topo_out->cmpt_rank_fromid(topoi_gend[id] - 1 + shift[id], id);
        } else {
            erank[id] = srank[id] - 1;
        }
        // accumulate the number of chunks over the dimensions
        // both ranks are included!
        const int n_ranks = (erank[id] - srank[id] + 1);
        n_chunks[0] *= m_max(n_ranks, 0);
    }
    FLUPS_INFO("sranks = %d %d %d - eranks = %d %d %d -> nchunks = %d", srank[0], srank[1], srank[2], erank[0], erank[1], erank[2], n_chunks[0]);
    *chunks = reinterpret_cast<MemChunk*>(m_calloc(n_chunks[0] * sizeof(MemChunk)));

    //--------------------------------------------------------------------------
    /** - fill the chunks only if we have some chunks */
    //--------------------------------------------------------------------------
    int chunk_counter = 0;
    for (int ir2 = srank[2]; ir2 <= erank[2]; ++ir2) {
        for (int ir1 = srank[1]; ir1 <= erank[1]; ++ir1) {
            for (int ir0 = srank[0]; ir0 <= erank[0]; ++ir0) {
                // get the current chunk
                MemChunk* cchunk = *chunks + chunk_counter;
                // store the current rank in a XYZ format
                const int irank[3] = {ir0, ir1, ir2};
                for (int id = 0; id < 3; ++id) {
                    // get the start and end index of the topo OUT and the topo IN in the global reference
                    const int topoo_start = topo_out->cmpt_start_id_from_rank(irank[id], id);
                    const int topoo_end   = topo_out->cmpt_start_id_from_rank(irank[id] + 1, id);

                    FLUPS_INFO("chunk in dim %d topoo_start = %d, topoo_end=%d", id, topoo_start, topoo_end);
                    
                    // Compute the start and end of the chunk in the global frame of reference
                    // make sure that the chunks belongs to the topo_in
                    const int cchunk_gstart = m_max(topoo_start - shift[id], topoi_gstart[id]);
                    const int cchunk_gend = m_min(topoo_end - shift[id], topoi_gend[id]);
                    
                    // The chunk start in the current topo must be in the local frame of the input topo
                    cchunk->istart[id] = cchunk_gstart - topo_in->cmpt_start_id(id);
                    cchunk->isize[id]  = cchunk_gend - cchunk_gstart; 

                    // // we have to now switch the start to the local indexing:
                    // cchunk->istart[id] = cchunk->istart[id];
                    FLUPS_CHECK(cchunk->istart[id] >= 0, "the start id = %d should be >= 0", cchunk->istart[id]);
                    FLUPS_CHECK(cchunk->isize[id] > 0, "The size of the chunk is %d which is not expecte", cchunk->isize[id]);
                }
                // get the destination rank, no translation is required here as we imposed identical communicators
                cchunk->dest_rank   = rankindex(irank, topo_out);
                cchunk->nda         = (size_t)topo_in->lda();
                cchunk->nf          = (size_t)topo_in->nf();
                cchunk->size_padded = get_ChunkPaddedSize(topo_in->nf(), cchunk);

                // fill the axis
                cchunk->axis      = topo_in->axis();
                cchunk->dest_axis = topo_out->axis();

                FLUPS_CHECK(topo_in->nf() == topo_out->nf(), "the 2 topo must have matching nfs: %d vs %d", topo_in->nf(), topo_out->nf());
                FLUPS_INFO("chunks going from %d %d %d with size %d %d %d and destination rank %d", cchunk->istart[0], cchunk->istart[1], cchunk->istart[2], cchunk->isize[0], cchunk->isize[1], cchunk->isize[2], cchunk->dest_rank);

                // increasee the chunk counter
                chunk_counter += 1;
            }
        }
    }
    FLUPS_CHECK(chunk_counter == n_chunks[0], "the chunk counter = %d must be = %d", chunk_counter, n_chunks[0]);
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

    // plan the real or complex plan
    // the nf is driven by the OUT topology ALWAYS
    if (!iscomplex) {
        chunk->shuffle = fftw_plan_guru_r2r(0, NULL, 2, dims, chunk->data, chunk->data, NULL, FFTW_FLAG);
        // FLUPS_CHECK(chunk->shuffle != NULL, "Plan has not been setup");
    } else {
        chunk->shuffle = fftw_plan_guru_dft(0, NULL, 2, dims, (fftw_complex*)chunk->data, (fftw_complex*)chunk->data, FLUPS_FORWARD, FFTW_FLAG);
        // FLUPS_CHECK(chunk->shuffle != NULL, "Plan has not been setup");
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

void DoShuffleChunk(MemChunk* chunk) {
    BEGIN_FUNC;
    //--------------------------------------------------------------------------
    // only the master call the fftw_execute which is executed in multithreading
    for (int ida = 0; ida < chunk->nda; ++ida) {
        opt_double_ptr data_ptr = chunk->data + ida * chunk->size_padded;
        if (chunk->nf == 1) {
            // we execute it on the different directions of the field but the memory has the same properties (alignement, lenght, etc)
            fftw_execute_r2r(chunk->shuffle, data_ptr, data_ptr);
        } else {
            fftw_execute_dft(chunk->shuffle, (opt_complex_ptr)(data_ptr), (opt_complex_ptr)(data_ptr));
        }
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
void CopyChunk2Data(const MemChunk* chunk, const int nmem[3], opt_double_ptr data) {
    BEGIN_FUNC;
    FLUPS_CHECK(FLUPS_ALIGNMENT == M_ALIGNMENT, "This is only temporary, the alignement should not be in H3LPR");
    //--------------------------------------------------------------------------
    // get the current ax as the topo_in one (otherwise the copy doesn't make sense)
    const int nf         = chunk->nf;
    const int ax0        = chunk->axis;
    const int ax[3]      = {ax0, (ax0 + 1) % 3, (ax0 + 2) % 3};
    const int listart[3] = {chunk->istart[ax[0]], chunk->istart[ax[1]], chunk->istart[ax[2]]};

    // get the indexes to copy
    const size_t n_loop    = chunk->isize[ax[1]] * chunk->isize[ax[2]];
    const size_t nmax_byte = chunk->isize[ax[0]] * nf * sizeof(double);

    FLUPS_INFO("copying data at %d %d %d", chunk->istart[0], chunk->istart[1], chunk->istart[2]);
    // FLUPS_INFO("nf = %d",nf);
    //     if (nf == 2) {
    //         for (int id2 = 0; id2 < chunk->isize[ax[2]]; ++id2) {
    //             for (int id1 = 0; id1 < chunk->isize[ax[1]]; ++id1) {
    //                 for (int id0 = 0; id0 < chunk->isize[ax[0]]; ++id0) {
    //                     const int id = localIndex(ax[0], id0, id1, id2, ax0, chunk->isize, nf, 0);
    //                     printf("(%f %f)", chunk->data[id], chunk->data[id + 1]);
    //                 }
    //                 printf("\n");
    //             }
    //         }
    //     }

#pragma omp parallel proc_bind(close)
    for (int lia = 0; lia < chunk->nda; ++lia) {
        // get the starting address for the chunk, taking into account the padding
        opt_double_ptr src_data = chunk->data + chunk->size_padded * lia;
        opt_double_ptr trg_data = data + localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, lia);

        // the chunk must be aligned all the time
        FLUPS_CHECK(m_isaligned(src_data), "The chunk memory should be aligned");

#pragma omp for schedule(static)
        for (int il = 0; il < n_loop; ++il) {
            // get the local indexes (we cannot used the collaspedIndex one!!!)
            const int i2 = il / (chunk->isize[ax[1]]);
            const int i1 = il % (chunk->isize[ax[1]]);
            // get the starting adddress for the memcpy
            const double* __restrict vsrc = src_data + localIndex(ax0, 0, i1, i2, ax0, chunk->isize, nf, 0);
            double* __restrict vtrg       = trg_data + localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);
            memcpy(vtrg, vsrc, nmax_byte);
        }
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

void CopyData2Chunk(const int nmem[3], const opt_double_ptr data, MemChunk* chunk) {
    BEGIN_FUNC;
    FLUPS_CHECK(FLUPS_ALIGNMENT == M_ALIGNMENT, "This is only temporary, the alignement should not be in H3LPR");
    //--------------------------------------------------------------------------
    // get the current ax as the topo_in one (otherwise the copy doesn't make sense)
    const int nf         = chunk->nf;
    const int ax0        = chunk->axis;
    const int ax[3]      = {ax0, (ax0 + 1) % 3, (ax0 + 2) % 3};
    const int listart[3] = {chunk->istart[ax[0]], chunk->istart[ax[1]], chunk->istart[ax[2]]};

    // get the indexes to copy
    const size_t n_loop    = chunk->isize[ax[1]] * chunk->isize[ax[2]];
    const size_t nmax_byte = chunk->isize[ax[0]] * nf * sizeof(double);

    FLUPS_CHECK((chunk->istart[0] + chunk->isize[0]) <= nmem[0],"istart = %d + size = %d must be smaller than the local size %d",chunk->istart[0], chunk->isize[0], nmem[0]);
    FLUPS_CHECK((chunk->istart[1] + chunk->isize[1]) <= nmem[1],"istart = %d + size = %d must be smaller than the local size %d",chunk->istart[1], chunk->isize[1], nmem[1]);
    FLUPS_CHECK((chunk->istart[2] + chunk->isize[2]) <= nmem[2],"istart = %d + size = %d must be smaller than the local size %d",chunk->istart[2], chunk->isize[2], nmem[2]);

#pragma omp parallel proc_bind(close)
    for (int lia = 0; lia < chunk->nda; ++lia) {
        // get the starting address for the chunk, taking into account the padding
        opt_double_ptr trg_data = chunk->data + chunk->size_padded * lia;
        opt_double_ptr src_data = data + localIndex(ax[0], listart[0], listart[1], listart[2], ax[0], nmem, nf, lia);

        // we alwas know that the chunk memory is aligned
        FLUPS_CHECK(m_isaligned(trg_data), "The chunk memory should be aligned, size_padded = %ld",chunk->size_padded);
        FLUPS_INFO("pointers are %p and %p",trg_data,src_data);
        FLUPS_INFO("copy %d %d %d from data to chunk",chunk->isize[0],chunk->isize[1],chunk->isize[2]);
        FLUPS_INFO("copy %zu bytes in %zu loops",nmax_byte,n_loop);
        FLUPS_INFO("local memory is %d %d %d, listart is %d %d %d", nmem[0], nmem[1], nmem[2], chunk->istart[0], chunk->istart[1], chunk->istart[2]);

#pragma omp for schedule(static)
        for (int il = 0; il < n_loop; ++il) {
            // get the local indexes (we cannot used the collaspedIndex one!!!)
            const int i2                  = il / (chunk->isize[ax[1]]);
            const int i1                  = il % (chunk->isize[ax[1]]);
            // get the starting adddress for the memcpy
            const double* __restrict vsrc = src_data + localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);
            double* __restrict vtrg       = trg_data + localIndex(ax0, 0, i1, i2, ax0, chunk->isize, nf, 0);
            FLUPS_INFO("offset src = %ld and trg = %ld",localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0),localIndex(ax0, 0, i1, i2, ax0, chunk->isize, nf, 0));
            std::memcpy(vtrg, vsrc, nmax_byte);
        }
    }
    //--------------------------------------------------------------------------
    END_FUNC;
}

