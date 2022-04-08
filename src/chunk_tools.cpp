#include "chunk_tools.hpp"

void PopulateChunk(const int shift[3], const Topology* topo_in, const Topology* topo_out, int* n_chunks, MemChunk* chunks) {
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

    //--------------------------------------------------------------------------
    /** - get the number of chunks to build */
    //--------------------------------------------------------------------------
    // get the real start and end index that will exist in both topologies
    int topoi_start[3], topoi_end[3];
    topo_in->cmpt_intersect_id(shift, topo_out, topoi_start, topoi_end);

    n_chunks[0] = 1;
    int srank[3], erank[3];
    for (int id = 0; id < 3; ++id) {
        // find the rank that will own the starting/ending point of the input topology
        // the ending point is taken with -1 to be sure to include the end rank
        srank[id] = topo_out->cmpt_rank_fromid(topoi_start[id] + shift[id],id);
        erank[id] = topo_out->cmpt_rank_fromid(topoi_end[id] - 1 + shift[id],id);

        // accumulate the number of chunks over the dimensions
        // both ranks are included!
        n_chunks[0] *= (erank[id] - srank[id] + 1);
    }
    chunks = reinterpret_cast<MemChunk*>(m_calloc(n_chunks[0] * sizeof(MemChunk)));

    //--------------------------------------------------------------------------
    /** - fill the chunks */
    //--------------------------------------------------------------------------
    int chunk_counter = 0;
    for (int ir2 = srank[2]; ir2 <= erank[2]; ++ir2) {
        for (int ir1 = srank[1]; ir1 <= erank[1]; ++ir1) {
            for (int ir0 = srank[0]; ir0 <= erank[0]; ++ir0) {
                // get the current chunk
                MemChunk* cchunk = chunks + chunk_counter;
                // store the current rank in a XYZ format
                const int irank[3] = {ir0, ir1, ir2};
                for (int id = 0; id < 3; ++id) {
                    // get the start and end index in topo OUT
                    const int topoo_start = topo_out->cmpt_start_id_from_rank(irank[id], id);
                    const int topoo_end   = topo_out->cmpt_start_id_from_rank(irank[id] + 1, id);

                    // make sure that the chunks belongs to the topo_in
                    cchunk->istart[id] = m_max(topoo_start - shift[id], topoi_start[id]);
                    cchunk->isize[id]  = m_min(topoo_end - shift[id], topoi_end[id]) - cchunk->istart[id];

                    // get the destination rank, no translation is required here as we imposed identical communicators
                    cchunk->dest_rank = rankindex(irank, topo_out);

                    cchunk->nda         = topo_in->lda();
                    cchunk->size_padded = get_ChunkPaddedSize(topo_in->nf(), cchunk);

                    // fill the axis
                    cchunk->iaxis = topo_in->axis();
                    cchunk->oaxis = topo_out->axis();
                }

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

    int iaxis[3] = {chunk->iaxis, (chunk->iaxis + 1) % 3, (chunk->iaxis + 2) % 3};
    int oaxis[3] = {chunk->oaxis, (chunk->oaxis + 1) % 3, (chunk->oaxis + 2) % 3};

    // note: this is a bit magical and I forgot the exact why we do that
    // loop over the input axis and update the dim[0] and dim[1] data structures accordingly
    for (int id = 0; id < 3; id++) {
        // if the axis we are currently viewing is not the output axis, update the size and the input stride
        if (iaxis[id] != chunk->oaxis) {
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
        if (oaxis[id] != chunk->iaxis) {
            dims[0].n  = dims[0].n * chunk->isize[oaxis[id]];
            dims[1].os = dims[1].os * chunk->isize[oaxis[id]];
        } else {
            break;
        }
    }
    // display some info
    FLUPS_INFO("shuffle: setting up the shuffle form %d to %d", chunk->iaxis, chunk->oaxis);
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
}