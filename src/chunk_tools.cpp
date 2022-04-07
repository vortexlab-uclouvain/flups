

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
    int srank[3], erank[3] for (int id = 0; id < 3; ++id) {
        // find the rank that will own the starting/ending point of the input topology
        // the ending point is taken with -1 to be sure to include the end rank
        srank[id] = topo_out->cmpt_rank_fromid(topoi_start[id] + shift[id]);
        erank[id] = topo_out->cmpt_rank_fromid(topoi_end[id] - 1 + shift[id]);

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
 * @param nf the maximum nf between
 * @param chunk
 */
void PlanShuffleChunk(const bool iscomplex, MemChunk* chunk, double* data) {
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
    FLUPS_INFO_3("shuffle: setting up the shuffle form %d to %d", topo_in->axis(), topo_out->axis());
    FLUPS_INFO_3("shuffle: nf = %d, blocksize = %d %d %d", nf, bSize[0], bSize[1], bSize[2]);
    FLUPS_INFO_3("shuffle: DIM 0: n = %d, is=%d, os=%d", dims[0].n, dims[0].is, dims[0].os);
    FLUPS_INFO_3("shuffle: DIM 1: n = %d, is=%d, os=%d", dims[1].n, dims[1].is, dims[1].os);

    // plan the real or complex plan
    // the nf is driven by the OUT topology ALWAYS
    if (!iscomplex) {
        *shuffle = fftw_plan_guru_r2r(0, NULL, 2, dims, data, data, NULL, FFTW_FLAG);
        FLUPS_CHECK(*shuffle != NULL, "Plan has not been setup");
    } else {
        *shuffle = fftw_plan_guru_dft(0, NULL, 2, dims, (fftw_complex*)data, (fftw_complex*)data, FLUPS_FORWARD, FFTW_FLAG);
        FLUPS_CHECK(*shuffle != NULL, "Plan has not been setup");
    }
}