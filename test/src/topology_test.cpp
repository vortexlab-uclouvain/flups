#include "gtest/gtest.h"

#include <mpi.h>

#include <string>
#include <vector>

#include "flups.h"

//==============================================================================
// TOPOLOGY DESCRIPTORS
//==============================================================================

/**
 * @brief one set of parameters used to build both a real Topology and its descriptor
 */
struct TopoParam {
    int  axis;
    int  lda;
    bool isComplex;
    int  nglob[3];
    int  nproc[3];
    int  axproc[3];
};

static const int topo_alignment_ = 16;

// the global sizes are chosen so that nglob%nproc != 0 in at least one direction
static const int topo_nglob_[3] = {17, 23, 9};
// with more than 2 ranks in one direction, some ranks own no pencil at all
static const int topo_nglob_small_[3] = {5, 3, 2};

/**
 * @brief returns a human readable description of a parameter set, used in the assert messages
 */
static std::string param_name(const TopoParam &p) {
    return "axis = " + std::to_string(p.axis) + ", lda = " + std::to_string(p.lda) +
           ", isComplex = " + std::to_string(p.isComplex) +
           ", nglob = " + std::to_string(p.nglob[0]) + " " + std::to_string(p.nglob[1]) + " " + std::to_string(p.nglob[2]) +
           ", nproc = " + std::to_string(p.nproc[0]) + " " + std::to_string(p.nproc[1]) + " " + std::to_string(p.nproc[2]) +
           ", axproc = " + std::to_string(p.axproc[0]) + " " + std::to_string(p.axproc[1]) + " " + std::to_string(p.axproc[2]);
}

/**
 * @brief build every parameter set to be tested for the given number of procs
 *
 * Both a size which is not a multiple of the number of procs and a size smaller than the number of
 * procs (some ranks then own nothing, which is allowed) are covered, as well as the scalar and the
 * vector cases.
 */
static std::vector<TopoParam> topo_params(const int comm_size) {
    const int nglob_list[2][3]  = {{topo_nglob_[0], topo_nglob_[1], topo_nglob_[2]},
                                   {topo_nglob_small_[0], topo_nglob_small_[1], topo_nglob_small_[2]}};
    const int nproc_list[3][3]  = {{comm_size, 1, 1}, {1, comm_size, 1}, {1, 1, comm_size}};
    const int axproc_list[2][3] = {{0, 1, 2}, {0, 2, 1}};

    std::vector<TopoParam> params;
    for (int ig = 0; ig < 2; ig++) {
        for (int ip = 0; ip < 3; ip++) {
            for (int ia = 0; ia < 2; ia++) {
                for (int axis = 0; axis < 3; axis++) {
                    for (int lda = 1; lda <= 3; lda += 2) {
                        for (int ic = 0; ic < 2; ic++) {
                            TopoParam p;
                            p.axis      = axis;
                            p.lda       = lda;
                            p.isComplex = (ic == 1);
                            for (int id = 0; id < 3; id++) {
                                p.nglob[id]  = nglob_list[ig][id];
                                p.nproc[id]  = nproc_list[ip][id];
                                p.axproc[id] = axproc_list[ia][id];
                            }
                            params.push_back(p);
                        }
                    }
                }
            }
        }
    }
    return params;
}

/**
 * @brief serialize every rank-explicit answer of a topology into a flat array of integers
 *
 * Only the questions which a descriptor is allowed to answer are asked, so the array of a real
 * topology and the array of its descriptor must be identical, whatever the rank building them.
 *
 * @param t the topology (real or descriptor)
 * @param comm_size the number of ranks of the decomposition
 */
static std::vector<int> topo_answers(const FLUPS_Topology *t, const int comm_size) {
    std::vector<int> ans;
    // the decomposition itself
    ans.push_back(flups_topo_get_axis(t));
    ans.push_back(flups_topo_get_lda(t));
    ans.push_back((int)flups_topo_get_isComplex(t));
    for (int id = 0; id < 3; id++) {
        ans.push_back(flups_topo_get_nglob(t, id));
        ans.push_back(flups_topo_get_nproc(t, id));
        ans.push_back(flups_topo_get_axproc(t, id));
    }
    // the rank owning every global id and the start id of every rank, in every direction
    for (int id = 0; id < 3; id++) {
        for (int ig = 0; ig < flups_topo_get_nglob(t, id); ig++) {
            ans.push_back(flups_topo_cmpt_rank_fromid(t, ig, id));
        }
        // the out of bound rank is included: it gives the end of the last rank
        for (int ir = 0; ir <= flups_topo_get_nproc(t, id); ir++) {
            ans.push_back(flups_topo_cmpt_start_id_from_rank(t, ir, id));
        }
    }
    // the rank decomposition, in both directions
    for (int ir = 0; ir < comm_size; ir++) {
        int rankd[3];
        flups_topo_ranksplit(t, ir, rankd);
        ans.push_back(rankd[0]);
        ans.push_back(rankd[1]);
        ans.push_back(rankd[2]);
        ans.push_back(flups_topo_rankindex(t, rankd));
    }
    return ans;
}

/**
 * @brief a descriptor must answer exactly what the real topology answers, on every rank
 */
TEST(Topology, DescriptorMatchesTopology) {
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const std::vector<TopoParam> params = topo_params(comm_size);
    ASSERT_FALSE(params.empty()) << "no valid parameter set for " << comm_size << " ranks";

    for (const TopoParam &p : params) {
        const std::string name = param_name(p);

        FLUPS_Topology *real = flups_topo_new(p.axis, p.lda, p.nglob, p.nproc, p.isComplex, p.axproc, topo_alignment_, MPI_COMM_WORLD);
        FLUPS_Topology *desc = flups_topo_new_desc(p.axis, p.lda, p.nglob, p.nproc, p.isComplex, p.axproc, topo_alignment_);

        ASSERT_FALSE(flups_topo_is_desc(real)) << name;
        ASSERT_TRUE(flups_topo_is_desc(desc)) << name;

        const std::vector<int> ans_real = topo_answers(real, comm_size);
        const std::vector<int> ans_desc = topo_answers(desc, comm_size);
        ASSERT_EQ(ans_real.size(), ans_desc.size()) << name;
        for (size_t i = 0; i < ans_real.size(); i++) {
            ASSERT_EQ(ans_real[i], ans_desc[i]) << name << " -- answer " << i;
        }

        // the local information of the real topology must match what the descriptor predicts for my rank
        int rank, rankd[3], istart[3];
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        flups_topo_ranksplit(desc, rank, rankd);
        flups_topo_get_istartGlob(real, istart);
        for (int id = 0; id < 3; id++) {
            const int start = flups_topo_cmpt_start_id_from_rank(desc, rankd[id], id);
            const int end   = flups_topo_cmpt_start_id_from_rank(desc, rankd[id] + 1, id);
            ASSERT_EQ(istart[id], start) << name << " -- istart[" << id << "]";
            ASSERT_EQ(flups_topo_get_nloc(real, id), end - start) << name << " -- nloc[" << id << "]";
        }

        flups_topo_free(real);
        flups_topo_free(desc);
    }
}

/**
 * @brief the ranks which are not members of a sub-communicator must answer the same as its members
 */
TEST(Topology, DescriptorOnSubComm) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size < 2) {
        GTEST_SKIP() << "at least 2 ranks are needed to build a strict sub-communicator";
    }

    //--------------------------------------------------------------------------
    /** - the topology lives on the first world_size-1 ranks only */
    //--------------------------------------------------------------------------
    const int  sub_size  = world_size - 1;
    const bool is_member = (world_rank < sub_size);
    MPI_Comm   sub_comm;
    MPI_Comm_split(MPI_COMM_WORLD, is_member ? 0 : MPI_UNDEFINED, world_rank, &sub_comm);
    ASSERT_EQ(is_member, sub_comm != MPI_COMM_NULL);

    const std::vector<TopoParam> params = topo_params(sub_size);
    ASSERT_FALSE(params.empty()) << "no valid parameter set for " << sub_size << " ranks";

    for (const TopoParam &p : params) {
        const std::string name = param_name(p);

        // every rank of the world can build the descriptor, member or not
        FLUPS_Topology  *desc     = flups_topo_new_desc(p.axis, p.lda, p.nglob, p.nproc, p.isComplex, p.axproc, topo_alignment_);
        std::vector<int> ans_desc = topo_answers(desc, sub_size);

        // the members compute the reference answers and share them with everybody (rank 0 is a member)
        std::vector<int> ans_ref(ans_desc.size(), 0);
        if (is_member) {
            FLUPS_Topology *real = flups_topo_new(p.axis, p.lda, p.nglob, p.nproc, p.isComplex, p.axproc, topo_alignment_, sub_comm);
            ans_ref              = topo_answers(real, sub_size);
            flups_topo_free(real);
        }
        ASSERT_EQ(ans_ref.size(), ans_desc.size()) << name;
        MPI_Bcast(ans_ref.data(), (int)ans_ref.size(), MPI_INT, 0, MPI_COMM_WORLD);

        for (size_t i = 0; i < ans_ref.size(); i++) {
            ASSERT_EQ(ans_ref[i], ans_desc[i]) << name << " -- answer " << i;
        }
        flups_topo_free(desc);
    }

    if (sub_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&sub_comm);
    }
}

//==============================================================================
// PENCIL LAYOUT
//==============================================================================

/**
 * @brief one set of boundary conditions, one per direction and per side
 */
struct BcSet {
    const char        *name;
    FLUPS_BoundaryType bc[3][2];
};

static const BcSet pencil_bcs_[] = {
    {"per", {{PER, PER}, {PER, PER}, {PER, PER}}},
    {"unb", {{UNB, UNB}, {UNB, UNB}, {UNB, UNB}}},
    {"even", {{EVEN, EVEN}, {EVEN, EVEN}, {EVEN, EVEN}}},
    {"mixed", {{PER, PER}, {UNB, UNB}, {EVEN, EVEN}}},
    {"mixunb", {{UNB, EVEN}, {PER, PER}, {EVEN, EVEN}}},
};

static const double pencil_L_[3] = {1.0, 1.0, 1.0};

/**
 * @brief allocate the FLUPS bc arrays (one value per component) for a given bc set
 */
static void alloc_bc(const BcSet &set, const int lda, FLUPS_BoundaryType *bc[3][2]) {
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            bc[id][is] = (FLUPS_BoundaryType *)flups_malloc(lda * sizeof(int));
            for (int lia = 0; lia < lda; lia++) {
                bc[id][is][lia] = set.bc[id][is];
            }
        }
    }
}

static void free_bc(FLUPS_BoundaryType *bc[3][2]) {
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(bc[id][is]);
        }
    }
}

/**
 * @brief the description of one case, used in the assert messages
 */
static std::string pencil_name(const BcSet &set, const FLUPS_CenterType center, const int lda, const int nproc[3]) {
    return std::string(set.name) + (center == NODE_CENTER ? "_node" : "_cell") + "_lda" + std::to_string(lda) +
           "_nproc" + std::to_string(nproc[0]) + std::to_string(nproc[1]) + std::to_string(nproc[2]);
}

/**
 * @brief the layout predicted without any communication must be the one the solver actually builds
 */
TEST(PencilLayout, MatchesSolverTopology) {
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nproc_list[2][3] = {{1, 1, comm_size}, {1, comm_size, 1}};

    for (const BcSet &set : pencil_bcs_) {
        for (int ictr = 0; ictr < 2; ictr++) {
            const FLUPS_CenterType center         = (ictr == 0) ? NODE_CENTER : CELL_CENTER;
            const FLUPS_CenterType center_type[3] = {center, center, center};
            const int              nglob[3]       = {(center == NODE_CENTER) ? 49 : 48,
                                                     (center == NODE_CENTER) ? 33 : 32,
                                                     (center == NODE_CENTER) ? 41 : 40};
            double                 h[3];
            for (int id = 0; id < 3; id++) {
                h[id] = pencil_L_[id] / (nglob[id] - (center == NODE_CENTER));
            }

            for (int lda = 1; lda <= 3; lda += 2) {
                FLUPS_BoundaryType *bc[3][2];
                alloc_bc(set, lda, bc);

                for (int ip = 0; ip < 2; ip++) {
                    const int         nproc[3] = {nproc_list[ip][0], nproc_list[ip][1], nproc_list[ip][2]};
                    const std::string name     = pencil_name(set, center, lda, nproc);

                    //----------------------------------------------------------
                    /** - the pure prediction, which needs no solver at all */
                    //----------------------------------------------------------
                    int axis_pred, nglob_pred[3], nproc_pred[3], axproc_pred[3];
                    flups_pencil_layout(nglob, nproc, bc, h, pencil_L_, center_type, lda, &axis_pred, nglob_pred, nproc_pred, axproc_pred);

                    //----------------------------------------------------------
                    /** - the real thing */
                    //----------------------------------------------------------
                    FLUPS_Topology *topo   = flups_topo_new(0, lda, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, MPI_COMM_WORLD);
                    FLUPS_Solver   *solver = flups_init(topo, bc, h, pencil_L_, NOD, center_type);
                    flups_set_greenType(solver, CHAT_2);
                    flups_setup(solver, false);

                    const FLUPS_Topology *inner = flups_get_innerTopo_physical(solver);
                    ASSERT_FALSE(flups_topo_get_isComplex(inner)) << name << " -- the inner physical topo must be real";
                    ASSERT_EQ(flups_topo_get_axis(inner), axis_pred) << name << " -- axis";
                    for (int id = 0; id < 3; id++) {
                        ASSERT_EQ(flups_topo_get_nglob(inner, id), nglob_pred[id]) << name << " -- nglob[" << id << "]";
                        ASSERT_EQ(flups_topo_get_nproc(inner, id), nproc_pred[id]) << name << " -- nproc[" << id << "]";
                        ASSERT_EQ(flups_topo_get_axproc(inner, id), axproc_pred[id]) << name << " -- axproc[" << id << "]";
                    }

                    flups_cleanup_solver(solver);
                    flups_topo_free(topo);
                }
                free_bc(bc);
            }
        }
    }
    flups_cleanup_backend();
}

/**
 * @brief a rank outside of the solver communicator must describe the pencil topology as its members do
 */
TEST(PencilLayout, DescriptorOnSubComm) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size < 2) {
        GTEST_SKIP() << "at least 2 ranks are needed to build a strict sub-communicator";
    }

    //--------------------------------------------------------------------------
    /** - the solver lives on the first world_size-1 ranks only */
    //--------------------------------------------------------------------------
    const int  sub_size  = world_size - 1;
    const bool is_member = (world_rank < sub_size);
    MPI_Comm   sub_comm;
    MPI_Comm_split(MPI_COMM_WORLD, is_member ? 0 : MPI_UNDEFINED, world_rank, &sub_comm);
    ASSERT_EQ(is_member, sub_comm != MPI_COMM_NULL);

    const int nproc[3] = {1, 1, sub_size};

    for (const BcSet &set : pencil_bcs_) {
        for (int ictr = 0; ictr < 2; ictr++) {
            const FLUPS_CenterType center         = (ictr == 0) ? NODE_CENTER : CELL_CENTER;
            const FLUPS_CenterType center_type[3] = {center, center, center};
            const int              nglob[3]       = {(center == NODE_CENTER) ? 49 : 48,
                                                     (center == NODE_CENTER) ? 33 : 32,
                                                     (center == NODE_CENTER) ? 41 : 40};
            double                 h[3];
            for (int id = 0; id < 3; id++) {
                h[id] = pencil_L_[id] / (nglob[id] - (center == NODE_CENTER));
            }

            for (int lda = 1; lda <= 3; lda += 2) {
                const std::string   name = pencil_name(set, center, lda, nproc);
                FLUPS_BoundaryType *bc[3][2];
                alloc_bc(set, lda, bc);

                //--------------------------------------------------------------
                /** - every rank of the world predicts the layout and builds the descriptor */
                //--------------------------------------------------------------
                int axis_pred, nglob_pred[3], nproc_pred[3], axproc_pred[3];
                flups_pencil_layout(nglob, nproc, bc, h, pencil_L_, center_type, lda, &axis_pred, nglob_pred, nproc_pred, axproc_pred);
                FLUPS_Topology *desc = flups_topo_new_desc(axis_pred, lda, nglob_pred, nproc_pred, false, axproc_pred, FLUPS_ALIGNMENT);

                std::vector<int> ans_desc = topo_answers(desc, sub_size);
                std::vector<int> ans_ref(ans_desc.size(), 0);

                //--------------------------------------------------------------
                /** - only the members build the solver, and share its answers with everybody */
                //--------------------------------------------------------------
                if (is_member) {
                    FLUPS_Topology *topo   = flups_topo_new(0, lda, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, sub_comm);
                    FLUPS_Solver   *solver = flups_init(topo, bc, h, pencil_L_, NOD, center_type);
                    flups_set_greenType(solver, CHAT_2);
                    flups_setup(solver, false);

                    // no assert before the Bcast below: it is collective on the world
                    ans_ref = topo_answers(flups_get_innerTopo_physical(solver), sub_size);

                    flups_cleanup_solver(solver);
                    flups_topo_free(topo);
                }
                // rank 0 is always a member of the sub-communicator
                MPI_Bcast(ans_ref.data(), (int)ans_ref.size(), MPI_INT, 0, MPI_COMM_WORLD);

                //--------------------------------------------------------------
                /** - members and non-members must answer exactly the same */
                //--------------------------------------------------------------
                ASSERT_EQ(ans_ref[0], axis_pred) << name << " -- axis";
                ASSERT_EQ(ans_ref[1], lda) << name << " -- lda";
                ASSERT_EQ(ans_ref[2], 0) << name << " -- the inner physical topo must be real";
                for (int id = 0; id < 3; id++) {
                    ASSERT_EQ(ans_ref[3 + 3 * id + 0], nglob_pred[id]) << name << " -- nglob[" << id << "]";
                    ASSERT_EQ(ans_ref[3 + 3 * id + 1], nproc_pred[id]) << name << " -- nproc[" << id << "]";
                    ASSERT_EQ(ans_ref[3 + 3 * id + 2], axproc_pred[id]) << name << " -- axproc[" << id << "]";
                }
                for (size_t i = 0; i < ans_ref.size(); i++) {
                    ASSERT_EQ(ans_ref[i], ans_desc[i]) << name << " -- answer " << i;
                }

                flups_topo_free(desc);
                free_bc(bc);
            }
        }
    }

    if (sub_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&sub_comm);
    }
    // only a member initialised the backend in this test (flups_init), so only a member may
    // release it: backend_finalize destroys the global handle unconditionally
    if (is_member) {
        flups_cleanup_backend();
    }
}
