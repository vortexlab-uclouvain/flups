#include <mpi.h>
#include "gtest/gtest.h"
#include "flups.h"

#include "tools_validation.hpp"

#define DOUBLE_TOL 1e-13
#define CONV_TOL   0.05

double Error_Inf(int N, FLUPS_GreenType kernel, FLUPS_BoundaryType bc_in){
    //-------------------------------------------------------------------------
    //Definition of the problem
    //-------------------------------------------------------------------------
    const int     nglob[3] = {N+1, N+1, N+1};
    const int     nproc[3] = {1, 1, 1};
    const double  L[3]     = {1., 1., 1.};

    const double h[3] = {L[0] / (nglob[0]-1), L[1] / (nglob[1]-1), L[2] / (nglob[2]-1)};

    const FLUPS_CenterType center_type[3] = {NODE_CENTER, NODE_CENTER, NODE_CENTER};

    FLUPS_BoundaryType *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            mybc[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
            mybc[id][is][0] = bc_in;
        }
    }

    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------
    // create a real topology
    FLUPS_Topology *topo      = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, MPI_COMM_WORLD);

    // solver creation and init
    FLUPS_Solver *mysolver = flups_init(topo, mybc, h, L, NOD, center_type);
    flups_set_greenType(mysolver, kernel);
    flups_setup(mysolver, true);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);

    const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;

    const double cos_cons = 2.0 * M_PI * 3;
    
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2) ; i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1) ; i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                double       x     = (istart[ax0] + i0) * h[ax0];
                double       y     = (istart[ax1] + i1) * h[ax1];
                double       z     = (istart[ax2] + i2) * h[ax2];
                const size_t id    = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                // second derivative of a cosine
                rhs[id] = -(cos_cons * cos_cons) * cos(cos_cons * x);
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - solve the equations */
    //-------------------------------------------------------------------------
    for (int is = 0; is < 1; is++) {
        flups_solve(mysolver, field, rhs, STD);
    }

    double err = std::numeric_limits<double>::lowest();
    for (int i2 = 1; i2 < flups_topo_get_nloc(topo,ax2) -1   ; i2++) {
        for (int i1 = 1; i1 < flups_topo_get_nloc(topo,ax1) -1; i1++) {
            for (int i0 = 1; i0 < flups_topo_get_nloc(topo,ax0)-1; i0++) {
                double       x     = (istart[ax0] + i0) * h[ax0];
                const size_t id    = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                err = std::max(err, abs(field[id] - cos(cos_cons * x)));
            }
        }
    }

    // //-------------------------------------------------------------------------
    // /** - CLEAN */
    // //-------------------------------------------------------------------------

    flups_free(rhs);
    flups_free(field);

    flups_topo_free(topo);
    flups_cleanup(mysolver);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }
    return err;
}

class NodeCenteredTest : public testing::TestWithParam<int> {
   public: 
        void Actual_test(FLUPS_GreenType kernel, FLUPS_BoundaryType boundary){
                int Nmin = 128;
                if (kernel == CHAT_2 || kernel == HEJ_0) {
                    double err = Error_Inf(Nmin, kernel, boundary);
                    printf(" The convergence of %s is spectral. You obtained a error of %e \n \n", kname[(int) kernel].c_str(), err);
                    ASSERT_NEAR(err, 0.0, DOUBLE_TOL);
                } else {
                    double erri[2];
                    erri[0]               = Error_Inf(Nmin, kernel, boundary);
                    erri[1]               = Error_Inf(2 * Nmin, kernel, boundary);
                    double expected_order = KernelOrder(kernel);
                    double computed_order = -log(erri[1] / erri[0]) / log(2.0);
                    printf(" The choosen kernel is %s has an order of %2.0f. You obtained a convergence of %2.3f \n \n", kname[(int)kernel].c_str(), expected_order, computed_order);
                    ASSERT_GE(computed_order, expected_order - CONV_TOL);
                }
        }

   protected:
    std::string kname[8]  = {"CHAT_2", "LGF_2", "HEJ_2", "HEJ_4", "HEJ_6", "HEJ_8", "HEJ_10", "HEJ_0"};
    void SetUp() override{};
    void TearDown() override{};
};

TEST_P(NodeCenteredTest, EvenEven){
    ASSERT_LT(GetParam(), 8);
    FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
    NodeCenteredTest::Actual_test(kernel, EVEN);
}

TEST_P(NodeCenteredTest, PerPer){
    ASSERT_LT(GetParam(), 8);
    FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
    NodeCenteredTest::Actual_test(kernel, PER);
}

INSTANTIATE_TEST_SUITE_P(NodeCentered,
                         NodeCenteredTest,
                         testing::Range(0, 8));