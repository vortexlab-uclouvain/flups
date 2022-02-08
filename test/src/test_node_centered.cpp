#include <mpi.h>
#include "base_convergence_test.hpp"


class NodeConvergenceTest : public BaseConvergenceTest {
    protected:
     void SetProblemDef_(int N) override {
         for (int dir = 0; dir < 3; dir++) {
             nglob_[dir] = N + 1;
             h_[dir] = L_[dir] / (nglob_[dir] - 1);
         }
     };

    public: 
    explicit NodeConvergenceTest(): BaseConvergenceTest(0.0, NODE_CENTER){};
    ~NodeConvergenceTest(){};
};


TEST_P(NodeConvergenceTest, AllBoundaryConditions){
    // ASSERT_LT(GetParam(), NMIXUNB);
    NodeConvergenceTest::PerformTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(NodeCentered,
                         NodeConvergenceTest,
                         testing::Range(0, 10));
// TEST_P(NodeConvergenceTest, Spectral){
//     ASSERT_LT(GetParam(), NSPECTRAL);
//     const int case_id = GetParam();
//     printf("The boundary conditions for this case are %d %d %d %d %d %d \n", spectral_bc[case_id][0],spectral_bc[case_id][1],spectral_bc[case_id][2],spectral_bc[case_id][3],spectral_bc[case_id][4],spectral_bc[case_id][5]);
    
//     // printf("I will perform %d different cases \n", n_case);
//     for(int i = 0; i < 8; i ++){
//         FLUPS_GreenType kernel = (FLUPS_GreenType) i;
//         printf("Testing kernel %s \n", NodeConvergenceTest::kname[(int)kernel].c_str());
//         NodeConvergenceTest::ActualTest(kernel, spectral_bc[case_id]);
//     }

// }

// TEST_P(NodeConvergenceTest, MixUnbounded){
//     ASSERT_LT(GetParam(), NMIXUNB);
//     const int case_id = GetParam();
//     printf("The boundary conditions for this case are %d %d %d %d %d %d \n", mix_unbounded_bc[case_id][0],mix_unbounded_bc[case_id][1],mix_unbounded_bc[case_id][2],mix_unbounded_bc[case_id][3],mix_unbounded_bc[case_id][4],mix_unbounded_bc[case_id][5]);
    
//     // printf("I will perform %d different cases \n", n_case);
//     for(int i = 3; i < 8; i ++){
//         FLUPS_GreenType kernel = (FLUPS_GreenType) i;
//         printf("Testing kernel %s \n", NodeConvergenceTest::kname[(int)kernel].c_str());
//         NodeConvergenceTest::ActualTest(kernel, mix_unbounded_bc[case_id]);
//     }

// }

// TEST_P(NodeConvergenceTest, PerPer){
//     ASSERT_LT(GetParam(), 8);
//     FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
//     NodeConvergenceTest::ActualTest(kernel, PER);
// }

// TEST_P(NodeConvergenceTest, OddOdd){
//     ASSERT_LT(GetParam(), 8);
//     FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
//     NodeConvergenceTest::ActualTest(kernel, ODD);
// }

