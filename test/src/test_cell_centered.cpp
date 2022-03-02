#include <mpi.h>

#include "base_convergence_test.hpp"


class CellConvergenceTest : public BaseConvergenceTest {
   protected:
    void SetProblemDef_(int N) override {
        for (int dir = 0; dir < 3; dir++) {
            nglob_[dir] = N;
            h_[dir]     = L_[dir] / (nglob_[dir]);
        }
    };

   public:
    explicit CellConvergenceTest() : BaseConvergenceTest(0.5, CELL_CENTER){};
    ~CellConvergenceTest(){};
};

TEST_P(CellConvergenceTest, AllBoundaryConditions){
    // ASSERT_LT(GetParam(), NMIXUNB);
    CellConvergenceTest::PerformTest(GetParam());
}

INSTANTIATE_TEST_SUITE_P(CellCentered,
                         CellConvergenceTest,
                         testing::Range(0, 1000));

// TEST_P(CellConvergenceTest, Spectral){
//     ASSERT_LT(GetParam(), NSPECTRAL);
//     const int case_id = GetParam();
//     printf("The boundary conditions for this case are %d %d %d %d %d %d \n", spectral_bc[case_id][0],spectral_bc[case_id][1],spectral_bc[case_id][2],spectral_bc[case_id][3],spectral_bc[case_id][4],spectral_bc[case_id][5]);
    
//     // printf("I will perform %d different cases \n", n_case);
//     for(int i = 0; i < 8; i ++){
//         FLUPS_GreenType kernel = (FLUPS_GreenType) i;
//         printf("Testing kernel %s \n", CellConvergenceTest::kname[(int)kernel].c_str());
//         CellConvergenceTest::ActualTest(kernel, spectral_bc[case_id]);
//     }

// }


// TEST_P(CellConvergenceTest, MixUnbounded){
//     ASSERT_LT(GetParam(), NMIXUNB);
//     const int case_id = GetParam();
//     printf("The boundary conditions for this case are %d %d %d %d %d %d \n", mix_unbounded_bc[case_id][0],mix_unbounded_bc[case_id][1],mix_unbounded_bc[case_id][2],mix_unbounded_bc[case_id][3],mix_unbounded_bc[case_id][4],mix_unbounded_bc[case_id][5]);
    
//     // printf("I will perform %d different cases \n", n_case);
//     for(int i = 3; i < 8; i ++){
//         FLUPS_GreenType kernel = (FLUPS_GreenType) i;
//         printf("Testing kernel %s \n", CellConvergenceTest::kname[(int)kernel].c_str());
//         CellConvergenceTest::ActualTest(kernel, mix_unbounded_bc[case_id]);
//     }

// }




    // CellConvergenceTest::InitBoundaryConditions(case_id);

    // CellConvergenceTest::PrintBcs();
    
    // // printf("I will perform %d different cases \n", n_case);
    // for(int i = 0; i < 7; i ++){
    //     FLUPS_GreenType kernel = (FLUPS_GreenType) i;
    //     printf("Testing kernel %s \n", CellConvergenceTest::kname[(int)kernel].c_str());
    //     CellConvergenceTest::ActualTest(kernel);
    // }
    // CellConvergenceTest::KillBoundaryConditions();


// TEST_P(CellConvergenceTest, FullyUnbounded){
//     const int n_case = 1;
//     ASSERT_LT(GetParam(), n_case);
//     const int case_id = GetParam();
//     printf("The boundary conditions for this case are %d %d %d %d %d %d \n", fully_unbounded_bc[case_id][0],fully_unbounded_bc[case_id][1],fully_unbounded_bc[case_id][2],fully_unbounded_bc[case_id][3],fully_unbounded_bc[case_id][4],fully_unbounded_bc[case_id][5]);
    
//     // printf("I will perform %d different cases \n", n_case);
//     for(int i = 0; i < 8; i ++){
//         FLUPS_GreenType kernel = (FLUPS_GreenType) i;
//         printf("Testing kernel %s \n", CellConvergenceTest::kname[(int)kernel].c_str());
//         CellConvergenceTest::ActualTest(kernel, fully_unbounded_bc[case_id]);
//     }

// }

// TEST_P(CellConvergenceTest, EvenEven){
//     ASSERT_LT(GetParam(), 8);
//     FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
//     CellConvergenceTest::ActualTest(kernel, EVEN);
// }

// TEST_P(CellConvergenceTest, PerPer){
//     ASSERT_LT(GetParam(), 8);
//     FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
//     CellConvergenceTest::ActualTest(kernel, PER);
// }

// TEST_P(CellConvergenceTest, OddOdd){
//     ASSERT_LT(GetParam(), 8);
//     FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
//     CellConvergenceTest::ActualTest(kernel, ODD);
// }

