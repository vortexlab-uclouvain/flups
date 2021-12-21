#include <mpi.h>


#include "tools_validation.hpp"


class NodeConvergenceTest : public BaseConvergenceTest {
    protected:
     void InitFlups_(int N, FLUPS_GreenType kernel, FLUPS_BoundaryType bc_in) override{
        //  -------------------------------------------------------------------------
        //  Definition of the problem
        //  -------------------------------------------------------------------------
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
     };

    public: 
    explicit NodeConvergenceTest(): BaseConvergenceTest(0.0){};
    ~NodeConvergenceTest(){};
};

// class CellConvergenceTest : public BaseConvergenceTest{

// };

TEST_P(NodeConvergenceTest, EvenEven){
    ASSERT_LT(GetParam(), 8);
    FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
    NodeConvergenceTest::ActualTest(kernel, EVEN);
}

TEST_P(NodeConvergenceTest, PerPer){
    ASSERT_LT(GetParam(), 8);
    FLUPS_GreenType kernel = (FLUPS_GreenType) GetParam();
    NodeConvergenceTest::ActualTest(kernel, PER);
}

INSTANTIATE_TEST_SUITE_P(NodeCentered,
                         NodeConvergenceTest,
                         testing::Range(0, 8));