#include <mpi.h>


#include "tools_validation.hpp"


class NodeConvergenceTest : public BaseConvergenceTest {
    protected:
     void InitFlups_(int N, FLUPS_GreenType kernel, const FLUPS_BoundaryType bc_in[6]) override{
        //  -------------------------------------------------------------------------
        //  Definition of the problem
        //  -------------------------------------------------------------------------
         const int     nglob[3] = {N+1, N+1, N+1};
         const int     nproc[3] = {1, 1, 1};
         
         for(int dir = 0 ; dir < 3; dir ++ ){
             h_[dir] = L_[dir] / (nglob[dir]-1);
         }

         const FLUPS_CenterType center_type[3] = {NODE_CENTER, NODE_CENTER, NODE_CENTER};

         for (int id = 0; id < 3; id++) {
             for (int is = 0; is < 2; is++) {
                 mybc_[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
                 mybc_[id][is][0] = bc_in[is + 2*id];
             }
         }

         //-------------------------------------------------------------------------
         /** - Initialize FLUPS */
         //-------------------------------------------------------------------------
         // create a real topology
         topo_ = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, MPI_COMM_WORLD);

         // solver creation and init
         mysolver_ = flups_init(topo_, mybc_, h_, L_, NOD, center_type);
         flups_set_greenType(mysolver_, kernel);
         flups_setup(mysolver_, true);
     };

    public: 
    explicit NodeConvergenceTest(): BaseConvergenceTest(0.0){};
    ~NodeConvergenceTest(){};
};


TEST_P(NodeConvergenceTest, Spectral){
    ASSERT_LT(GetParam(), NSPECTRAL);
    const int case_id = GetParam();
    printf("The boundary conditions for this case are %d %d %d %d %d %d \n", spectral_bc[case_id][0],spectral_bc[case_id][1],spectral_bc[case_id][2],spectral_bc[case_id][3],spectral_bc[case_id][4],spectral_bc[case_id][5]);
    
    // printf("I will perform %d different cases \n", n_case);
    for(int i = 0; i < 8; i ++){
        FLUPS_GreenType kernel = (FLUPS_GreenType) i;
        printf("Testing kernel %s \n", NodeConvergenceTest::kname[(int)kernel].c_str());
        NodeConvergenceTest::ActualTest(kernel, spectral_bc[case_id]);
    }

}

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

INSTANTIATE_TEST_SUITE_P(NodeCentered,
                         NodeConvergenceTest,
                         testing::Range(0, NSPECTRAL));