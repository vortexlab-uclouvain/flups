#include "gtest/gtest.h"

#include <limits>
#include <mpi.h>
#include <map>

#include "analytical_field.hpp"

#define ZERO_TOL 1000.0 * std::numeric_limits<double>::epsilon() 
#define N_TEST 48

// Specific log function for the test library 
#define test_log(format, ...)                              \
    ({                                                     \
        int test_log_rank_;                                \
        MPI_Comm_rank(MPI_COMM_WORLD, &test_log_rank_);    \
        char test_log_msg_[1024];                          \
        sprintf(test_log_msg_, format, ##__VA_ARGS__);     \
        if (test_log_rank_ == 0)                           \
        {                                                  \
            fprintf(stdout, "[test] %s \n", test_log_msg_); \
        }                                                  \
    })

// List of all the compatible boundary conditions for 1 direction
static const std::array<FLUPS_BoundaryType,2> EV_EV = {EVEN, EVEN};
static const std::array<FLUPS_BoundaryType,2> OD_OD = {ODD, ODD};
static const std::array<FLUPS_BoundaryType,2> UN_UN = {UNB, UNB};
static const std::array<FLUPS_BoundaryType,2> PE_PE = {PER, PER};
static const std::array<FLUPS_BoundaryType,2> EV_OD = {EVEN, ODD};
static const std::array<FLUPS_BoundaryType,2> OD_EV = {ODD, EVEN};
static const std::array<FLUPS_BoundaryType,2> EV_UN = {EVEN, UNB};
static const std::array<FLUPS_BoundaryType,2> UN_EV = {UNB, EVEN};
static const std::array<FLUPS_BoundaryType,2> OD_UN = {ODD, UNB};
static const std::array<FLUPS_BoundaryType,2> UN_OD = {UNB, ODD};
static const std::array<FLUPS_BoundaryType,2> NO_NO = {NONE, NONE};

// Output for log
// static const std::map<FLUPS_CenterType, std::string> cname = {{NODE_CENTER, "NODE_CENTER"}, {CELL_CENTER, "CELL_CENTER"}};
// static const std::map<FLUPS_GreenType, std::string> kname = {{CHAT_2, "CHAT_2"}, {LGF_2, "LGF_2"}, {HEJ_2, "HEJ_2"}, {HEJ_4, "HEJ_4"}, {HEJ_6, "HEJ_6"}, {HEJ_8, "HEJ_8"}, {HEJ_10, "HEJ_10"}, {HEJ_0, "HEJ_0"}};
// static const std::map<FLUPS_BoundaryType, std::string> bname = {{PER, "PER"}, {UNB, "UNB"}, {ODD, "ODD"}, {EVEN, "EVEN"}}; 
static const std::map<FLUPS_CenterType, std::string> cname = {{NODE_CENTER, "node"}, {CELL_CENTER, "cell"}};
static const std::map<FLUPS_GreenType, std::string> kname = {{CHAT_2, "chat2"}, {LGF_2, "lgf2"}, {LGF_4, "lgf4"}, {LGF_6, "lgf6"}, {LGF_8, "lgf8"}, {HEJ_2, "hej2"}, {HEJ_4, "hej4"}, {HEJ_6, "hej6"}, {HEJ_8, "hej8"}, {HEJ_10, "hej10"}, {HEJ_0, "hej0"}};
static const std::map<FLUPS_BoundaryType, std::string> bname = {{EVEN, "0"}, {ODD, "1"}, {PER, "3"}, {UNB, "4"}}; 


/**
 * @brief The base class for a convergence test. 
 * Both the Node Centered and the Cell Centered test uses the exact same functions to perform the tests. 
 * However, the function to compute the number of points must be overwritten
 * 
 */

using ParamType = std::tuple<FLUPS_CenterType, FLUPS_GreenType,
                                std::array<FLUPS_BoundaryType, 2>, 
                                std::array<FLUPS_BoundaryType, 2>, 
                                std::array<FLUPS_BoundaryType ,2> >;

std::string TestNameGenerator(const ::testing::TestParamInfo<ParamType>& info) {
    const ParamType param = info.param;
    FLUPS_CenterType center = std::get<0>(param);
    FLUPS_GreenType green = std::get<1>(param);
    FLUPS_BoundaryType bdy[6] = {std::get<2>(param)[0], std::get<2>(param)[1], std::get<3>(param)[0], std::get<3>(param)[1], std::get<4>(param)[0], std::get<4>(param)[1]};

    std::string test_name = "";
    test_name += cname.at(center);
    test_name += "_" + kname.at(green);
    test_name += "_" + bname.at(bdy[0]) + bname.at(bdy[1]) + "_" + bname.at(bdy[2]) + bname.at(bdy[3]) + "_" + bname.at(bdy[4]) + bname.at(bdy[5]);
    return test_name;
}

class ConvergenceTest : public testing::TestWithParam<ParamType>{
protected:
    // Default variable for all the tests
    // const int    nproc_[3]  = {4, 4, 4};
    const int    nproc_[3]  = {1, 1, 1};
    const double L_[3]      = {1., 1., 1.};

    // Variable specific to a test_case
    FLUPS_BoundaryType *mybc_[3][2];
    FLUPS_GreenType green_;
    FLUPS_CenterType center_type_[3];
    double           shift_;
    
    void SetUp() override {
        //--------------------------------------------------------------------
        // Node or cell center
        for(int idim = 0; idim < 3; ++idim){
            center_type_[idim] = std::get<0>(GetParam());
        }        
        shift_ = (std::get<0>(GetParam()) == CELL_CENTER) ? 0.5 : 0.0;

        //..........................................
        // Kernel 
        green_ = std::get<1>(GetParam());
        
        //..........................................
        // Boundary conditions
        for(int idim = 0; idim < 3; ++idim ){
            for(int ibc = 0; ibc < 2; ++ibc){
                mybc_[idim][ibc] = (FLUPS_BoundaryType*)(flups_malloc(1*sizeof(int)));
            }
        }
        // x boundary condition
        *(mybc_[0][0]) = std::get<2>(GetParam())[0];
        *(mybc_[0][1]) = std::get<2>(GetParam())[1];
        
        // y boundary condition
        *(mybc_[1][0]) = std::get<3>(GetParam())[0];
        *(mybc_[1][1]) = std::get<3>(GetParam())[1];
        
        // z boundary condition
        *(mybc_[2][0]) = std::get<4>(GetParam())[0];
        *(mybc_[2][1]) = std::get<4>(GetParam())[1];

        //..........................................
        test_log("===========================================================");
        test_log("TESTING: data type = %s; kernel = %s and bc = %d %d - %d %d - %d %d", 
        (center_type_[0] == NODE_CENTER)? "node-centered" : "cell-centered", kname.at(green_).c_str(),
        *mybc_[0][0], *mybc_[0][1], *mybc_[1][0], *mybc_[1][1], *mybc_[2][0], *mybc_[2][1]);
        test_log("===========================================================");
        //--------------------------------------------------------------------
    };

    void TearDown() override {
        //--------------------------------------------------------------------
        for(int idim = 0; idim < 3; ++idim ){
            for(int ibc = 0; ibc < 2; ++ibc){
                flups_free(mybc_[idim][ibc]);
            }
        }
        //--------------------------------------------------------------------
    };
};

TEST_P(ConvergenceTest, AllBoundaryConditions){
    //--------------------------------------------------------------------
    // Check if we deal with a spectral case or an mix cases
    bool isXSpectral = (*(mybc_[0][0]) != UNB) && (*(mybc_[0][1]) != UNB); 
    bool isYSpectral = (*(mybc_[1][0]) != UNB) && (*(mybc_[1][1]) != UNB);
    bool isZSpectral = (*(mybc_[2][0]) != UNB) && (*(mybc_[2][1]) != UNB);
    int  n_spectral = isXSpectral + isYSpectral + isZSpectral;
    bool is_green_spectral = (n_spectral == 3) && (CHAT_2 == green_ || HEJ_0 == green_); // Chat2 and Hej0 have a spectral accuracy with spectral bc
    bool is_green_invalid = (HEJ_0 == green_) && (n_spectral != 3); // Hej0 cannot handle unbounded bcs.
    is_green_invalid |= (LGF_2 == green_) && (n_spectral == 1); // LGF_2 cannot handle one spectral direction
    is_green_invalid |= (LGF_4 == green_) && (n_spectral == 1); // LGF_4 cannot handle one or two spectral directions
    is_green_invalid |= (LGF_6 == green_) && (n_spectral == 1); // LGF_6 cannot handle one or two spectral directions
    is_green_invalid |= (LGF_8 == green_) && (n_spectral == 1); // LGF_8 cannot handle one or two spectral directions
    
    //  -------------------------------------------------------------------------
    // Perform the test
    double erri[2]; 
    int    ntest  = (is_green_invalid) ? 0 : (is_green_spectral ? 1 : 2) ;
    
    for(int i = 0; i < ntest; ++i){
        //  -------------------------------------------------------------------------
        //  Definition of the problem
        int nglob[3]; double h[3];
        for (int dir = 0; dir < 3; dir++) {
            nglob[dir] = (i+1) * N_TEST + (NODE_CENTER == center_type_[0]);  // Need to have an additional data if we have node-centered data
            h[dir] = L_[dir] / (nglob[dir] - (NODE_CENTER == center_type_[0])); 
         }
        
        //-------------------------------------------------------------------------
        // Initialize FLUPS
        test_log("Setting up solver");
        FLUPS_Topology* mytopo = flups_topo_new(0, 1, nglob, nproc_, false, NULL, FLUPS_ALIGNMENT, MPI_COMM_WORLD);
        FLUPS_Solver* mysolver = flups_init(mytopo, mybc_, h, L_, NOD, center_type_);
        flups_set_greenType(mysolver, green_);
        flups_setup(mysolver, true);
        

        //-------------------------------------------------------------------------
        // allocate rhs and solution 
        double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(mytopo));
        double *sol   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(mytopo));
        double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(mytopo));

        int istart[3];
        flups_topo_get_istartGlob(mytopo, istart);
        const int nmem[3] = {flups_topo_get_nmem(mytopo, 0), flups_topo_get_nmem(mytopo, 1), flups_topo_get_nmem(mytopo, 2)};
        const int ax0 = flups_topo_get_axis(mytopo);
        const int ax1 = (ax0 + 1) % 3;
        const int ax2 = (ax0 + 2) % 3;

        //-------------------------------------------------------------------------
        // find the correct analytical functions 
        test_log("Allocate analytical Field ");
        AnalyticalField analytics[3];
        for (int dir = 0; dir < 3; ++dir){
            analytics[dir].SetParam(*mybc_[dir][0], *mybc_[dir][1]);
        }
        analytics[0].SetFreq(1.);
        analytics[1].SetFreq(2.);
        analytics[2].SetFreq(4.);

        for (int i2 = 0; i2 < flups_topo_get_nloc(mytopo, ax2); i2++){
            for (int i1 = 0; i1 < flups_topo_get_nloc(mytopo, ax1); i1++){
                for (int i0 = 0; i0 < flups_topo_get_nloc(mytopo, ax0); i0++){
                    const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                    sol[id] = 1.0;
                }
            }
        }

        for (int i2 = 0; i2 < flups_topo_get_nloc(mytopo, ax2); i2++){
            for (int i1 = 0; i1 < flups_topo_get_nloc(mytopo, ax1); i1++){
                for (int i0 = 0; i0 < flups_topo_get_nloc(mytopo, ax0); i0++){
                    const size_t id = flups_locID(0, i0, i1, i2, 0, ax0, nmem, 1);
                    const double x[3] = {(istart[ax0] + i0 + shift_) * h[ax0],
                                         (istart[ax1] + i1 + shift_) * h[ax1],
                                         (istart[ax2] + i2 + shift_) * h[ax2]};
                    for (int dir0 = 0; dir0 < 3; dir0++){
                        const int dir1 = (dir0 + 1) % 3;
                        const int dir2 = (dir0 + 2) % 3;
                        sol[id] *= analytics[dir0].Sol(x[dir0], L_[dir0]);
                        rhs[id] += analytics[dir0].Rhs(x[dir0], L_[dir0]) * analytics[dir1].Sol(x[dir1], L_[dir1]) * analytics[dir2].Sol(x[dir2], L_[dir2]);
                    }
                }
            }
        }

        //-------------------------------------------------------------------------
        // solve the equations
        test_log("Solve the equation");
        for (int is = 0; is < 1; is++){
            flups_solve(mysolver, field, rhs, STD);
        }

        double local_err = std::numeric_limits<double>::lowest();
        double global_err = std::numeric_limits<double>::lowest();
        for (int i2 = 0; i2 < flups_topo_get_nloc(mytopo, ax2); i2++){
            for (int i1 = 0; i1 < flups_topo_get_nloc(mytopo, ax1); i1++){
                for (int i0 = 0; i0 < flups_topo_get_nloc(mytopo, ax0); i0++){
                    const size_t id = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                    local_err = std::max(local_err, abs(field[id] - sol[id]));
                }
            }
        }

        MPI_Allreduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        erri[i] = global_err;
        test_log("Test %d: err = %1.3e", i, erri[i]);
        
        flups_free(rhs);
        flups_free(sol);
        flups_free(field);
        
        flups_cleanup(mysolver);
        flups_topo_free(mytopo);
    }

    if(ntest == 0){
        test_log("These boundary conditions are not yet implementedfor this kernel: %s", kname.at(green_).c_str());
    }else if(ntest == 1){
        test_log(" The convergence of %s is spectral. You obtained a error of %e ", kname.at(green_).c_str(), erri[0]);
        ASSERT_NEAR(erri[0], 0.0, ZERO_TOL);
    }else{
        double expected_order = KernelOrder(green_);
        double computed_order = -log(erri[1] / erri[0]) / log(2.0);
        test_log(" The choosen kernel is %s has an order of %2.0f. You obtained a convergence of %2.3f", kname.at(green_).c_str(), expected_order, computed_order);
        ASSERT_GE(computed_order, 0.9*expected_order);
    }
    //--------------------------------------------------------------------
}

INSTANTIATE_TEST_SUITE_P(AllTest,
                         ConvergenceTest,
                         testing::Combine(testing::Values(NODE_CENTER, CELL_CENTER), 
                                          testing::Values(CHAT_2, LGF_2, LGF_4, LGF_6, LGF_8, HEJ_2, HEJ_4, HEJ_6, HEJ_8, HEJ_10, HEJ_0),
                                          testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
                                          testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
                                          testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD)),  // z boundary conditions
                         TestNameGenerator);

// INSTANTIATE_TEST_SUITE_P(Cell,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER), 
//                                           testing::Values(CHAT_2, LGF_2, HEJ_2, HEJ_4, HEJ_6, HEJ_8, HEJ_10, HEJ_0),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

//...........................................................................................................................................................
// Instantiate the test for cell-centered data
// INSTANTIATE_TEST_SUITE_P(CellCHAT2,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER), 
//                                           testing::Values(CHAT_2),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions


// INSTANTIATE_TEST_SUITE_P(CellLGF2,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(LGF_2),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(CellHEJ2,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(HEJ_2),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(CellHEJ4,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(HEJ_4),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(CellHEJ6,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(HEJ_6),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(CellHEJ8,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(HEJ_8),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions
// INSTANTIATE_TEST_SUITE_P(CellHEJ10,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(HEJ_10),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(CellHEJ0,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(CELL_CENTER),
//                                           testing::Values(HEJ_0),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions
//...........................................................................................................................................................
// Instantiate the test for node-centered data
// INSTANTIATE_TEST_SUITE_P(NodeCHAT2,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER), 
//                                           testing::Values(CHAT_2),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions


// INSTANTIATE_TEST_SUITE_P(NodeLGF2,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(LGF_2),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(NodeHEJ2,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(HEJ_2),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(NodeHEJ4,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(HEJ_4),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(NodeHEJ6,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(HEJ_6),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(NodeHEJ8,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(HEJ_8),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions
// INSTANTIATE_TEST_SUITE_P(NodeHEJ10,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(HEJ_10),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions

// INSTANTIATE_TEST_SUITE_P(NodeHEJ0,
//                          ConvergenceTest,
//                          testing::Combine(testing::Values(NODE_CENTER),
//                                           testing::Values(HEJ_0),
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // x boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD),   // y boundary conditions
//                                           testing::Values(PE_PE, EV_EV, OD_OD, UN_UN, EV_UN, UN_EV, EV_OD, OD_EV, OD_UN, UN_OD))); // z boundary conditions