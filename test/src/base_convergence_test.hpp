#ifndef _SRC_BASECONVERGENCE_TEST_
#define _SRC_BASECONVERGENCE_TEST_

#include <limits>
#include "gtest/gtest.h"
#include "analytical_field.hpp"
#include "tools_test.hpp"
#include "flups.h"

#define ZERO_TOL 1000.0 * std::numeric_limits<double>::epsilon()
#define CONVERGENCE_TOL 0.4 //Great tolerance because of the slow convergence of the kernels.

#define test_log(format, ...)                              \
    ({                                                     \
        int test_log_rank_;                                \
        MPI_Comm_rank(MPI_COMM_WORLD, &test_log_rank_);    \
        char test_log_msg_[1024];                          \
        sprintf(test_log_msg_, format, ##__VA_ARGS__);     \
        if (test_log_rank_ == 0)                           \
        {                                                  \
            fprintf(stdout, "[murphy] %s", test_log_msg_); \
        }                                                  \
    })

// #define NSPECTRAL 7
// static const FLUPS_BoundaryType spectral_bc[NSPECTRAL][6] = {
//     {PER, PER, PER, PER, PER, PER},
//     {EVEN, EVEN, EVEN, EVEN, EVEN},
//     {ODD, ODD, ODD, ODD, ODD, ODD},
//     {EVEN, EVEN, ODD, ODD, EVEN, EVEN},
//     {PER, PER, EVEN, EVEN, ODD, ODD},
//     {PER, PER, ODD, EVEN, PER, PER},
//     {EVEN, ODD, EVEN, EVEN, EVEN, EVEN}
// };

// static const FLUPS_BoundaryType fully_unbounded_bc[1][6] = {
//     {UNB, UNB, UNB, UNB, UNB, UNB},
// };

// #define NMIXUNB 6
// static const FLUPS_BoundaryType mix_unbounded_bc[NMIXUNB][6] = {
//     {EVEN, EVEN, UNB, UNB, EVEN, EVEN},
//     {UNB, UNB, UNB, UNB, UNB, UNB},
//     {ODD, ODD, ODD, UNB, ODD, ODD},
//     {EVEN, EVEN, EVEN, UNB, EVEN, EVEN},
//     {ODD, ODD, UNB, ODD, ODD, ODD},
//     {EVEN, EVEN, UNB, EVEN, EVEN, EVEN}
    
// };

/**
 * @brief Return the order of the kernel 
 * 
 * @param kernel the asked kernel
 * @return double the order of the kernel in the arguments
 */
static double KernelOrder(FLUPS_GreenType kernel){
    switch (kernel) {
        case CHAT_2:
            // printf("Warning, you shouldn't ask for the convergence order of a spectral kernel\n");
            return 2.0;
            break;
        case LGF_2:
            return 2.0;
            break;
        case HEJ_2:
            return 2.0;
            break;
        case HEJ_4:
            return 4.0;
            break;
        case HEJ_6: 
            return 6.0;
            break;
        case HEJ_8:
            return 8.0;
            break;
        case HEJ_10:
            return 10.0;
            break;
        case HEJ_0:
            printf("Warning, you shouldn't ask for the convergence order of a spectral kernel\n");
            return std::numeric_limits<double>::max();
            break;
    }
    return 0; 
}


/**
 * @brief The base class for a convergence test. 
 * Both the Node Centered and the Cell Centered test uses the exact same functions to perform the tests. 
 * However, the function to compute the number of points must be overwritten
 * 
 */
class BaseConvergenceTest : public testing::TestWithParam<int> {
   protected:
    // Google test functions 
    void        SetUp() override{};
    void        TearDown() override{};

    // Vector 
    double          errinf_[2];
    double         *rhs_   = nullptr;
    double         *sol_   = nullptr;
    double         *field_ = nullptr;

    // Variable specific to the test_case
    const FLUPS_CenterType center_type_[3] = {CELL_CENTER, CELL_CENTER, CELL_CENTER};
    const double           shift_         = 0.5;

    // Flups related variables
    const int    nproc_[3] = {2, 2, 2};
    const double L_[3]     = {1., 1., 1.};
    int          nglob_[3] = {32, 32, 32};
    double       h_[3]     = {1. / 32., 1. / 32., 1. / 32.};

    FLUPS_BoundaryType* mybc_[3][2];
    FLUPS_Topology*     topo_     = nullptr;
    FLUPS_Solver*       mysolver_ = nullptr;
    std::string         kname[8]  = {"CHAT_2", "LGF_2", "HEJ_2", "HEJ_4", "HEJ_6", "HEJ_8", "HEJ_10", "HEJ_0"};

    // Functions used in the test - To be overriden
    bool   IsBcSpectral();
    void   InitFlups_(int N, FLUPS_GreenType kernel);
    double ComputeErr_();
    void   KillFlups_();

    virtual void SetProblemDef_(int N) = 0;

   public:
    BaseConvergenceTest(double shift, FLUPS_CenterType center) : shift_(shift), center_type_{center, center, center}{};

    void PerformTest(int case_id);

    // Functions needed to perform the tests
    void InitBoundaryConditions(int case_id);
    void PrintBcs();
    void DoMagic(FLUPS_GreenType kernel);
    void KillBoundaryConditions();
};

#endif //_SRC_BASECONVERGENCE_TEST_