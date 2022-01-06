#ifndef _SRC_BASECONVERGENCE_TEST_
#define _SRC_BASECONVERGENCE_TEST_

#include <limits>
#include "gtest/gtest.h"
#include "analytical_field.hpp"
#include "tools_test.hpp"
#include "flups.h"

#define DOUBLE_TOL 5*1e-15
#define CONV_TOL   0.15

#define NSPECTRAL 7
static const FLUPS_BoundaryType spectral_bc[NSPECTRAL][6] = {
    {PER, PER, PER, PER, PER, PER},
    {EVEN, EVEN, EVEN, EVEN, EVEN},
    {ODD, ODD, ODD, ODD, ODD, ODD},
    {EVEN, EVEN, ODD, ODD, EVEN, EVEN},
    {PER, PER, EVEN, EVEN, ODD, ODD},
    {PER, PER, ODD, EVEN, PER, PER},
    {EVEN, ODD, EVEN, EVEN, EVEN, EVEN}
};

static const FLUPS_BoundaryType fully_unbounded_bc[1][6] = {
    {UNB, UNB, UNB, UNB, UNB, UNB},
};

#define NMIXUNB 6
static const FLUPS_BoundaryType mix_unbounded_bc[NMIXUNB][6] = {
    {EVEN, EVEN, UNB, UNB, EVEN, EVEN},
    {UNB, UNB, UNB, UNB, UNB, UNB},
    {ODD, ODD, ODD, UNB, ODD, ODD},
    {EVEN, EVEN, EVEN, UNB, EVEN, EVEN},
    {ODD, ODD, UNB, ODD, ODD, ODD},
    {EVEN, EVEN, UNB, EVEN, EVEN, EVEN}
    
};


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
    const int    nproc_[3] = {1, 1, 1};
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

    void InitBoundaryConditions(int case_id){
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                mybc_[id][is]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
            }
        }
        InitBcFromId(case_id, mybc_);
    };

    void PrintBcs(){
        printf("The boundary conditions for this case are ");
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                printf("%d ", *mybc_[id][is]);
            }
        }
        printf("\n");
    }

    void DoMagic(FLUPS_GreenType kernel) {
        int Nmin = 64;
        if ((kernel == CHAT_2 || kernel == HEJ_0) && (IsBcSpectral())) {
            InitFlups_(Nmin, kernel);

            double err = ComputeErr_();
            printf(" The convergence of %s is spectral. You obtained a error of %e \n \n", kname[(int)kernel].c_str(), err);
            ASSERT_NEAR(err, 0.0, DOUBLE_TOL);
            KillFlups_();
        } else if (!IsBcSpectral() && (kernel == HEJ_0 || kernel == LGF_2)) {
            printf("Unbounded boundary conditions are not implemented yet for this kernel:  %s \n", kname[(int)kernel].c_str());
        } else { 
            double erri[2];
            for(int i = 1; i <= 2; i++){
                InitFlups_(Nmin*i, kernel);
                erri[i-1] = ComputeErr_();
                KillFlups_();
            }
            double expected_order = KernelOrder(kernel);
            double computed_order = -log(erri[1] / erri[0]) / log(2.0);
            printf(" The choosen kernel is %s has an order of %2.0f. You obtained a convergence of %2.3f\n \n", kname[(int)kernel].c_str(), expected_order, computed_order);
            ASSERT_GE(computed_order, expected_order - CONV_TOL);
        }
    }

    void KillBoundaryConditions() {
#pragma unroll 3
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                flups_free(mybc_[id][is]);
            }
        }
    }


    void PerformTest(int case_id){
        InitBoundaryConditions(case_id);
        PrintBcs();

        // printf("I will perform %d different cases \n", n_case);
        for (int i = 0; i < 8; i++) {
            FLUPS_GreenType kernel = (FLUPS_GreenType)i;
            printf("Testing kernel %s \n", kname[(int)kernel].c_str());
            DoMagic(kernel);
        }
        KillBoundaryConditions();
    }
};

#endif //_SRC_BASECONVERGENCE_TEST_