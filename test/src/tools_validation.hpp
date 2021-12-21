#ifndef _SRC_TOOLS_VALIDATION
#define _SRC_TOOLS_VALIDATION

#include <limits>
#include "gtest/gtest.h"
#include "flups.h"

#define DOUBLE_TOL 1e-15
#define CONV_TOL   0.05


static const double c_1opi     = 1.0 / (1.0 * M_PI);
static const double c_1o2pi    = 1.0 / (2.0 * M_PI);
static const double c_1o4pi    = 1.0 / (4.0 * M_PI);
static const double c_1osqrtpi = 1.0 / sqrt(M_PI);
static const double c_1o2      = 1.0 / 2.0;
static const double c_1o4      = 1.0 / 4.0;
static const double c_7o4      = 7. / 4.;
static const double c_19o8     = 19. / 8;
static const double c_2o3      = 2. / 3;
static const double c_1o24     = 1. / 24;

static const double c_1osqrt2 = 1.0 / M_SQRT2;

static const double c_2pi = 2.0 * M_PI;

static double KernelOrder(FLUPS_GreenType kernel){
    switch (kernel) {
        case CHAT_2:
            printf("Warning, you shouldn't ask for the convergence order of a spectral kernel\n");
            return std::numeric_limits<double>::max();
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


class AnalyticalField {
    protected:
     int                lia_     = 0;
     double             freq_    = 1;  // an integer or 0.5
     double             sign_[2] = {0., 0.};
     double             center_  = 0.5;
     double             sigma_   = .5;  // sigma for the compact Gaussian
     FLUPS_BoundaryType bc_[2]   = {PER, PER};


    public: 
    explicit AnalyticalField() = default; 
    explicit AnalyticalField(FLUPS_BoundaryType bc[2]): bc_{bc[0], bc[1]} {}; 
    explicit AnalyticalField(FLUPS_BoundaryType bc[2], int lda, double freq, double sign[2], double center, double sigma):bc_{bc[0], bc[1]}, lia_(lda), freq_(freq), sign_{sign[0], sign[1]}, center_(center), sigma_(sigma){}; 
    ~AnalyticalField(){};

    void SetParam(FLUPS_BoundaryType bc[2]){
        bc_[0] = bc[0];
        bc_[1] = bc[1];
    }
    
    inline double freq() const {return freq_;};

    template <FLUPS_BoundaryType BC0, FLUPS_BoundaryType BC1>
    inline double RhsSpe(const double x, const double L){
        
        // FLUPS_CHECK(false, "There is no specialisation for your case", LOCATION);
    };

    template <FLUPS_BoundaryType BC0, FLUPS_BoundaryType BC1>
    inline double SolSpe(const double x, const double L){
        // FLUPS_CHECK(false, "There is no specialisation for your case", LOCATION);
    };
    double Rhs(const double x, const double L); 
    double Sol(const double x, const double L);
};

class BaseConvergenceTest : public testing::TestWithParam<int> {
   protected:
    // Google test functions 
    void        SetUp() override{};
    void        TearDown() override{};

    // Test case variables
    const double    shift_ = 0.5;
    double          errinf_[2];
    double         *rhs_   = nullptr;
    double         *sol_   = nullptr;
    double         *field_ = nullptr;

    // Flups related variables
    const double  L_[3]     = {1., 1., 1.};
    double h_[3]            = {1./32., 1./32., 1./32.};
    FLUPS_BoundaryType *mybc_[3][2];
    FLUPS_Topology     *topo_     = nullptr;
    FLUPS_Solver       *mysolver_ = nullptr;
    std::string         kname[8]  = {"CHAT_2", "LGF_2", "HEJ_2", "HEJ_4", "HEJ_6", "HEJ_8", "HEJ_10", "HEJ_0"};

    void KillFlups_();
    double ComputeErr_(); 

    // Functions used in the test - To be overriden
    virtual void InitFlups_(int N, FLUPS_GreenType kernel, FLUPS_BoundaryType bc_in) = 0;
   public:
    BaseConvergenceTest(double shift) : shift_(shift){};

    void ActualTest(FLUPS_GreenType kernel, FLUPS_BoundaryType boundary) {
        int Nmin = 128;
        if (kernel == CHAT_2 || kernel == HEJ_0) {
            InitFlups_(Nmin, kernel, boundary);

            double err = ComputeErr_();
            printf(" The convergence of %s is spectral. You obtained a error of %e \n \n", kname[(int)kernel].c_str(), err);
            ASSERT_NEAR(err, 0.0, DOUBLE_TOL);
        } else {
            double erri[2];
            for(int i = 1; i <= 2; i++){
                InitFlups_(Nmin*i, kernel, boundary);
                erri[i] = ComputeErr_();
                KillFlups_();
            }
            double expected_order = KernelOrder(kernel);
            double computed_order = -log(erri[1] / erri[0]) / log(2.0);
            printf(" The choosen kernel is %s has an order of %2.0f. You obtained a convergence of %2.3f \n \n", kname[(int)kernel].c_str(), expected_order, computed_order);
            ASSERT_GE(computed_order, expected_order - CONV_TOL);
        }
    }

};

template<> 
inline double AnalyticalField::RhsSpe<PER, PER>(const double x, const double L){
    return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
};

template<> 
inline double AnalyticalField::SolSpe<PER, PER>(const double x, const double L){
    return sin((c_2pi / L * freq_) * x);
};

template<> 
inline double AnalyticalField::RhsSpe<ODD, ODD>(const double x, const double L){
    return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
};

template<> 
inline double AnalyticalField::SolSpe<ODD, ODD>(const double x, const double L){
    return sin((c_2pi / L * freq_) * x);
};

template<> 
inline double AnalyticalField::RhsSpe<EVEN, EVEN>(const double x, const double L){
    return  -(M_PI / L * freq_) * (M_PI / L * freq_) * cos((M_PI / L * freq_) * x);
};

template<> 
inline double AnalyticalField::SolSpe<EVEN, EVEN>(const double x, const double L){
    return sin((M_PI / L * (freq_ + .5)) * x);
};

template<> 
inline double AnalyticalField::RhsSpe<ODD, EVEN>(const double x, const double L){
    return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
};

template<> 
inline double AnalyticalField::SolSpe<ODD, EVEN>(const double x, const double L){
    return -(M_PI / L * (freq_ + .5)) * (M_PI / L * (freq_ + .5)) * sin((M_PI / L * (freq_ + .5)) * x);
};

#endif //_SRC_TOOLS_VALIDATION