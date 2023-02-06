#ifndef _SRC_ANALYTICAL_FIELD_
#define _SRC_ANALYTICAL_FIELD_

#include "h3lpr/profiler.hpp"
#include "flups.h"

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
        case LGF_4:
            return 4.0;
            break;
        case LGF_6:
            return 6.0;
            break;
        case LGF_8:
            return 8.0;
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

class AnalyticalField;
typedef double (*anal_fn)(const double, const double, const AnalyticalField* fid);

class AnalyticalField {
    protected:
     int                lia_     = 0;
     double             freq_    = 1;  // an integer or 0.5
     double             sign_[2] = {0., 0.};
     double             center_  = 0.5;
     double             sigma_   = .5;  // sigma for the compact Gaussian
     FLUPS_BoundaryType bc_[2]   = {PER, PER};
     anal_fn            sol_fn_  = nullptr;
     anal_fn            rhs_fn_  = nullptr;
     
     void SetAnalFn_();

    public: 
    explicit AnalyticalField() = default; 
    explicit AnalyticalField(FLUPS_BoundaryType bc[2]);
    explicit AnalyticalField(FLUPS_BoundaryType bc[2], int lda, double freq, double sign[2], double center, double sigma);//:bc_{bc[0], bc[1]}, lia_(lda), freq_(freq), sign_{sign[0], sign[1]}, center_(center), sigma_(sigma){}; 
    ~AnalyticalField(){};

    
    void SetParam(FLUPS_BoundaryType bc0, FLUPS_BoundaryType bc1){
        bc_[0] = bc0;
        bc_[1] = bc1;
        SetAnalFn_();
    }

    // Getter 
    double freq() const { return freq_; };
    double sign(const int id) const { return sign_[id]; };
    double center() const { return center_; };
    double sigma() const { return sigma_; };

    // Setter
    void SetFreq(const double freq){freq_ = freq;};
    void SetCenter(const double center){center_ = center;};
    void SetSign(const int id, const double sign){sign_[id] = sign;};

    double Rhs(const double x, const double L){
        return rhs_fn_(x, L, this);
    }; 
    double Sol(const double x, const double L){
        return sol_fn_(x, L, this);
    };
};



static inline double fOddOdd(const double x, const double L, const AnalyticalField* fid) {
    return sin((c_2pi / L * fid->freq()) * x);
}
static inline double d2dx2_fOddOdd(const double x, const double L, const AnalyticalField* fid) {
    return -(c_2pi / L * fid->freq()) * (c_2pi / L * fid->freq()) * sin((c_2pi / L * fid->freq()) * x);
}

static inline double fEvenEven(const double x, const double L, const AnalyticalField* fid) {
    return cos((M_PI / L * fid->freq()) * x);
}
static inline double d2dx2_fEvenEven(const double x, const double L, const AnalyticalField* fid) {
    return -(M_PI / L * fid->freq()) * (M_PI / L * fid->freq()) * cos((M_PI / L * fid->freq()) * x);
}

static inline double fOddEven(const double x, const double L, const AnalyticalField* fid) {
    return sin((M_PI / L * (fid->freq()+.5)) * x);
}
static inline double d2dx2_fOddEven(const double x, const double L, const AnalyticalField* fid) {
    return -(M_PI / L * (fid->freq()+.5)) * (M_PI / L * (fid->freq()+.5)) * sin((M_PI / L * (fid->freq()+.5)) * x);
}

static inline double fEvenOdd(const double x, const double L, const AnalyticalField* fid) {
    return cos((M_PI / L * (fid->freq()+.5)) * x);
}
static inline double d2dx2_fEvenOdd(const double x, const double L, const AnalyticalField* fid) {
    return -(M_PI / L * (fid->freq()+.5)) * (M_PI / L * (fid->freq()+.5)) * cos((M_PI / L * (fid->freq()+.5)) * x);
}

static inline double fUnb(const double x, const double L, const AnalyticalField* fid) {
    const double x0 = (x -       fid->center()  * L) / fid->sigma();
    const double x1 = (x +       fid->center()  * L) / fid->sigma();
    const double x2 = (x - (2. - fid->center()) * L) / fid->sigma();
    return                   exp(-x0*x0) + \
            fid->sign(0) * exp(-x1*x1) + \
            fid->sign(1) * exp(-x2*x2) ;
}
static inline double ddx_fUnb(const double x, const double L, const AnalyticalField* fid) {
    const double x0 = (x -       fid->center()  * L) / fid->sigma();
    const double x1 = (x +       fid->center()  * L) / fid->sigma();
    const double x2 = (x - (2. - fid->center()) * L) / fid->sigma();

    return                   -2.0/(fid->sigma()) * exp(-x0*x0) * x0 + \
            fid->sign(0)   * -2.0/(fid->sigma()) * exp(-x1*x1) * x1 + \
            fid->sign(1)   * -2.0/(fid->sigma()) * exp(-x2*x2) * x2;
}
static inline double d2dx2_fUnb(const double x, const double L, const AnalyticalField* fid) {
    const double x0 = (x -       fid->center()  * L) / fid->sigma();
    const double x1 = (x +       fid->center()  * L) / fid->sigma();
    const double x2 = (x - (2. - fid->center()) * L) / fid->sigma();
    return                   -2. / (fid->sigma() * fid->sigma()) * exp(-x0*x0 ) * (1. - 2. * ( x0 * x0)) + \
            fid->sign(0) * -2. / (fid->sigma() * fid->sigma()) * exp(-x1*x1 ) * (1. - 2. * ( x1 * x1)) + \
            fid->sign(1) * -2. / (fid->sigma() * fid->sigma()) * exp(-x2*x2 ) * (1. - 2. * ( x2 * x2)) ;
}

static const double c_C = 10.;

static inline double fUnbSpietz(const double x, const double L, const AnalyticalField* fid) {
    const double x0 = (x -       fid->center()  * L) / (fid->sigma() * L);
    const double x1 = (x +       fid->center()  * L) / (fid->sigma() * L);
    const double x2 = (x - (2. - fid->center()) * L) / (fid->sigma() * L);
    return    (fabs(x0)>=1. ? 0.0 :                  exp(c_C * (1. - 1. / (1. - x0 * x0)))) + \
              (fabs(x1)>=1. ? 0.0 : fid->sign(0) * exp(c_C * (1. - 1. / (1. - x1 * x1)))) + \
              (fabs(x2)>=1. ? 0.0 : fid->sign(1) * exp(c_C * (1. - 1. / (1. - x2 * x2))));
}
static inline double ddx_fUnbSpietz(const double x, const double L, const AnalyticalField* fid) {
    const double x0 = (x -       fid->center()  * L) / (fid->sigma() * L);
    const double x1 = (x +       fid->center()  * L) / (fid->sigma() * L);
    const double x2 = (x - (2. - fid->center()) * L) / (fid->sigma() * L);

    return    (fabs(x0)>=1. ? 0.0 :                  (-2.0) * c_C * x0 * exp(c_C * (1.0 - 1.0 / (1.0 - x0*x0))) / (fid->sigma()*L * pow((1.0 - x0*x0),2.0))) + \
              (fabs(x1)>=1. ? 0.0 : fid->sign(0) *   (-2.0) * c_C * x1 * exp(c_C * (1.0 - 1.0 / (1.0 - x1*x1))) / (fid->sigma()*L * pow((1.0 - x1*x1),2.0))) + \
              (fabs(x2)>=1. ? 0.0 : fid->sign(1) *   (-2.0) * c_C * x2 * exp(c_C * (1.0 - 1.0 / (1.0 - x2*x2))) / (fid->sigma()*L * pow((1.0 - x2*x2),2.0))) ;
}
static inline double d2dx2_fUnbSpietz(const double x, const double L, const AnalyticalField* fid) {
    const double x0sq = pow((x -       fid->center()  * L) / (fid->sigma() * L),2);
    const double x1sq = pow((x +       fid->center()  * L) / (fid->sigma() * L),2);
    const double x2sq = pow((x - (2. - fid->center()) * L) / (fid->sigma() * L),2);
    return    (fabs(x0sq)>=1. ? 0.0 :                  exp(c_C * (1. - 1. / (1. - x0sq))) * (2. * c_C * (2. * (c_C - 1.) * x0sq + 3. * x0sq * x0sq - 1.)) / pow(x0sq - 1., 4) / (fid->sigma()*L*fid->sigma()*L)) + \
              (fabs(x1sq)>=1. ? 0.0 : fid->sign(0) * exp(c_C * (1. - 1. / (1. - x1sq))) * (2. * c_C * (2. * (c_C - 1.) * x1sq + 3. * x1sq * x1sq - 1.)) / pow(x1sq - 1., 4) / (fid->sigma()*L*fid->sigma()*L)) + \
              (fabs(x2sq)>=1. ? 0.0 : fid->sign(1) * exp(c_C * (1. - 1. / (1. - x2sq))) * (2. * c_C * (2. * (c_C - 1.) * x2sq + 3. * x2sq * x2sq - 1.)) / pow(x2sq - 1., 4) / (fid->sigma()*L*fid->sigma()*L)) ;
}

static inline double fCst(const double x, const double L, const AnalyticalField* fid) {
    return 1.0;
}
static inline double fZero(const double x, const double L, const AnalyticalField* fid) {
    return 0.0;
}


#endif //_SRC_ANALYTICAL_FIELD_