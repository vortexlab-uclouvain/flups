#include "tools_validation.hpp"


 void BaseConvergenceTest::KillFlups_() {
        flups_free(rhs_);
        flups_free(sol_);
        flups_free(field_);

        flups_topo_free(topo_);
        flups_cleanup(mysolver_);

#pragma unroll 3
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                flups_free(mybc_[id][is]);
            }
        }
}

double BaseConvergenceTest::ComputeErr_() {
    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    rhs_   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo_));
    sol_   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo_));
    field_ = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo_));
    memset(rhs_, 0, sizeof(double) * flups_topo_get_memsize(topo_));
    memset(sol_, 0, sizeof(double) * flups_topo_get_memsize(topo_));
    memset(field_, 0, sizeof(double) * flups_topo_get_memsize(topo_));


    int istart[3];
    flups_topo_get_istartGlob(topo_, istart);

    const int nmem[3] = {flups_topo_get_nmem(topo_, 0), flups_topo_get_nmem(topo_, 1), flups_topo_get_nmem(topo_, 2)};
    const int ax0     = flups_topo_get_axis(topo_);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;

    // const double cos_cons = 2.0 * M_PI * 3;
    AnalyticalField analytics[3];
    for (int dir = 0; dir < 3; dir++) {
        analytics[dir].SetParam(*mybc_[dir][0], *mybc_[dir][1]);
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo_, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo_, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo_, ax0); i0++) {
                const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                sol_[id]         = 1.0;
            }
        }
    }

    const double cos_cons = 2.5 * M_PI;
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo_, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo_, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo_, ax0); i0++) {
                const size_t id   = flups_locID(0, i0, i1, i2, 0, ax0, nmem, 1);
                const double x[3] = {(istart[ax0] + i0 + shift_) * h_[ax0],
                                     (istart[ax1] + i1 + shift_) * h_[ax1],
                                     (istart[ax2] + i2 + shift_) * h_[ax2]};

                // second derivative of a cosine
                // rhs_[id] =  cos(cos_cons * x[0]);
                for (int dir0 = 0; dir0 < 3; dir0++) {
                    const int dir1 = (dir0 + 1) % 3;
                    const int dir2 = (dir0 + 2) % 3;
                    sol_[id] *= analytics[dir0].Sol(x[dir0], L_[dir0]);
                    rhs_[id] += analytics[dir0].Rhs(x[dir0], L_[dir0]) * analytics[dir1].Sol(x[dir1], L_[dir1]) * analytics[dir2].Sol(x[dir2], L_[dir2]);
                }
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - solve the equations */
    //-------------------------------------------------------------------------
    for (int is = 0; is < 1; is++) {
        flups_solve(mysolver_, field_, rhs_, STD);
    }

    double err = std::numeric_limits<double>::lowest();
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo_, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo_, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo_, ax0); i0++) {
                const size_t id = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                err             = std::max(err, abs(field_[id] - sol_[id]));
            }
        }
    }
    return err;
}


void AnalyticalField::SetAnalFn_(){
    if ((bc_[0] == PER && bc_[1] == PER) || (bc_[0] == ODD && bc_[1] == ODD)) {
        sol_fn_ = &fOddOdd;
        rhs_fn_ = &d2dx2_fOddOdd;

    } else if (bc_[0] == EVEN && bc_[1] == EVEN) {
        sol_fn_ = &fEvenEven;
        rhs_fn_ = &d2dx2_fEvenEven;
    } else if (bc_[0] == EVEN && bc_[1] == ODD) {
        sol_fn_ = &fEvenOdd;
        rhs_fn_ = &d2dx2_fEvenOdd;
    } else if (bc_[0] == ODD && bc_[1] == EVEN) {
        sol_fn_ = &fOddEven;
        rhs_fn_ = &d2dx2_fOddEven;
    }  else if (bc_[0] == UNB && bc_[1] == UNB) {
        sol_fn_ = &fOddEven;
        rhs_fn_ = &d2dx2_fOddEven;
    }
    if (sol_fn_ == nullptr || rhs_fn_ == nullptr) {
        printf("[WARNING] - You do not have any Rhs --> bc_[0] = %d -- bc_[1] == %d \n", bc_[0], bc_[1]);
    }
}; 

// AnalyticalField::AnalyticalField(FLUPS_BoundaryType bc[2], int lda, double freq, double sign[2], double center, double sigma): AnalyticalField(bc), lia_(lda), freq_(freq), sign_{sign[0], sign[1]}, center_(center), sigma_(sigma){}; 


// double AnalyticalField::Rhs(const double x, const double L) {
//     if((bc_[0] == PER && bc_[1] == PER) || (bc_[0] == ODD && bc_[1] == ODD )){
//         return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
    
//     }else if(bc_[0] == EVEN && bc_[1] == EVEN){
//         return  -(M_PI / L * freq_) * (M_PI / L * freq_) * cos((M_PI / L * freq_) * x);
    
//     }else if(bc_[0] == EVEN && bc_[1] == ODD){
//         return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
//     }

//     printf("[WARNING] - You do not have any Rhs \n");
//     return std::numeric_limits<double>::max();
// }

// double AnalyticalField::Sol(const double x, const double L) {
//     if((bc_[0] == PER && bc_[1] == PER) || (bc_[0] == ODD && bc_[1] == ODD )){
//         return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
    
//     }else if(bc_[0] == EVEN && bc_[1] == EVEN){
//         return  -(M_PI / L * freq_) * (M_PI / L * freq_) * cos((M_PI / L * freq_) * x);
    
//     }else if(bc_[0] == EVEN && bc_[1] == ODD){
//         return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
//     }

//     printf("[WARNING] - You do not have any Rhs \n");
//     return std::numeric_limits<double>::max();
// }

// template<> 
// inline double AnalyticalField::SolSpe<PER, PER>(const double x, const double L){
//     return sin((c_2pi / L * freq_) * x);
// };


// template<> 
// inline double AnalyticalField::SolSpe<ODD, ODD>(const double x, const double L){
//     return sin((c_2pi / L * freq_) * x);
// };

// template<> 
// inline double AnalyticalField::RhsSpe<EVEN, EVEN>(const double x, const double L){
    
// };

// template<> 
// inline double AnalyticalField::SolSpe<EVEN, EVEN>(const double x, const double L){
//     return sin((M_PI / L * (freq_ + .5)) * x);
// };

// template<> 
// inline double AnalyticalField::RhsSpe<ODD, EVEN>(const double x, const double L){
//     return -(c_2pi / L * freq_) * (c_2pi / L * freq_) * sin((c_2pi / L * freq_) * x);
// };

// template<> 
// inline double AnalyticalField::SolSpe<ODD, EVEN>(const double x, const double L){
//     return -(M_PI / L * (freq_ + .5)) * (M_PI / L * (freq_ + .5)) * sin((M_PI / L * (freq_ + .5)) * x);
// };
// }