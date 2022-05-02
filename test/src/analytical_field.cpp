#include "analytical_field.hpp"

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
    } else if (bc_[0] == UNB && bc_[1] == UNB) {
        sol_fn_ = &fUnbSpietz;
        rhs_fn_ = &d2dx2_fUnbSpietz;
    } else if (bc_[0] == UNB && bc_[1] != UNB){
        if(bc_[1] == ODD){
            center_  = 0.7;
            sign_[1] = -1;
        } else if (bc_[1] == EVEN){
            center_  = 0.7;
            sign_[1] = +1;
        }
        sol_fn_ = &fUnbSpietz;
        rhs_fn_ = &d2dx2_fUnbSpietz;
    } else if (bc_[1] == UNB && bc_[0] != UNB){
        if(bc_[0] == ODD){
            center_  = 0.3;
            sign_[0] = -1;
        } else if (bc_[0] == EVEN){
            center_  = 0.3;
            sign_[0] = +1;
        }
        sol_fn_ = &fUnbSpietz;
        rhs_fn_ = &d2dx2_fUnbSpietz;
    }
    
    if (sol_fn_ == nullptr || rhs_fn_ == nullptr) {
        printf("[WARNING] - You do not have any Rhs --> bc_[0] = %d -- bc_[1] == %d \n", bc_[0], bc_[1]);
    }
} 

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