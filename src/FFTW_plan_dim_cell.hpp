/**
 * @file FFTW_plan_dim_cell.hpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef FFTW_PLAN_DIM_CELL_HPP
#define FFTW_PLAN_DIM_CELL_HPP

#include "FFTW_plan_dim.hpp"

/**
 * @brief A FFTW plan in one dimension
 * 
 */
class FFTW_plan_dim_cell : public FFTW_plan_dim {
    protected:

    public: 
    FFTW_plan_dim_cell(const int lda, const int dimID, const double h[3], const double L[3], BoundaryType* mybc[2], const int sign, const bool isGreen):\
        FFTW_plan_dim(lda, dimID, h, L, mybc, sign, isGreen){};
    ~FFTW_plan_dim_cell(){}; 

        /**
     * @name Initialization
     */
    /**@{ */
    void init_real2real_(const int size[3], bool isComplex)    override;
    void init_mixunbounded_(const int size[3], bool isComplex) override;
    void init_periodic_(const int size[3], bool isComplex)     override;
    void init_unbounded_(const int size[3], bool isComplex)    override;
    void init_empty_(const int size[3], bool isComplex) override {};
    std::string disp_data_center() const override {return ("cell-centered");};
    /**@} */
     
};
# endif //FFTW_PLAN_DIM_CELL_HPP
