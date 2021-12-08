/**
 * @file FFTW_plan_dim.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef FFTW_PLAN_DIM_NODE_HPP
#define FFTW_PLAN_DIM_NODE_HPP

#include "FFTW_plan_dim.hpp"

/**
 * @brief A FFTW plan in one dimension
 * 
 */
class FFTW_plan_dim_node : public FFTW_plan_dim {
    protected:

    public: 
    FFTW_plan_dim_node(const int lda, const int dimID, const double h[3], const double L[3], BoundaryType* mybc[2], const int sign, const bool isGreen):\
        FFTW_plan_dim(lda, dimID, h, L, mybc, sign, isGreen){};
    ~FFTW_plan_dim_node(){}; 

        /**
     * @name Initialization
     */
    /**@{ */
    void init_real2real_(const int size[3], bool isComplex)    override;
    void init_mixunbounded_(const int size[3], bool isComplex) override;
    void init_periodic_(const int size[3], bool isComplex)     override;
    void init_unbounded_(const int size[3], bool isComplex)    override;
    void init_empty_(const int size[3], bool isComplex)        override;
    std::string disp_data_center() const override {return ("node-centered");} ;
    /**@} */
     
};
# endif //FFTW_PLAN_DIM_NODE_HPP
