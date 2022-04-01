#ifndef SRC_TOOLS_INTERFACE_H_
#define SRC_TOOLS_INTERFACE_H_
#include "Solver.hpp"
#include "FFTW_plan_dim_cell.hpp"
#include "FFTW_plan_dim_node.hpp"


int hint_proc_repartition(const int lda, const double h[3], const double L[3], BoundaryType* bc[3][2], const CenterType center_type[3]){
    FFTW_plan_dim* plan[3];
    for (int id = 0; id < 3; id++) {
        if (CELL_CENTER == center_type[id]) {
            plan[id]  = new FFTW_plan_dim_cell(lda, id, h, L, bc[id], FLUPS_FORWARD, false);
        } else if (NODE_CENTER == center_type[id]) {
            plan[id]  = new FFTW_plan_dim_node(lda, id, h, L, bc[id], FLUPS_FORWARD, false);
        } else {
            FLUPS_CHECK(false, "The type of data you asked is not supported");
        }
    }
    
    int ridx = sort_plans(plan); 

    for (int id = 0; id < 3; id++) {
        delete plan[id];
    }
    return ridx;
}

#endif 