#include "base_convergence_test.hpp"


bool BaseConvergenceTest::IsBcSpectral(){
    bool isXSpectral = (*(mybc_[0][0]) != UNB) && (*(mybc_[0][1]) != UNB); 
    bool isYSpectral = (*(mybc_[1][0]) != UNB) && (*(mybc_[1][1]) != UNB);
    bool isZSpectral = (*(mybc_[2][0]) != UNB) && (*(mybc_[2][1]) != UNB);
    return (isXSpectral && isYSpectral && isZSpectral);
};

void BaseConvergenceTest::InitFlups_(int N, FLUPS_GreenType kernel){
    //  -------------------------------------------------------------------------
    //  Definition of the problem
    //  -------------------------------------------------------------------------
    SetProblemDef_(N);

    //-------------------------------------------------------------------------
    /** - Initialize FLUPS */
    //-------------------------------------------------------------------------
    // create a real topology
    topo_ = flups_topo_new(0, 1, nglob_, nproc_, false, NULL, FLUPS_ALIGNMENT, MPI_COMM_WORLD);

    // solver creation and init
    mysolver_ = flups_init(topo_, mybc_, h_, L_, NOD, center_type_);
    flups_set_greenType(mysolver_, kernel);
    flups_setup(mysolver_, true);
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
    analytics[0].SetFreq(1.);
    analytics[1].SetFreq(2.);
    analytics[2].SetFreq(4.);

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

void BaseConvergenceTest::KillFlups_() {
    flups_free(rhs_);
    flups_free(sol_);
    flups_free(field_);

    flups_topo_free(topo_);
    flups_cleanup(mysolver_);
}
