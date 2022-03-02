#include "base_convergence_test.hpp"

/**
 * @brief Return true if the green functions are fully spectral
 * 
 * @return true 
 * @return false 
 */
bool BaseConvergenceTest::IsBcSpectral(){
    bool isXSpectral = (*(mybc_[0][0]) != UNB) && (*(mybc_[0][1]) != UNB); 
    bool isYSpectral = (*(mybc_[1][0]) != UNB) && (*(mybc_[1][1]) != UNB);
    bool isZSpectral = (*(mybc_[2][0]) != UNB) && (*(mybc_[2][1]) != UNB);
    return (isXSpectral && isYSpectral && isZSpectral);
}

/**
 * @brief Perform the test based on an id 
 *        There is 1000 tests. Base on the test number, we first compute the corresponding 
 *        boundary conditions. 
 *        The tests are performed on all the kernels 
 *        We use the same tricks as for the validation. See caprace2020 for the description of the equations
 * 
 * @param case_id 
 */
void BaseConvergenceTest::PerformTest(int case_id){
    InitBoundaryConditions(case_id);
    PrintBcs();
    
    // Perform the test for the Chat_2 kernel
    for (int i = 0; i < 1; i++) {
        FLUPS_GreenType kernel = (FLUPS_GreenType)i;
        test_log("Testing kernel %s \n", kname[(int)kernel].c_str());
        DoMagic(kernel);
    }
    KillBoundaryConditions();
}

/**
 * @brief Init the boundary condition based on the id of the test
 * 
 * @param case_id 
 */
void BaseConvergenceTest::InitBoundaryConditions(int case_id){
    for (int id = 0; id < 3; id++)
    {
        for (int is = 0; is < 2; is++)
        {
            mybc_[id][is] = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * 1);
        }
    }
    InitBcFromId(case_id, mybc_);
}

/**
 * @brief Print the test boundary conditions. 
 * 
 */
void BaseConvergenceTest::PrintBcs(){
    test_log("The boundary conditions for this case are ");
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            test_log("%d ", *mybc_[id][is]);
        }
    }
    test_log("\n");
}

/**
 * @brief Do the actual test
 * 
 * @param kernel 
 */
void BaseConvergenceTest::DoMagic(FLUPS_GreenType kernel) {
    int Nmin = 64;
    if ((kernel == CHAT_2 || kernel == HEJ_0) && (IsBcSpectral())) {
        InitFlups_(Nmin, kernel);
        double err = ComputeErr_();
        test_log(" The convergence of %s is spectral. You obtained a error of %e \n \n", kname[(int)kernel].c_str(), err);
        ASSERT_NEAR(err, 0.0, ZERO_TOL);
        KillFlups_();
    } else if (!IsBcSpectral() && (kernel == HEJ_0 || kernel == LGF_2)) {
        test_log("Unbounded boundary conditions are not implemented yet for this kernel:  %s \n", kname[(int)kernel].c_str());
    } else { 
        double erri[2];
        for(int i = 1; i <= 2; i++){
            InitFlups_(Nmin*i, kernel);
            erri[i-1] = ComputeErr_();
            KillFlups_();
        }
        double expected_order = KernelOrder(kernel);
        double computed_order = -log(erri[1] / erri[0]) / log(2.0);
        test_log(" The choosen kernel is %s has an order of %2.0f. You obtained a convergence of %2.3f\n \n", kname[(int)kernel].c_str(), expected_order, computed_order);
        ASSERT_GE(computed_order, 0.9*expected_order);
    }
}

/**
 * @brief Free the vector of boundary conditions
 * 
 */
void BaseConvergenceTest::KillBoundaryConditions() {
#pragma unroll 3
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc_[id][is]);
        }
    }
}

/**
 * @brief Init all the structures required by flups
 * 
 * @param N 
 * @param kernel 
 */
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


/**
 * @brief Compute the error. 
 *        Here we solve \Nabla^2 sol = rhs. 
 *        The expression of the sol and the rhs are computed as explained in Caprace 2020, in validation
 * 
 * @return double 
 */
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

    double local_err = std::numeric_limits<double>::lowest();
    double global_err = std::numeric_limits<double>::lowest();
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo_, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo_, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo_, ax0); i0++) {
                const size_t id = flups_locID(0, i0, i1, i2, 0, 0, nmem, 1);
                local_err             = std::max(local_err, abs(field_[id] - sol_[id]));
            }
        }
    }

    MPI_Allreduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_err;
}

/**
 * @brief Free all the structures initiated by flups
 * 
 */
void BaseConvergenceTest::KillFlups_() {
    flups_free(rhs_);
    flups_free(sol_);
    flups_free(field_);

    flups_topo_free(topo_);
    flups_cleanup(mysolver_);
}
