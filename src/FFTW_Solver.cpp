/**
 * @file FFTW_Solver.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "FFTW_Solver.hpp"

/**
 * @brief Construct a fftw Poisson solver
 * 
 * @param topo the current topology of the data
 * @param mybc the boundary conditions of the solver
 * @param h the grid spacing
 * @param L 
 */
FFTW_Solver::FFTW_Solver(const Topology *topo, const BoundaryType mybc[DIM][2], const double h[3], const double L[3]) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Store the field size */
    //-------------------------------------------------------------------------
    // for(int id=0; id<DIM; id++) _size_field[id] = topo->nglob(id);
    // topo->disp();

    //-------------------------------------------------------------------------
    /** - For each dim, create the plans and sort them type */
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++)
        _hgrid[id] = h[id];

    for (int id = 0; id < DIM; id++) {
        _plan_forward[id]  = new FFTW_plan_dim(id, h, L, mybc[id], UP_FORWARD, false);
        _plan_backward[id] = new FFTW_plan_dim(id, h, L, mybc[id], UP_BACKWARD, false);
        _plan_green[id]    = new FFTW_plan_dim(id, h, L, mybc[id], UP_FORWARD, true);
    }

    _sort_plan(_plan_forward);
    _sort_plan(_plan_backward);
    _sort_plan(_plan_green);

    //-------------------------------------------------------------------------
    /** - Initialise the plans and get the sizes */
    //-------------------------------------------------------------------------
    _init_plan(topo, _topo_hat, _switchtopo, _plan_forward, false);
    _init_plan(topo, NULL, NULL, _plan_backward, false);
    _init_plan(topo, _topo_green, _switchtopo_green, _plan_green, true);

    //-------------------------------------------------------------------------
    /** - Get the factors #_normfact, #_volfact, #_shiftgreen and #_nbr_imult */
    //-------------------------------------------------------------------------
    _normfact  = 1.0;
    _volfact   = 1.0;
    _nbr_imult = 0;
    for (int ip = 0; ip < 3; ip++) {
        _normfact *= _plan_forward[ip]->normfact();
        _volfact *= _plan_forward[ip]->volfact();

        _shiftgreen[_plan_forward[ip]->dimID()] = _plan_forward[ip]->shiftgreen();

        if (_plan_forward[ip]->imult())
            _nbr_imult++;
        if (_plan_backward[ip]->imult())
            _nbr_imult--;
        if (_plan_green[ip]->imult())
            _nbr_imult++;
    }
}

/**
 * @brief Sets up the Solver
 * 
 * After this function the parameter of the solver (size etc) cannot be changed anymore
 * 
 * -------------------------------------------
 * We do the following operations
 */
void FFTW_Solver::setup() {
    //-------------------------------------------------------------------------
    /** - allocate the data for the field and Green */
    //-------------------------------------------------------------------------
    _allocate_data(_topo_hat, &_data);
    _allocate_data(_topo_green, &_green);

    //-------------------------------------------------------------------------
    /** - allocate the plans forward and backward for the field */
    //-------------------------------------------------------------------------
    _allocate_plan(_topo_hat, _plan_forward, _data);
    _allocate_plan(_topo_hat, _plan_backward, _data);

    //-------------------------------------------------------------------------
    /** - allocate the plan and comnpute the Green's function */
    //-------------------------------------------------------------------------
    _allocate_plan(_topo_green, _plan_green, _green);
    _cmptGreenFunction(_topo_green, _green, _plan_green);

    //-------------------------------------------------------------------------
    /** - delete the useless data for Green */
    //-------------------------------------------------------------------------
    _delete_plan(_plan_green);
    _delete_switchtopo(_switchtopo_green);
}

/**
 * @brief Destroy the fftw solver
 * 
 */
FFTW_Solver::~FFTW_Solver() {
    BEGIN_FUNC
    // delete plans
    _delete_plan(_plan_forward);
    _delete_plan(_plan_backward);

    // delete datas
    if (_green != NULL)
        fftw_free(_green);
    if (_data != NULL)
        fftw_free(_data);

    // delete switchtopo
    for (int id = 0; id < 3; id++) {
        if (_switchtopo[id] != NULL)
            delete _switchtopo[id];
    }

    //cleanup
    fftw_cleanup();
}
/**
 * @brief delete the FFTW_plan_dim stored in planmap
 * 
 * @param planmap 
 */
void FFTW_Solver::_delete_plan(FFTW_plan_dim *planmap[3]) {
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++) {
        delete planmap[ip];
        planmap[ip] = NULL;
    }
}

void FFTW_Solver::_delete_switchtopo(SwitchTopo *switchtopo[3]) {
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++) {
        delete switchtopo[ip];
        switchtopo[ip] = NULL;
    }
}
void FFTW_Solver::_delete_topology(Topology *topo[3]) {
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++) {
        delete topo[ip];
        topo[ip] = NULL;
    }
}

void FFTW_Solver::_sort_plan(FFTW_plan_dim *plan[3]) {
    int priority[3];

    for (int id = 0; id < 3; id++)
        priority[id] = plan[id]->type();

    // do the sort by hand...
    if (priority[0] > priority[1]) {
        FFTW_plan_dim *temp_plan = plan[1];
        plan[1]                  = plan[0];
        plan[0]                  = temp_plan;
    }
    // we are sure to have the first to sorted
    // if the last one is smaller than the second one, we change, if not, it is also bigger than the first one
    if (priority[1] > priority[2]) {
        FFTW_plan_dim *temp_plan = plan[2];
        plan[2]                  = plan[1];
        plan[1]                  = temp_plan;
        if (priority[0] > priority[1]) {
            FFTW_plan_dim *temp_plan = plan[1];
            plan[1]                  = plan[0];
            plan[0]                  = temp_plan;
        }
    }
}

/**
 * @brief Initialize a set of 3 plans by doing a dry run through the plans
 * 
 * @param topo the starting topology
 * @param topomap the topology array to go through each dim ( may be NULL) it corresponds to the topology AFTER the plan
 * @param switchtopo the switchtopoing array to switch between topologies (may be NULL)
 * @param planmap the plan that will be created
 * @param isGreen indicates if the plans are for Green
 */
void FFTW_Solver::_init_plan(const Topology *topo, Topology *topomap[3], SwitchTopo *switchtopo[3], FFTW_plan_dim *planmap[3], bool isGreen) {
    BEGIN_FUNC

    //-------------------------------------------------------------------------
    /** - Store the current topology */
    //-------------------------------------------------------------------------
    const Topology *current_topo = topo;

    //-------------------------------------------------------------------------
    /** - Get the sizes to start with */
    //-------------------------------------------------------------------------
    int size_tmp[DIM];
    for (int id = 0; id < DIM; id++)
        size_tmp[id] = topo->nglob(id);

    //-------------------------------------------------------------------------
    /** - create the plan (and the topologies if not Green) */
    //-------------------------------------------------------------------------
    bool isComplex = false;
    int  nproc[3];
    for (int ip = 0; ip < 3; ip++) {
        // initialize the plan
        planmap[ip]->init(size_tmp, isComplex);
        // update the size_tmp variable and get the complex information
        planmap[ip]->get_outsize(size_tmp);
        // virtually execute the plan
        planmap[ip]->get_isNowComplex(&isComplex);

        // we store a new topology BEFORE the plan is executed
        if (!isGreen && topomap != NULL && switchtopo != NULL) {
            // get the fastest rotating index
            int dimID = planmap[ip]->dimID();  // store the correspondance of the transposition
            // get the proc repartition
            _pencil_nproc(dimID, nproc, topo->comm_size());
            // create the new topology in the output layout (size and isComplex)
            topomap[ip] = new Topology(dimID, size_tmp, nproc, isComplex);
            // get the fieldstart = the point where the old topo has to begin in the new
            int fieldstart[3] = {0};
            planmap[ip]->get_fieldstart(fieldstart);
            // compute the transfert between the current topo and the new one
            // if the topo was real before the plan and is now complex
            if (planmap[ip]->isr2c()) {
                topomap[ip]->switch2real();
                switchtopo[ip] = new SwitchTopo(current_topo, topomap[ip], fieldstart);
                topomap[ip]->switch2complex();
            } else {
                // create the switchtopoMPI to change topology
                switchtopo[ip] = new SwitchTopo(current_topo, topomap[ip], fieldstart);
            }
            // update the current topo to the new one
            current_topo = topomap[ip];
        }

        planmap[ip]->disp();
    }

    //-------------------------------------------------------------------------
    /** - Do the magic trick  */
    //-------------------------------------------------------------------------
    // reset the correct input size for Green if the plan is R2C
    // only the R2C changes the size between INPUT and OUTPUT
    // this change has to be done ONLY if the first transform is real <=> is not R2C
    if (isGreen) {
        for (int ip = 0; ip < 3; ip++) {
            int dimID = planmap[ip]->dimID();
            if (!planmap[0]->isr2c() && planmap[ip]->isr2c()) size_tmp[dimID] *= 2;
        }
    }

    //-------------------------------------------------------------------------
    /** - For Green we need to compute the topologies using the full size of the domain  */
    //-------------------------------------------------------------------------
    isComplex = false;
    if (isGreen && topomap != NULL && switchtopo != NULL) {
        for (int ip = 0; ip < 3; ip++) {
            // virtually execute the plan
            planmap[ip]->get_isNowComplex(&isComplex);
            // get the fastest rotating index
            int dimID = planmap[ip]->dimID();  // store the correspondance of the transposition
            // get the proc repartition
            _pencil_nproc(dimID, nproc, topo->comm_size());
            // create the new topology in the output layout (size and isComplex)
            topomap[ip] = new Topology(dimID, size_tmp, nproc, isComplex);
            // get the fieldstart = the point where the old topo has to begin in the new
            int fieldstart[3] = {0};
            planmap[ip]->get_fieldstart(fieldstart);
            // compute the transfert between the current topo and the new one
            // if the topo was real before the plan and is now complex
            if (planmap[ip]->isr2c()) {
                topomap[ip]->switch2real();
                switchtopo[ip] = new SwitchTopo(current_topo, topomap[ip], fieldstart);
                topomap[ip]->switch2complex();
            } else {
                // create the switchtopoMPI to change topology
                switchtopo[ip] = new SwitchTopo(current_topo, topomap[ip], fieldstart);
            }
            // update the current topo to the new one
            current_topo = topomap[ip];
        }
    }

    //-------------------------------------------------------------------------
    /** - reset the topologies to real if needed  */
    //-------------------------------------------------------------------------
    for (int ip = 0; ip < 3; ip++) {
        if (planmap[ip]->isr2c() && topomap != NULL) {
            topomap[ip]->switch2real();
        }
    }
}
void FFTW_Solver::_allocate_plan(const Topology *const topo[3], FFTW_plan_dim *planmap[3], double *data) {
    BEGIN_FUNC

    for (int ip = 0; ip < 3; ip++) {
        UP_CHECK2(!(planmap[ip]->isr2c() && topo[ip]->isComplex()), "The topology %d need to be reset to the state BEFORE the plan to have the correct sizes for allocation (isComplex=%d)", ip, topo[ip]->isComplex());
        int size_plan[3] = {topo[ip]->nloc(0), topo[ip]->nloc(1), topo[ip]->nloc(2)};
        planmap[ip]->allocate_plan(size_plan, data);
    }
}

/**
 * @brief 
 * 
 * @param topo_hat the topologies of the pencils
 * @param data 
 */
void FFTW_Solver::_allocate_data(const Topology *const topo_hat[3], double **data) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0((*data) == NULL, "Pointer has to be NULL for allocation");

    //-------------------------------------------------------------------------
    /** - Do the memory allocation */
    //-------------------------------------------------------------------------
    // the bigger size will be in the pencils
    size_t size_tot = 1;
    for (int id = 0; id < 3; id++)
        size_tot = std::max(topo_hat[id]->locmemsize(), size_tot);

    INFOLOG2("Complex memory allocation, size = %ld\n", size_tot);
    (*data) = (double *)fftw_malloc(size_tot * sizeof(fftw_complex));

    //-------------------------------------------------------------------------
    /** - Check memory alignement */
    //-------------------------------------------------------------------------
    UP_CHECK1(UP_ISALIGNED(*data), "FFTW alignement not compatible with UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
}

/**
 * @brief compute the Green's function
 * 
 * @param topo 
 * @param green 
 * @param planmap 
 * 
 * -----------------------------------
 * We do the following operations
 */
void FFTW_Solver::_cmptGreenFunction(Topology *topo[3], double *green, FFTW_plan_dim *planmap[3]) {
    BEGIN_FUNC

    //-------------------------------------------------------------------------
    /** - get the direction where we need to do spectral diff and count them */
    //-------------------------------------------------------------------------
    bool dospectral[3] = {false};

    double hfact[3];
    double kfact[3];
    int    symstart[3];

    for (int ip = 0; ip < 3; ip++) {
        const int dimID = planmap[ip]->dimID();
        // get usefull datas
        dospectral[dimID] = planmap[ip]->dospectral();
        symstart[dimID]   = planmap[ip]->symstart();
        hfact[dimID]      = _hgrid[dimID];
        kfact[dimID]      = 0.0;

        if (dospectral[dimID]) {
            hfact[dimID] = 0.0;
            kfact[dimID] = planmap[ip]->kfact();
        }
    }

    // count the number of spectral dimensions
    int nbr_spectral = 0;
    for (int id = 0; id < 3; id++)
        if (dospectral[id])
            nbr_spectral++;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //-------------------------------------------------------------------------
    /** - get the expression of Green in the full domain*/
    //-------------------------------------------------------------------------
    if (nbr_spectral == 0) {
        INFOLOG(">> using Green function 3 dir unbounded\n");
        if (DIM == 3) {
            Green_3D_3dirunbounded_0dirspectral(topo[0], hfact, symstart, green, _typeGreen, _alphaGreen);
        }
    } else if (nbr_spectral == 1) {
        INFOLOG(">> using Green function 2 dir unbounded - 1 dir spectral\n");
        // _compute_Green_2dirunbounded_1dirspectral
        UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _typeGreen, nbr_spectral);
    } else if (nbr_spectral == 2) {
        INFOLOG(">> using Green function 1 dir unbounded - 2 dir spectral\n");
        // _compute_Green_1dirunbounded_2dirspectral
        UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _typeGreen, nbr_spectral);
    } else if (nbr_spectral == 3) {
        INFOLOG(">> using Green function 3 dir spectral\n");
        // _compute_Green_0dirunbounded_3dirspectral
        UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _typeGreen, nbr_spectral);
    }

    hdf5_dump(topo[0], "green", green);

    //-------------------------------------------------------------------------
    /** - scale the Green data using #_volfact */
    //-------------------------------------------------------------------------
    _scaleGreenFunction(topo[0], green);

    //-------------------------------------------------------------------------
    /** - compute a symmetry and do the forward transform*/
    //-------------------------------------------------------------------------
    for (int ip = 0; ip < 3; ip++) {
        // go to the topology for the plan
        if (ip > 0) {
            _switchtopo_green[ip]->execute(_green, UP_FORWARD);
        }
        // execute the plan
        _plan_green[ip]->execute_plan();

        if (_plan_green[ip]->isr2c()) {
            topo[ip]->switch2complex();
        }
    }
}

/**
 * @brief scales the Green's function given the #_volfact factor
 * 
 * @param topo the current topo
 * @param data the Green's function
 */
void FFTW_Solver::_scaleGreenFunction(const Topology *topo, double *data) {
    // the symmetry is done along the fastest rotating index
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                const size_t id = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
                data[id]        = data[id] * _volfact;
            }
        }
    }
}

/**
 * @brief Solve the Poisson equation
 * 
 * @param field 
 * @param rhs 
 * 
 * -----------------------------------------------
 * We perform the following operations:
 */
void FFTW_Solver::solve(const Topology *topo, double *field, double *rhs, const SolverType type) {
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0(field != NULL, "field is NULL");
    UP_CHECK0(rhs != NULL, "rhs is NULL");
    UP_CHECK1(UP_ISALIGNED(field), "pointer no aligned to UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(rhs), "pointer no aligned to UP_ALIGNMENT (=%d)", UP_ALIGNMENT);

    opt_double_ptr       myfield = field;
    opt_double_ptr       mydata  = (double *)_data;
    const opt_double_ptr myrhs   = rhs;

    //-------------------------------------------------------------------------
    /** - clean the data memory */
    //-------------------------------------------------------------------------
    // reset at the max size
    size_t size_tot = topo->locmemsize();
    for (int id = 0; id < 3; id++)
        size_tot = std::max(_topo_hat[id]->locmemsize(), size_tot);
    std::memset(mydata, 0, sizeof(double) * size_tot);

    //-------------------------------------------------------------------------
    /** - copy the rhs in the correct order */
    //-------------------------------------------------------------------------
    // INFOLOG("------------------------------------------\n");
    // INFOLOG("## memory information\n")
    // INFOLOG4("- size field   = %d %d %d\n", _size_field[0], _size_field[1], _size_field[2]);
    // INFOLOG4("- size hat     = %d %d %d\n", _size_hat[0], _size_hat[1], _size_hat[2]);
    // INFOLOG4("- dim order    = %d %d %d\n", _dimorder[0], _dimorder[1], _dimorder[2]);
    // INFOLOG4("- field start  = %d %d %d\n", _fieldstart[0], _fieldstart[1], _fieldstart[2]);
    // INFOLOG4("- dim multfact = %d %d %d\n", _dim_multfact[0], _dim_multfact[1], _dim_multfact[2]);
    // INFOLOG2("- offset       = %ld\n", _offset);
    // INFOLOG("------------------------------------------\n");

    int ax0 = topo->axis();
    int ax1 = (ax0 + 1) % 3;
    int ax2 = (ax0 + 2) % 3;
    if (topo->nf() == 1) {
        for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                    // comnpute the index permutation
                    const size_t id = localindex_ao(i0, i1, i2, topo);
                    // put the data
                    mydata[id] = myrhs[id];
                }
            }
        }
    } else if (topo->nf() == 2) {
        for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                    // comnpute the index permutation
                    const size_t id = localindex_ao(i0, i1, i2, topo);
                    // put the data
                    mydata[id + 0] = myrhs[id + 0];
                    mydata[id + 1] = myrhs[id + 1];
                }
            }
        }
    } else {
        UP_CHECK0(false, "size of Topological nf not supported");
    }

    hdf5_dump(topo, "rhs", mydata);

    //-------------------------------------------------------------------------
    /** - go to Fourier */
    //-------------------------------------------------------------------------
    for (int ip = 0; ip < 3; ip++) {
        // go to the correct topo
        _switchtopo[ip]->execute(mydata, UP_FORWARD);
        // run the FFT
        _plan_forward[ip]->execute_plan();
        // get if we are now complex
        if (_plan_forward[ip]->isr2c()) {
            _topo_hat[ip]->switch2complex();
        }
    }

    hdf5_dump(_topo_hat[2], "rhs_h", mydata);

    //-------------------------------------------------------------------------
    /** - Perform the magic */
    //-------------------------------------------------------------------------
    if (type == UP_SRHS) {
        if (!_topo_hat[2]->isComplex()) {
            // if (_nbr_imult == 0)
            //     dothemagic_rhs_real();
            // else
            UP_CHECK1(false, "the number of imult = %d is not supported", _nbr_imult);
        } else {
            if (_nbr_imult == 0)
                dothemagic_rhs_complex_nmult0();
            // else if(_nbr_imult == 1) dothemagic_rhs_complex_nmult1();
            // else if(_nbr_imult == 2) dothemagic_rhs_complex_nmult2();
            // else if(_nbr_imult == 3) dothemagic_rhs_complex_nmult3();
            else
                UP_CHECK1(false, "the number of imult = %d is not supported", _nbr_imult);
        }
    } else {
        UP_CHECK1(false, "type of solver %d not implemented", type);
    }

    hdf5_dump(_topo_hat[2], "sol_h", mydata);

    //-------------------------------------------------------------------------
    /** - go back to reals */
    //-------------------------------------------------------------------------
    for (int ip = 2; ip >= 0; ip--) {
        _plan_backward[ip]->execute_plan();
        // get if we are now complex
        if (_plan_forward[ip]->isr2c()) {
            _topo_hat[ip]->switch2real();
        }
        _switchtopo[ip]->execute(mydata, UP_BACKWARD);
    }

    //-------------------------------------------------------------------------
    /** - copy the solution in the field */
    //-------------------------------------------------------------------------
    for (int i2 = 0; i2 < topo->nloc(ax2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++) {
                // comnpute the index permutation
                const size_t id = localindex_ao(i0, i1, i2, topo);
                // put the data
                myfield[id] = mydata[id];
            }
        }
    }

    hdf5_dump(topo, "sol", myfield);
}

/**
 * @brief perform the convolution for real to real cases
 * 
 */
void FFTW_Solver::dothemagic_rhs_real() {
    BEGIN_FUNC

    UP_CHECK0(_topo_hat[2]->axis() == _topo_green[2]->axis(), "field and Green must have the same axis");

    opt_double_ptr       mydata  = _data;
    const opt_double_ptr mygreen = _green;

    const int ax0 = _topo_hat[2]->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    for (int i2 = 0; i2 < _topo_hat[2]->nloc(ax2); ++i2) {
        for (int i1 = 0; i1 < _topo_hat[2]->nloc(ax1); ++i1) {
            size_t id       = localindex_ao(0, i1, i2, _topo_hat[2]);
            size_t id_green = localindex_ao(_shiftgreen[ax0], i1 + _shiftgreen[ax1], i2 + _shiftgreen[ax2], _topo_green[2]);
            for (int i0 = 0; i0 < _topo_hat[2]->nloc(ax0); ++i0) {
                mydata[id] *= _normfact * mygreen[id_green];

                ++id;
                ++id_green;
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult0() {
    BEGIN_FUNC

    printf("doing the dothemagic_rhs_complex_nmult0\n");

    UP_CHECK0(_topo_hat[2]->axis() == _topo_green[2]->axis(), "field and Green must have the same axis");

    opt_double_ptr       mydata  = _data;
    const opt_double_ptr mygreen = _green;

    const int ax0 = _topo_hat[2]->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    for (int i2 = 0; i2 < _topo_hat[2]->nloc(ax2); ++i2) {
        for (int i1 = 0; i1 < _topo_hat[2]->nloc(ax1); ++i1) {
            size_t id       = localindex_ao(0, i1, i2, _topo_hat[2]);
            size_t id_green = localindex_ao(_shiftgreen[ax0], i1 + _shiftgreen[ax1], i2 + _shiftgreen[ax2], _topo_green[2]);

            for (int i0 = 0; i0 < _topo_hat[2]->nloc(ax0); ++i0) {
                const double a = mydata[id + 0];
                const double b = mydata[id + 1];
                const double c = mygreen[id_green + 0];
                const double d = mygreen[id_green + 1];

                // update the values
                mydata[id + 0] = _normfact * (a * c - b * d);
                mydata[id + 1] = _normfact * (a * d + b * c);

                id += 2;
                id_green += 2;
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult1() {
    BEGIN_FUNC
    UP_CHECK0(false, "not implemented yet");
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-1)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult2() {
    BEGIN_FUNC
    UP_CHECK0(false, "not implemented yet");
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult3() {
    BEGIN_FUNC
    UP_CHECK0(false, "not implemented yet");
}
