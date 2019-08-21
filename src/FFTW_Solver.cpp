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
FFTW_Solver::FFTW_Solver(const Topology *topo, const BoundaryType mybc[DIM][2], const double h[3], const double L[3])
{
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

    for (int id = 0; id < DIM; id++)
    {
        _plan_forward[id] = new FFTW_plan_dim(id, h, L, mybc[id], FFTW_FORWARD, false);
        _plan_backward[id] = new FFTW_plan_dim(id, h, L, mybc[id], FFTW_BACKWARD, false);
        _plan_green[id] = new FFTW_plan_dim(id, h, L, mybc[id], FFTW_FORWARD, true);
    }

    _sort_plan(_plan_forward);
    _sort_plan(_plan_backward);
    _sort_plan(_plan_green);

    //-------------------------------------------------------------------------
    /** - Initialise the plans and get the sizes */
    //-------------------------------------------------------------------------
    _init_plan(topo, _topo_hat, _reorder, _plan_forward, false);
    _init_plan(topo, NULL, NULL, _plan_backward, false);
    _init_plan(topo, _topo_green, _reorder_green, _plan_green, true);

    //-------------------------------------------------------------------------
    /** - Get the factors #_normfact, #_volfact, #_shiftgreen and #_nbr_imult */
    //-------------------------------------------------------------------------
    _normfact = 1.0;
    _volfact = 1.0;
    _nbr_imult = 0;
    for (int ip = 0; ip < 3; ip++)
    {
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
void FFTW_Solver::setup()
{
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

    // // attention on doit recalculer les topos des Green sur les tailles FULL
    // // car on ne sait pas gagner comme on le fait pour les fields avec les 0

    // //-------------------------------------------------------------------------
    // /** - do the forward transfrom for green*/
    // //-------------------------------------------------------------------------
    // for(int ip=0; ip<3;ip++)
    // {
    //     // go to the topology for the plan
    //     _reorder_green[ip]->execute(_green,FFTW_FORWARD);

    //     // do a symmetry if needed
    //     cmptGreenSymmetry(_topo_green,_green);

    //     // execute the plan
    //     _plan_green[ip]->execute_plan();
    // }

    //-------------------------------------------------------------------------
    /** - delete the useless data for Green */
    //-------------------------------------------------------------------------
    _delete_plan(_plan_green);
    _delete_reorder(_reorder_green);
}

/**
 * @brief Destroy the fftw solver
 * 
 */
FFTW_Solver::~FFTW_Solver()
{
    BEGIN_FUNC
    // delete plans
    _delete_plan(_plan_forward);
    _delete_plan(_plan_backward);

    // delete datas
    if (_green != NULL)
        fftw_free(_green);
    if (_data != NULL)
        fftw_free(_data);

    // delete reorder
    for (int id = 0; id < 3; id++)
    {
        if (_reorder[id] != NULL)
            delete _reorder[id];
    }

    //cleanup
    fftw_cleanup();
}
/**
 * @brief delete the FFTW_plan_dim stored in planmap
 * 
 * @param planmap 
 */
void FFTW_Solver::_delete_plan(FFTW_plan_dim *planmap[3])
{
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++)
    {
        delete planmap[ip];
        planmap[ip] = NULL;
    }
}

void FFTW_Solver::_delete_reorder(Reorder_MPI *reorder[3])
{
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++)
    {
        delete reorder[ip];
        reorder[ip] = NULL;
    }
}
void FFTW_Solver::_delete_topology(Topology *topo[3])
{
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++)
    {
        delete topo[ip];
        topo[ip] = NULL;
    }
}

void FFTW_Solver::_sort_plan(FFTW_plan_dim *plan[3])
{
    int priority[3];

    for (int id = 0; id < 3; id++)
        priority[id] = plan[id]->type();

    // do the sort by hand...
    if (priority[0] > priority[1])
    {
        FFTW_plan_dim *temp_plan = plan[1];
        plan[1] = plan[0];
        plan[0] = temp_plan;
    }
    // we are sure to have the first to sorted
    // if the last one is smaller than the second one, we change, if not, it is also bigger than the first one
    if (priority[1] > priority[2])
    {
        FFTW_plan_dim *temp_plan = plan[2];
        plan[2] = plan[1];
        plan[1] = temp_plan;
        if (priority[0] > priority[1])
        {
            FFTW_plan_dim *temp_plan = plan[1];
            plan[1] = plan[0];
            plan[0] = temp_plan;
        }
    }
}

/**
 * @brief Initialize a set of 3 plans by doing a dry run through the plans
 * 
 * @param topo the starting topology
 * @param topomap the topology array to go through each dim ( may be NULL) it corresponds to the topology AFTER the plan
 * @param reorder the reordering array to switch between topologies (may be NULL)
 * @param planmap the plan that will be created
 * @param isGreen indicates if the plans are for Green
 */
void FFTW_Solver::_init_plan(const Topology *topo, Topology *topomap[3], Reorder_MPI *reorder[3], FFTW_plan_dim *planmap[3], bool isGreen)
{
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
    int nproc[3];
    for (int ip = 0; ip < 3; ip++)
    {
        // initialize the plan
        planmap[ip]->init(size_tmp, isComplex);
        // update the size_tmp variable and get the complex information
        planmap[ip]->get_outsize(size_tmp);
        // we store the temp topologies for non-Green plans
        if (!isGreen && topomap != NULL && reorder != NULL)
        {
            // get the fastest rotating index
            int dimID = planmap[ip]->dimID(); // store the correspondance of the transposition
            // get the proc repartition
            _pencil_nproc(dimID, nproc, topo->comm_size());
            // create the new topology
            topomap[ip] = new Topology(dimID, size_tmp, nproc, isComplex);
            // get the fieldstart = the point where the old topo has to begin in the new
            int fieldstart[3] = {0};
            planmap[ip]->get_fieldstart(fieldstart);
            // create the reorderMPI to change topology
            reorder[ip] = new Reorder_MPI(current_topo, topomap[ip], fieldstart);
            // virtually execute the plan...
            // if the data is now complex, update and change the topo
            if (planmap[ip]->switch2Complex())
            {
                isComplex = true;
                topomap[ip]->switch2complex();
            }
            // update the current topo
            current_topo = topomap[ip];
        }
    }

    //-------------------------------------------------------------------------
    /** - For Green we need to compute the topologies using the full size of the domain  */
    //-------------------------------------------------------------------------
    isComplex = false;
    if (isGreen && topomap != NULL && reorder != NULL)
    {
        for (int ip = 0; ip < 3; ip++)
        {
            // get the fastest rotating index
            int dimID = planmap[ip]->dimID(); // store the correspondance of the transposition
            // get the proc repartition
            _pencil_nproc(dimID, nproc, topo->comm_size());
            // create the new topology
            topomap[ip] = new Topology(dimID, size_tmp, nproc, isComplex);
            // get the fieldstart = the point where the old topo has to begin in the new
            int fieldstart[3] = {0};
            planmap[ip]->get_fieldstart(fieldstart);
            // create the reorderMPI to change topology
            reorder[ip] = new Reorder_MPI(current_topo, topomap[ip], fieldstart);
            // virtually execute the plan...
            // if the data is now complex, update and change the topo
            if (planmap[ip]->switch2Complex())
            {
                isComplex = true;
                topomap[ip]->switch2complex();
            }
            // update the current topo
            current_topo = topomap[ip];
        }
    }
    // if(isGreen && topomap != NULL)
    // {

    //     for(int ip=0; ip<3;ip++)
    //     {
    //         // // initialize the plan
    //         // planmap[ip]->init(size_tmp,isComplex);
    //         // // update the size
    //         // // planmap[ip]->get_outsize(size_tmp); // update the size for the next plans, keep order unchanged
    //         // planmap[ip]->get_isNowComplex(&isComplex);
    //         // get the fastest rotating index
    //         int dimID = planmap[ip]->dimID(); // store the correspondance of the transposition

    //         /**
    //          * @todo
    //          *
    //          * the plans in planmap have been created on the wrong sizes... is it an issue?
    //          */

    //         // get the number of procs
    //         int nproc[3];
    //         _pencil_nproc(dimID,nproc,topo->comm_size());
    //         UP_CHECK4(nproc[0]*nproc[1]*nproc[2] == topo->comm_size(),"the number of proc %d %d %d does not match the comm size %d",nproc[0],nproc[1],nproc[2],topo->comm_size());

    //         // create the new topology with the size_tmp = the final size
    //         topomap[ip] = new Topology(dimID, size_tmp, nproc);

    //         if (reorder != NULL)
    //         {
    //             // get the fieldstart = the point where the old topo has to begin in the new
    //             int fieldstart[3] = {0};
    //             planmap[ip]->get_fieldstart(fieldstart);
    //             // create the reorderMPI to change topology
    //             printf("\ncreating reorder for Green with fieldstart = %d %d %d\n\n\n",fieldstart[0],fieldstart[1],fieldstart[2]);
    //             reorder[ip] = new Reorder_MPI(1, current_topo, topomap[ip], fieldstart);
    //         }

    //         current_topo = topomap[ip];
    //     }
    // }
}
void FFTW_Solver::_allocate_plan(const Topology *const topo[3], FFTW_plan_dim *planmap[3], double *data)
{
    BEGIN_FUNC

    for (int ip = 0; ip < 3; ip++)
    {
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
void FFTW_Solver::_allocate_data(const Topology *const topo_hat[3], double **data)
{
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

    // if(_isComplex){
    //     INFOLOG2("Complex memory allocation, size = %ld\n",size_tot);
    //     (*data) =(double*) fftw_malloc(size_tot*sizeof(fftw_complex));
    // }
    // else{
    //     INFOLOG2("Real memory allocation, size = %ld\n",size_tot);
    //     (*data) =(double*) fftw_malloc(size_tot*sizeof(double));
    // }

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
void FFTW_Solver::_cmptGreenFunction(Topology * topo[3], double *green, FFTW_plan_dim *planmap[3])
{
    BEGIN_FUNC

    //-------------------------------------------------------------------------
    /** - get the direction where we need to do spectral diff and count them */
    //-------------------------------------------------------------------------
    int order[3];
    bool dospectral[3] = {false};

    double hfact[3];
    double kfact[3];
    int symstart[3];

    for (int ip = 0; ip < 3; ip++)
    {
        const int dimID = planmap[ip]->dimID();
        // get usefull datas
        dospectral[dimID] = planmap[ip]->dospectral();
        symstart[dimID] = planmap[ip]->symstart();
        hfact[dimID] = _hgrid[dimID];
        kfact[dimID] = 0.0;

        if (dospectral[dimID])
        {
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
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //-------------------------------------------------------------------------
    /** - set the Green function to real if needed */
    //-------------------------------------------------------------------------
    bool firstTopoIsComplex = topo[0]->isComplex();
    if(firstTopoIsComplex) topo[0]->switch2real();

    //-------------------------------------------------------------------------
    /** - get the expression of Green in the full domain*/
    //-------------------------------------------------------------------------
    if (nbr_spectral == 0)
    {
        INFOLOG(">> using Green function 3 dir unbounded\n");
        // if(DIM==2){
        //     if     (_greenorder == CHAT_2) Green_3D_3dirunbounded_0dirspectral_chat2(topo[0],hfact,green);
        //     else if(_greenorder == HEJ_2 ) Green_3D_3dirunbounded_0dirspectral_hej2 (topo[0],hfact,green,_greenalpha);
        //     else if(_greenorder == HEJ_4 ) Green_3D_3dirunbounded_0dirspectral_hej4 (topo[0],hfact,green,_greenalpha);
        //     else UP_CHECK2(false,"Green Function = %d  unknow for nbr_spectral = %d",_greenorder,nbr_spectral);
        // }
        if (DIM == 3)
        {
            // if (_greenorder == CHAT_2)
            //     Green_3D_3dirunbounded_0dirspectral_chat2(topo[0], hfact, symstart, green);
            // else
            if (_greenorder == HEJ_2)
                Green_3D_3dirunbounded_0dirspectral_hej2(topo[0], hfact, symstart, green, _greenalpha);
            // else if (_greenorder == HEJ_4)
            //     // Green_3D_3dirunbounded_0dirspectral_hej4(topo[0], hfact, green, _greenalpha);
            else
                UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _greenorder, nbr_spectral);
        }
    }
    else if (nbr_spectral == 1)
    {
        INFOLOG(">> using Green function 2 dir unbounded - 1 dir spectral\n");
        // _compute_Green_2dirunbounded_1dirspectral
        UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _greenorder, nbr_spectral);
    }
    else if (nbr_spectral == 2)
    {
        INFOLOG(">> using Green function 1 dir unbounded - 2 dir spectral\n");
        // _compute_Green_1dirunbounded_2dirspectral
        UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _greenorder, nbr_spectral);
    }
    else if (nbr_spectral == 3)
    {
        INFOLOG(">> using Green function 3 dir spectral\n");
        // _compute_Green_0dirunbounded_3dirspectral
        UP_CHECK2(false, "Green Function = %d  unknow for nbr_spectral = %d", _greenorder, nbr_spectral);
    }

    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // printf("[%d] topo is of size %d %d %d in axis %d\n",rank,topo[0]->nglob(0),topo[0]->nglob(1),topo[0]->nglob(2),topo[0]->axis());
    // printf("[%d] topo is of size %d %d %d in axis %d\n",rank,topo[0]->nloc(0),topo[0]->nloc(1),topo[0]->nloc(2),topo[0]->axis());

    
    xmf_write(topo[0], "green", "data");
    hdf5_write(topo[0], "green", "data",green);

    //-------------------------------------------------------------------------
    /** - scale the Green data using #_volfact */
    //-------------------------------------------------------------------------
    _scaleGreenFunction(topo[0], green);

    //-------------------------------------------------------------------------
    /** - compute a symmetry and do the forward transform*/
    //-------------------------------------------------------------------------
    bool isComplex = false;
    for (int ip = 0; ip < 3; ip++)
    {
        // go to the topology for the plan
        if (ip > 0)
        {
            _reorder_green[ip]->execute(_green, FFTW_FORWARD);
        }

        if (ip == 0)
        {
            xmf_write(topo[0],  "green_t0", "data");
            hdf5_write(topo[0], "green_t0", "data", green);
        }
        else if (ip == 1)
        {
            xmf_write(topo[1],  "green_t1", "data");
            hdf5_write(topo[1], "green_t1", "data", green);
        }
        else if (ip == 2)
        {
            xmf_write(topo[2],  "green_t2", "data");
            hdf5_write(topo[2], "green_t2", "data", green);
        }

        // execute the plan
        _plan_green[ip]->execute_plan();

        if (ip == 0)
        {
            xmf_write(topo[0],  "green_t0_h", "data");
            hdf5_write(topo[0], "green_t0_h", "data", green);
        }
        else if (ip == 1)
        {
            xmf_write(topo[1],  "green_t1_h", "data");
            hdf5_write(topo[1], "green_t1_h", "data", green);
        }
        else if (ip == 2)
        {
            xmf_write(topo[2],  "green_t2_h", "data");
            hdf5_write(topo[2], "green_t2_h", "data", green);
        }
    }

    // //-------------------------------------------------------------------------
    // /** - do the symmetry */
    // //-------------------------------------------------------------------------

    // // {int mysize[2] = {dsize_green[0],dsize_green[1]};
    // // write_array(mysize,green,"green_ext");}

    // // check if we have to symmetrize a direction

    // int ax0 = topo[0]->axis();
    // int ax1 = (ax0+1)%3;
    // int ax2 = (ax0+2)%3;

    // for(int i2=0; i2<topo[0]->nloc(ax2); i2++){
    //     for(int i1=0; i1<topo[0]->nloc(ax1); i1++){
    //         for(int i0=0; i0<topo[0]->nloc(ax0); i0++){

    //             const size_t id = i0 + topo[0]->nloc(ax0)*(i1 + topo[0]->nloc(ax1)*i2);
    //             // we have to take the symmetry around symstart: symstart - (iy - symstart) = 2 symstart - iy
    //             // to use the abs we have to go back to integers
    //             const int is0 = (symstart[ax0]==0 || i0 <= symstart[ax0]) ? i0 : abs(2*(int)symstart[ax0]-i0);
    //             const int is1 = (symstart[ax1]==0 || i1 <= symstart[ax1]) ? i1 : abs(2*(int)symstart[ax1]-i1);
    //             const int is2 = (symstart[ax2]==0 || i2 <= symstart[ax2]) ? i2 : abs(2*(int)symstart[ax2]-i2);

    //             const size_t id_sym = is0 + topo[0]->nloc(ax0)*(is1 + topo[0]->nloc(ax1)*is2);

    //             UP_CHECK3(id >= 0,"ID is not positive => id = %d, %d, %d",i0,i1,i2);
    //             UP_CHECK3(id < topo[0]->nloc(ax0)*topo[0]->nloc(ax1)*topo[0]->nloc(ax2),"ID is greater than the max size => id = %d, %d, %d",i0,i1,i2);

    //             green[id] = green[id_sym];
    //         }
    //     }
    // }
}

/**
 * @brief compute the Green's function symmetry along the fastest rotating index
 * 
 * We take the symmetry around sym_idx: sym_idx - (i0 - sym_idx) = 2 sym_idx - i0
 * 
 * @warning The padding starts at sym_idx+1 due to the vertex-centered data of the Green's function
 * 
 * @param topo the topology we are currently on
 * @param sym_idx the index of symmetry: data(sym_idx+1) is padded
 * @param data the data storage
 * @param isComplex boolean to indicate if the data is
 */
// void FFTW_Solver::_cmptGreenSymmetry(const Topology *topo, const int sym_idx, double *data, const bool isComplex)
// {
//     BEGIN_FUNC

//     // if no symmetry is required, we stop
//     if (sym_idx == 0)
//         return;

//     // the symmetry is done along the fastest rotating index
//     const int ax0 = topo->axis();
//     const int ax1 = (ax0 + 1) % 3;
//     const int ax2 = (ax0 + 2) % 3;

//     printf("\n");
//     printf("doing the symmetry along axis %d around id %d\n",topo->axis(),sym_idx);

//     if (!isComplex)
//     {
//         for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
//         {
//             for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
//             {
//                 // this direction is continuous in memory hence no idstart is required
//                 for (int i0 = sym_idx + 1; i0 < topo->nloc(ax0); i0++)
//                 {

//                     const size_t id = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
//                     const size_t id_sym = abs(2 * (int)sym_idx - i0) + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);

//                     UP_CHECK3(id >= 0, "ID is not positive => id = %d, %d, %d", i0, i1, i2);
//                     UP_CHECK3(id < topo->nloc(ax0) * topo->nloc(ax1) * topo->nloc(ax2), "ID is greater than the max size => id = %d, %d, %d", i0, i1, i2);

//                     data[id] = data[id_sym];
//                 }
//             }
//         }
//     }
//     else
//     {
//         for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
//         {
//             for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
//             {
//                 for (int i0 = sym_idx + 1; i0 < topo->nloc(ax0); i0 += 2)
//                 {

//                     const size_t id = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
//                     const size_t id_sym = abs(2 * (int)sym_idx - i0) + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);

//                     UP_CHECK3(id >= 0, "ID is not positive => id = %d, %d, %d", i0, i1, i2);
//                     UP_CHECK3(id < topo->nloc(ax0) * topo->nloc(ax1) * topo->nloc(ax2), "ID is greater than the max size => id = %d, %d, %d", i0, i1, i2);

//                     data[id] = data[id_sym];
//                     data[id + 1] = data[id_sym + 1];
//                 }
//             }
//         }
//     }
// }
void FFTW_Solver::_scaleGreenFunction(const Topology *topo, double *data)
{
    // the symmetry is done along the fastest rotating index
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
    {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
        {
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++)
            {
                const size_t id = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
                data[id] = data[id] * _volfact;
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
void FFTW_Solver::solve(const Topology *topo, double *field, double *rhs, const SolverType type)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0(field != NULL, "field is NULL");
    UP_CHECK0(rhs != NULL, "rhs is NULL");
    UP_CHECK1(UP_ISALIGNED(field), "pointer no aligned to UP_ALIGNMENT (=%d)", UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(rhs), "pointer no aligned to UP_ALIGNMENT (=%d)", UP_ALIGNMENT);

    opt_double_ptr myfield = field;
    opt_double_ptr mydata = (double *)_data;
    const opt_double_ptr myrhs = rhs;

    //-------------------------------------------------------------------------
    /** - clean the data memory */
    //-------------------------------------------------------------------------
    // get the data pointer to double style
    size_t size_tot = topo->locmemsize();
    for (int id = 0; id < 3; id++)
        size_tot = std::max(_topo_hat[id]->locmemsize(), size_tot);
    std::memset(mydata, 0, sizeof(double) * size_tot);

    //-------------------------------------------------------------------------
    /** - copy the rhs in the correct order */
    //-------------------------------------------------------------------------
    INFOLOG("------------------------------------------\n");
    INFOLOG("## memory information\n")
    INFOLOG4("- size field   = %d %d %d\n", _size_field[0], _size_field[1], _size_field[2]);
    INFOLOG4("- size hat     = %d %d %d\n", _size_hat[0], _size_hat[1], _size_hat[2]);
    INFOLOG4("- dim order    = %d %d %d\n", _dimorder[0], _dimorder[1], _dimorder[2]);
    INFOLOG4("- field start  = %d %d %d\n", _fieldstart[0], _fieldstart[1], _fieldstart[2]);
    INFOLOG4("- dim multfact = %d %d %d\n", _dim_multfact[0], _dim_multfact[1], _dim_multfact[2]);
    INFOLOG2("- offset       = %ld\n", _offset);
    INFOLOG("------------------------------------------\n");

    int ax0 = topo->axis();
    int ax1 = (ax0 + 1) % 3;
    int ax2 = (ax0 + 2) % 3;
    if (topo->nf() == 1)
    {
        for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
        {
            for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
            {
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++)
                {
                    // comnpute the index permutation
                    const size_t id = localindex_ao(i0, i1, i2, topo);
                    // put the data
                    mydata[id] = myrhs[id];
                }
            }
        }
    }
    else if (topo->nf() == 2)
    {
        for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
        {
            for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
            {
                for (int i0 = 0; i0 < topo->nloc(ax0); i0++)
                {
                    // comnpute the index permutation
                    const size_t id = localindex_ao(i0, i1, i2, topo);
                    // put the data
                    mydata[id + 0] = myrhs[id + 0];
                    mydata[id + 1] = myrhs[id + 1];
                }
            }
        }
    }
    else
    {
        UP_CHECK0(false, "size of Topological nf not supported");
    }

    //-------------------------------------------------------------------------
    /** - go to Fourier */
    //-------------------------------------------------------------------------
    bool isComplex = false;
    for (int ip = 0; ip < 3; ip++)
    {
        // go to the correct topo
        _reorder[ip]->execute(mydata, FFTW_FORWARD);
        // run the FFT
        _plan_forward[ip]->execute_plan();
        // get if we are now complex
        _plan_forward[ip]->get_isNowComplex(&isComplex);
    }

    //-------------------------------------------------------------------------
    /** - Perform the magic */
    //-------------------------------------------------------------------------
    printf("doing the magic\n");
    if (type == UP_SRHS)
    {
        if (!isComplex)
        {
            if (_nbr_imult == 0)
                dothemagic_rhs_real();
            else
                UP_CHECK1(false, "the number of imult = %d is not supported", _nbr_imult);
        }
        else
        {
            if (_nbr_imult == 0)
                dothemagic_rhs_complex_nmult0();
            // else if(_nbr_imult == 1) dothemagic_rhs_complex_nmult1();
            // else if(_nbr_imult == 2) dothemagic_rhs_complex_nmult2();
            // else if(_nbr_imult == 3) dothemagic_rhs_complex_nmult3();
            else
                UP_CHECK1(false, "the number of imult = %d is not supported", _nbr_imult);
        }
    }
    else
    {
        UP_CHECK1(false, "type of solver %d not implemented", type);
    }

    //-------------------------------------------------------------------------
    /** - go back to reals */
    //-------------------------------------------------------------------------
    for (int ip = 2; ip >= 0; ip--)
    {
        _plan_backward[ip]->execute_plan();
        _reorder[ip]->execute(mydata, FFTW_BACKWARD);
    }

    //-------------------------------------------------------------------------
    /** - copy the solution in the field */
    //-------------------------------------------------------------------------
    for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
    {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
        {
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++)
            {
                // comnpute the index permutation
                const size_t id = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
                // put the data
                myrhs[id] = mydata[id];
            }
        }
    }
}

/**
 * @brief perform the convolution for real to real cases
 * 
 */
void FFTW_Solver::dothemagic_rhs_real()
{
    BEGIN_FUNC

    UP_CHECK0(_topo_hat[2]->axis() == _topo_green[2]->axis(),"field and Green must have the same axis");

    opt_double_ptr mydata = _data;
    const opt_double_ptr mygreen = _green;

    size_t id = 0;
    const int ax0 = _topo_hat[2]->axis();
    const int ax1 = (ax0+1)%3;
    const int ax2 = (ax0+2)%3;

    for (int i2 = 0; i2 < _topo_hat[2]->nloc(ax2); ++i2)
    {
        for (int i1 = 0; i1 < _topo_hat[2]->nloc(ax1); ++i1)
        {
            size_t id       = localindex_ao(0,i1,i2,_topo_hat[2]);
            size_t id_green = localindex_ao(_shiftgreen[ax0],i1 + _shiftgreen[ax1],i2 + _shiftgreen[ax2],_topo_green[2]);
            for (int i0 = 0; i0 < _topo_hat[2]->nloc(ax0); ++i0)
            {
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
void FFTW_Solver::dothemagic_rhs_complex_nmult0()
{
    BEGIN_FUNC

    UP_CHECK0(_topo_hat[2]->axis() == _topo_green[2]->axis(),"field and Green must have the same axis");

    opt_double_ptr mydata = _data;
    const opt_double_ptr mygreen = _green;

    const int ax0 = _topo_hat[2]->axis();
    const int ax1 = (ax0+1)%3;
    const int ax2 = (ax0+2)%3;

    for (int i2 = 0; i2 < _topo_hat[2]->nloc(ax2); ++i2)
    {
        for (int i1 = 0; i1 < _topo_hat[2]->nloc(ax1); ++i1)
        {
            size_t id       = localindex_ao(0,i1,i2,_topo_hat[2]);
            size_t id_green = localindex_ao(_shiftgreen[ax0],i1 + _shiftgreen[ax1],i2 + _shiftgreen[ax2],_topo_green[2]);

            for (int i0 = 0; i0 < _topo_hat[2]->nloc(ax0); ++i0)
            {
                const double a = mydata[id+0];
                const double b = mydata[id+1];
                const double c = mygreen[id_green+0];
                const double d = mygreen[id_green+1];

                // update the values
                mydata[id+0] = _normfact * (a * c - b * d);
                mydata[id+1] = _normfact * (a * d + b * c);

                id+= 2;
                id_green+= 2;
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult1()
{
    BEGIN_FUNC
    UP_CHECK0(false,"not implemented yet");
    // opt_complex_ptr mydata = (fftw_complex *)_data;
    // const opt_complex_ptr mygreen = (fftw_complex *)_green;

    // for (int iz = 0; iz < _size_hat[2]; iz++)
    // {
    //     for (int iy = 0; iy < _size_hat[1]; iy++)
    //     {
    //         for (int ix = 0; ix < _size_hat[0]; ix++)
    //         {

    //             const size_t id = ix + _size_hat[0] * (iy + _size_hat[1] * iz);
    //             const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

    //             const double temp_real = mydata[id][0] * mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
    //             const double temp_imag = mydata[id][1] * mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

    //             mydata[id][0] = (-1.0) * _normfact * (temp_imag);
    //             mydata[id][1] = _normfact * (temp_real);
    //         }
    //     }
    // }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-1)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult2()
{

    BEGIN_FUNC
    UP_CHECK0(false,"not implemented yet");

    // opt_complex_ptr mydata = (fftw_complex *)_data;
    // const opt_complex_ptr mygreen = (fftw_complex *)_green;

    // for (int iz = 0; iz < _size_hat[2]; iz++)
    // {
    //     for (int iy = 0; iy < _size_hat[1]; iy++)
    //     {
    //         for (int ix = 0; ix < _size_hat[0]; ix++)
    //         {

    //             const size_t id = ix + _size_hat[0] * (iy + _size_hat[1] * iz);
    //             const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

    //             const double temp_real = mydata[id][0] * mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
    //             const double temp_imag = mydata[id][1] * mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

    //             mydata[id][0] = (-1.0) * _normfact * (temp_real);
    //             mydata[id][1] = (-1.0) * _normfact * (temp_imag);
    //         }
    //     }
    // }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult3()
{

    BEGIN_FUNC
    UP_CHECK0(false,"not implemented yet");

    // opt_complex_ptr mydata = (fftw_complex *)_data;
    // const opt_complex_ptr mygreen = (fftw_complex *)_green;

    // for (int iz = 0; iz < _size_hat[2]; iz++)
    // {
    //     for (int iy = 0; iy < _size_hat[1]; iy++)
    //     {
    //         for (int ix = 0; ix < _size_hat[0]; ix++)
    //         {

    //             const size_t id = ix + _size_hat[0] * (iy + _size_hat[1] * iz);
    //             const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

    //             const double temp_real = mydata[id][0] * mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
    //             const double temp_imag = mydata[id][1] * mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

    //             mydata[id][0] = _normfact * (temp_real);
    //             mydata[id][1] = (-1.0) * _normfact * (temp_imag);
    //         }
    //     }
    // }
}