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
 * @brief Construct a new fftw solver::fftw solver object
 * 
 * @param topo the current #Topology object
 * @param mybc the boundary conditions asked on the solution
 */
FFTW_Solver::FFTW_Solver(const Topology* topo,const BoundaryType mybc[DIM][2])
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Store the field size */
    //-------------------------------------------------------------------------
    // for(int id=0; id<DIM; id++) _size_field[id] = topo->nglob(id);

    //-------------------------------------------------------------------------
    /** - For each dim, compute the plan and its type */
    //-------------------------------------------------------------------------
    double h[3] = {topo->h(0),topo->h(1),topo->h(2)};
    double L[3] = {topo->L(0),topo->L(1),topo->L(2)};

    for(int id=0; id<DIM; id++)
    {
        _plan_forward[id]   = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_FORWARD,false);
        _plan_backward[id]  = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_BACKWARD,false);
        _plan_green[id]     = new FFTW_plan_dim(id,h,L,mybc[id],FFTW_FORWARD,true);
    }

    _sort_plan(_plan_forward);
    _sort_plan(_plan_backward);
    _sort_plan(_plan_green);

    //-------------------------------------------------------------------------
    /** - Initialise the plans and get the sizes + dimorder */
    //-------------------------------------------------------------------------
    _init_plan(topo,_topo_hat,_reorder,_plan_forward,false);
    _init_plan(topo,NULL,NULL,_plan_backward,false);
    _init_plan(topo,_topo_green,_reorder_green,_plan_green,true);

    //-------------------------------------------------------------------------
    /** - Get the normalization factors #_normfact, #_volfact and the #_shiftgreen factor */
    //-------------------------------------------------------------------------
    _normfact  = 1.0;
    _volfact   = 1.0;
    for(int ip=0; ip<3;ip++){
        _normfact *= _plan_forward[ip]->normfact();
        _volfact  *= _plan_forward[ip]->volfact ();

        _shiftgreen[_plan_forward[ip]->dimID()] = _plan_forward[ip]->shiftgreen();
    }

    //-------------------------------------------------------------------------
    /** - Get the imult factor */
    //-------------------------------------------------------------------------
    _nbr_imult = 0;
    for(int ip=0; ip<3;ip++){
        if(_plan_forward[ip]->imult())  _nbr_imult ++;
        if(_plan_backward[ip]->imult()) _nbr_imult --;
        if(_plan_green[ip]->imult())    _nbr_imult ++;
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
void FFTW_Solver::setup(const SolverType mytype)
{
    //-------------------------------------------------------------------------
    /** - Store the solver type */
    //-------------------------------------------------------------------------
    _type = mytype;

    UP_CHECK1(mytype == UP_SRHS,"unsupported solver type %d\n",mytype);

    // //-------------------------------------------------------------------------
    // /** - Store some usefull factors in 'double' index calculus */
    // //-------------------------------------------------------------------------
    // // _dim_multfact is used to compute the tranposed index location = 
    // // _dim_multfact[0] * ix +  _dim_multfact[1] * iy + _dim_multfact[2] * iz
    // size_t acc = 1;
    // for(int id=0; id<3; id++){
    //     _dim_multfact[_dimorder[id]] = acc;
    //     acc *= _size_hat[id];
    // }
    // // if we are complex, we have to double the indexes to use them with doubles, except the first ones
    // if(_isComplex) for(int id=1; id<3; id++) _dim_multfact[_dimorder[id]] *=2;

    // if(_isComplex){
    //     _offset = _fieldstart[2]*_size_hat[1]*(_size_hat[0]*2)  + _fieldstart[1] * (_size_hat[0]*2) + _fieldstart[0];
    // }
    // else{
    //     _offset = _fieldstart[2]*_size_hat[1]*_size_hat[0]      + _fieldstart[1] * _size_hat[0]     + _fieldstart[0];
    // }

    //-------------------------------------------------------------------------
    /** - allocate the data */
    //-------------------------------------------------------------------------
    _allocate_data(_topo_hat,   &_data);
    _allocate_data(_topo_green, &_green);

    //-------------------------------------------------------------------------
    /** - allocate the plans */
    //-------------------------------------------------------------------------
    _allocate_plan(_topo_hat,_plan_forward ,_data);
    _allocate_plan(_topo_hat,_plan_backward,_data);

    //-------------------------------------------------------------------------
    /** - allocate the plan and comnpute the Green's function */
    //-------------------------------------------------------------------------
    _allocate_plan(_topo_green,_plan_green,_green);
    _compute_Green(_topo_green,_green,_plan_green);

    // attention on doit recalculer les topos des Green sur les tailles FULL
    // car on ne sait pas gagner comme on le fait pour les fields avec les 0

    //-------------------------------------------------------------------------
    /** - do the forward transfrom for green*/
    //-------------------------------------------------------------------------
    for(int ip=0; ip<3;ip++)
    {
        _reorder_green[ip]->execute(_green,FFTW_FORWARD);
        _plan_green[ip]->execute_plan();
    }

    //-------------------------------------------------------------------------
    /** - delete the plan */
    //-------------------------------------------------------------------------
    _delete_plan(_plan_green);
    _delete_reorder(_reorder_green);

}

/**
 * @brief Destroy the fftw solver
 * 
 */
FFTW_Solver::~FFTW_Solver(){
    BEGIN_FUNC
    // delete plans
    _delete_plan(_plan_forward);
    _delete_plan(_plan_backward);

    // delete datas
    if(_green != NULL) fftw_free(_green);
    if(_data  != NULL) fftw_free(_data);

    // delete reorder
    for(int id=0; id<3; id++){
        delete _reorder_green[id];
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
    }
}

void FFTW_Solver::_delete_reorder(Reorder_MPI *reorder[3])
{
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++)
    {
        delete reorder[ip];
    }
}
void FFTW_Solver::_delete_topology(Topology *topo[3])
{
    BEGIN_FUNC
    // deallocate the plans
    for (int ip = 0; ip < 3; ip++)
    {
        delete topo[ip];
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

void FFTW_Solver::_init_plan(const Topology *topo, Topology* topomap[3],Reorder_MPI* reorder[3], FFTW_plan_dim* planmap[3], bool isGreen)
{
    BEGIN_FUNC

    //-------------------------------------------------------------------------
    /** - Store the current topology */
    //-------------------------------------------------------------------------
    const Topology* current_topo = topo;
    const double h[3] = {current_topo->h(0),current_topo->h(1),current_topo->h(2)};

    //-------------------------------------------------------------------------
    /** - get the sizes to start with */
    //-------------------------------------------------------------------------
    int size_tmp[DIM];
    for(int id=0; id<DIM; id++) size_tmp[id] = topo->nglob(id);

    //-------------------------------------------------------------------------
    /** - create the plan and the topologies */
    //-------------------------------------------------------------------------
    bool isComplex = false;
    for(int ip=0; ip<3;ip++)
    {
        // initialize the plan
        planmap[ip]->init(size_tmp,&isComplex);
        // update the size
        planmap[ip]->get_outsize(size_tmp); // update the size for the next plans, keep order unchanged
        

        if (topomap != NULL && !isGreen)
        {
            // get the fastest rotating index
            int dimID = planmap[ip]->dimID(); // store the correspondance of the transposition

            // get the number of procs
            int nproc[3];
            int id1 = (dimID+1)%3;
            int id2 = (dimID+2)%3;

            _pencil_nproc(dimID,nproc,topo->comm_size());
            UP_CHECK4(nproc[0]*nproc[1]*nproc[2] == topo->comm_size,"the number of proc %d %d %d does not match the comm size %d",nproc[0],nproc[1],nproc[2],topo->comm_size);

            // create the new topology
            topomap[ip] = new Topology(dimID, size_tmp, nproc, h);

            if (reorder != NULL)
            {
                // get the fieldstart = the point where the old topo has to begin in the new
                int fieldstart[3] = {0};
                planmap[ip]->get_fieldstart(fieldstart);
                // create the reorderMPI to change topology
                reorder[ip] = new Reorder_MPI(1, current_topo, topomap[ip], fieldstart);
            }
        }
        // display it
        planmap[ip]->disp();
    }
    // if we are Green, we need to compute the topologies with the FULL size
    if(isGreen && topomap != NULL)
    {
        for(int ip=0; ip<3;ip++)
        {
            // get the fastest rotating index
            int dimID = planmap[ip]->dimID(); // store the correspondance of the transposition

            // get the number of procs
            int nproc[3];
            int id1 = (dimID+1)%3;
            int id2 = (dimID+2)%3;

            _pencil_nproc(dimID,nproc,topo->comm_size());
            UP_CHECK4(nproc[0]*nproc[1]*nproc[2] == topo->comm_size,"the number of proc %d %d %d does not match the comm size %d",nproc[0],nproc[1],nproc[2],topo->comm_size);

            // create the new topology with the size_tmp = the final size
            topomap[ip] = new Topology(dimID, size_tmp, nproc, h);

            if (reorder != NULL)
            {
                // get the fieldstart = the point where the old topo has to begin in the new
                int fieldstart[3] = {0};
                planmap[ip]->get_fieldstart(fieldstart);
                // create the reorderMPI to change topology
                reorder[ip] = new Reorder_MPI(1, current_topo, topomap[ip], fieldstart);
            }
            // display it
            planmap[ip]->disp();
        }
    }
}
void  FFTW_Solver::_allocate_plan(const Topology* const topo[3], FFTW_plan_dim* planmap[3], double* data)
{
    BEGIN_FUNC

    for(int ip=0; ip<3; ip++){
        int size_plan[3] = {topo[ip]->nloc(0),topo[ip]->nloc(1),topo[ip]->nloc(2)};
        planmap[ip]->allocate_plan(size_plan,data);
    }
    
}

/**
 * @brief 
 * 
 * @param topo_hat the topologies of the pencils
 * @param data 
 */
void FFTW_Solver::_allocate_data(const Topology* const topo_hat[3],double** data)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0((*data) == NULL,"Pointer has to be NULL for allocation");

    //-------------------------------------------------------------------------
    /** - Do the memory allocation */
    //-------------------------------------------------------------------------
    // the bigger size will be in the pencils
    size_t size_tot =1;
    for(int id=0; id<3; id++) size_tot = std::max(topo_hat[id]->locsize(),size_tot);

    INFOLOG2("Complex memory allocation, size = %ld\n",size_tot);
    (*data) =(double*) fftw_malloc(size_tot*sizeof(fftw_complex));

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
    UP_CHECK1(UP_ISALIGNED(*data),"FFTW alignement not compatible with UP_ALIGNMENT (=%d)",UP_ALIGNMENT);
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
void FFTW_Solver::_compute_Green(const Topology * const topo[3], double *green, FFTW_plan_dim *planmap[3])
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - Sanity checks */
    //-------------------------------------------------------------------------
    assert(green != NULL);

    //-------------------------------------------------------------------------
    /** - get the direction where we need to do spectral diff and count them */
    //-------------------------------------------------------------------------
    int  order[3];
    bool dospectral[3] = {false};

    double hfact[3];
    double kfact[3];
    int symstart[3];
    
    for(int ip=0; ip<3; ip++)
    {
        const int dimID   = planmap[ip]->dimID();

        dospectral[dimID] = planmap[ip]->dospectral();
        symstart[dimID]   = planmap[ip]->symstart();
        hfact[dimID]      = topo[ip]->h(dimID);

        if(dospectral[dimID]){
            hfact[dimID] = 0.0;
            kfact[dimID] = planmap[ip]->kfact();
        }
    }

    // count the number of spectral dimensions
    int nbr_spectral = 0;
    for(int id=0; id<3; id++) if(dospectral[id]) nbr_spectral++;

    //-------------------------------------------------------------------------
    /** - get the expression of Green in the full domain*/
    //-------------------------------------------------------------------------
    
    if(nbr_spectral == 0){
        INFOLOG(">> using Green function 3 dir unbounded\n");
        // if(DIM==2){
        //     if     (_greenorder == CHAT_2) Green_3D_3dirunbounded_0dirspectral_chat2(topo[0],hfact,green);
        //     else if(_greenorder == HEJ_2 ) Green_3D_3dirunbounded_0dirspectral_hej2 (topo[0],hfact,green,_greenalpha);
        //     else if(_greenorder == HEJ_4 ) Green_3D_3dirunbounded_0dirspectral_hej4 (topo[0],hfact,green,_greenalpha);
        //     else UP_CHECK2(false,"Green Function = %d  unknow for nbr_spectral = %d",_greenorder,nbr_spectral);
        // }
        if(DIM==3){
            if     (_greenorder == CHAT_2) Green_3D_3dirunbounded_0dirspectral_chat2(topo[0],hfact,green);
            else if(_greenorder == HEJ_2 ) Green_3D_3dirunbounded_0dirspectral_hej2 (topo[0],hfact,green,_greenalpha);
            else if(_greenorder == HEJ_4 ) Green_3D_3dirunbounded_0dirspectral_hej4 (topo[0],hfact,green,_greenalpha);
            else UP_CHECK2(false,"Green Function = %d  unknow for nbr_spectral = %d",_greenorder,nbr_spectral);
        }
    }
    else if(nbr_spectral == 1){
        INFOLOG(">> using Green function 2 dir unbounded - 1 dir spectral\n");
        // _compute_Green_2dirunbounded_1dirspectral
    } 
    else if(nbr_spectral == 2){
        INFOLOG(">> using Green function 1 dir unbounded - 2 dir spectral\n");
        // _compute_Green_1dirunbounded_2dirspectral
    }
    else if(nbr_spectral == 3){
        INFOLOG(">> using Green function 3 dir spectral\n");
        // _compute_Green_0dirunbounded_3dirspectral
    }

    //-------------------------------------------------------------------------
    /** - do the symmetry */
    //-------------------------------------------------------------------------

    // {int mysize[2] = {dsize_green[0],dsize_green[1]};
    // write_array(mysize,green,"green_ext");}

    // check if we have to symmetrize a direction

    int ax0 = topo[0]->axis();
    int ax1 = (ax0+1)%3;
    int ax2 = (ax0+2)%3;

    for(int i2=0; i2<topo[0]->nloc(ax2); i2++){
        for(int i1=0; i1<topo[0]->nloc(ax1); i1++){
            for(int i0=0; i0<topo[0]->nloc(ax0); i0++){

                const size_t id = i0 + topo[0]->nloc(ax0)*(i1 + topo[0]->nloc(ax1)*i2);
                // we have to take the symmetry around symstart: symstart - (iy - symstart) = 2 symstart - iy
                // to use the abs we have to go back to integers
                const int is0 = (symstart[ax0]==0 || i0 <= symstart[ax0]) ? i0 : abs(2*(int)symstart[ax0]-i0);
                const int is1 = (symstart[ax1]==0 || i1 <= symstart[ax1]) ? i1 : abs(2*(int)symstart[ax1]-i1);
                const int is2 = (symstart[ax2]==0 || i2 <= symstart[ax2]) ? i2 : abs(2*(int)symstart[ax2]-i2);

                const size_t id_sym = is0 + topo[0]->nloc(ax0)*(is1 + topo[0]->nloc(ax1)*is2);

                UP_CHECK3(id >= 0,"ID is not positive => id = %d, %d, %d",i0,i1,i2);
                UP_CHECK3(id < topo[0]->nloc(ax0)*topo[0]->nloc(ax1)*topo[0]->nloc(ax2),"ID is greater than the max size => id = %d, %d, %d",i0,i1,i2);

                green[id] = green[id_sym];
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - scale the Green data using #_volfact */
    //-------------------------------------------------------------------------
    for(int i2=0; i2<topo[0]->nloc(ax2); i2++){
        for(int i1=0; i1<topo[0]->nloc(ax1); i1++){
            for(int i0=0; i0<topo[0]->nloc(ax0); i0++){
                const size_t id = i0 + topo[0]->nloc(ax0)*(i1 + topo[0]->nloc(ax1)*i2);
                green[id] = green[id] * _volfact;
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
void FFTW_Solver::solve(const Topology* topo,double * field, double * rhs)
{
    BEGIN_FUNC
    //-------------------------------------------------------------------------
    /** - sanity checks */
    //-------------------------------------------------------------------------
    UP_CHECK0(field != NULL,"field is NULL");
    UP_CHECK0(rhs   != NULL,"rhs is NULL");
    UP_CHECK1(UP_ISALIGNED(field),"pointer no aligned to UP_ALIGNMENT (=%d)",UP_ALIGNMENT);
    UP_CHECK1(UP_ISALIGNED(rhs  ),"pointer no aligned to UP_ALIGNMENT (=%d)",UP_ALIGNMENT);
    

    opt_double_ptr myfield     = field;
    opt_double_ptr mydata      = (double*) _data;
    const opt_double_ptr myrhs = rhs;

    //-------------------------------------------------------------------------
    /** - clean the data memory */
    //-------------------------------------------------------------------------
    // get the data pointer to double style    
    size_t size_tot = topo->locsize();
    for(int id=0; id<3; id++) size_tot = std::max(_topo_hat[id]->locsize(),size_tot);
    std::memset(mydata,0,sizeof(double)*size_tot);

    //-------------------------------------------------------------------------
    /** - copy the rhs in the correct order */
    //-------------------------------------------------------------------------
    INFOLOG ("------------------------------------------\n");
    INFOLOG ("## memory information\n")
    INFOLOG4("- size field   = %d %d %d\n",_size_field[0],_size_field[1],_size_field[2]);
    INFOLOG4("- size hat     = %d %d %d\n",_size_hat[0],_size_hat[1],_size_hat[2]);
    INFOLOG4("- dim order    = %d %d %d\n",_dimorder[0],_dimorder[1],_dimorder[2]);
    INFOLOG4("- field start  = %d %d %d\n",_fieldstart[0],_fieldstart[1],_fieldstart[2]);
    INFOLOG4("- dim multfact = %d %d %d\n",_dim_multfact[0],_dim_multfact[1],_dim_multfact[2]);
    INFOLOG2("- offset       = %ld\n",_offset);
    INFOLOG ("------------------------------------------\n");

    int ax0 = topo->axis();
    int ax1 = (ax0+1)%3;
    int ax2 = (ax0+2)%3;
    for(int i2=0; i2<topo->nloc(ax2); i2++){
        for(int i1=0; i1<topo->nloc(ax1); i1++){
            for(int i0=0; i0<topo->nloc(ax0); i0++){
                // comnpute the index permutation
                const size_t id   = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
                // put the data
                mydata[id] = myrhs[id] ;
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - go to Fourier */
    //-------------------------------------------------------------------------
    for(int ip=0; ip<3; ip++)
    {
        _reorder[ip]->execute(mydata,FFTW_FORWARD);
        _plan_forward[ip]->execute_plan();
    }

    //-------------------------------------------------------------------------
    /** - Perform the magic */
    //-------------------------------------------------------------------------
    printf("doing the magic\n");
    // if(_type == UP_SRHS){
    //     if(!_isComplex){
    //         if     (_nbr_imult == 0) dothemagic_rhs_real();
    //         else UP_CHECK1(false,"the number of imult = %d is not supported\n",_nbr_imult);
    //     }
    //     else{
    //         if     (_nbr_imult == 0) dothemagic_rhs_complex_nmult0();
    //         else if(_nbr_imult == 1) dothemagic_rhs_complex_nmult1();
    //         else if(_nbr_imult == 2) dothemagic_rhs_complex_nmult2();
    //         else if(_nbr_imult == 3) dothemagic_rhs_complex_nmult3();
    //         else UP_CHECK1(false,"the number of imult = %d is not supported\n",_nbr_imult);
    //     }
    // }

    //-------------------------------------------------------------------------
    /** - go back to reals */
    //-------------------------------------------------------------------------
    for(int ip=2; ip>=0; ip--)
    {
        _plan_backward[ip]->execute_plan();
        _reorder[ip]->execute(mydata,FFTW_BACKWARD);
        
    }

    //-------------------------------------------------------------------------
    /** - copy the solution in the field */
    //-------------------------------------------------------------------------
    for(int i2=0; i2<topo->nloc(ax2); i2++){
        for(int i1=0; i1<topo->nloc(ax1); i1++){
            for(int i0=0; i0<topo->nloc(ax0); i0++){
                // comnpute the index permutation
                const size_t id   = i0 + topo->nloc(ax0) * (i1 + topo->nloc(ax1) * i2);
                // put the data
                myrhs[id]=mydata[id];
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
    opt_double_ptr mydata  = _data;
    const opt_double_ptr mygreen = _green;

    size_t id =0;
    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            int id_green = _shiftgreen[0] + _size_hat_green[0] * ( iy+_shiftgreen[1] + _size_hat_green[1] * (iz+_shiftgreen[2]));

            for(int ix=0; ix<_size_hat[0]; ix++){
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
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            size_t id       = _size_hat[0] * (iy + _size_hat[1] * iz);
            size_t id_green = _shiftgreen[0] + _size_hat_green[0] * (iy+_shiftgreen[1] + _size_hat_green[1] * (iz+_shiftgreen[2]));
            for(int ix=0; ix<_size_hat[0]; ix++){

                // const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                // const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                const double a = mydata[id][0];
                const double b = mydata[id][1];
                const double c = mygreen[id_green][0];
                const double d = mygreen[id_green][1];

                // update the values
                mydata[id][0] = _normfact * (a*c-b*d);
                mydata[id][1] = _normfact * (a*d+b*c);

                id++;
                id_green++;
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
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                const size_t id       = ix + _size_hat[0]       * (iy + _size_hat[1]       * iz);
                const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                mydata[id][0] = (-1.0) * _normfact * ( temp_imag );
                mydata[id][1] =          _normfact * ( temp_real );
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (-1)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult2()
{
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){

                const size_t id       = ix + _size_hat[0]       * (iy + _size_hat[1]       * iz);
                const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                mydata[id][0] = (-1.0) * _normfact * ( temp_real );
                mydata[id][1] = (-1.0) * _normfact * ( temp_imag );
            }
        }
    }
}

/**
 * @brief Do the convolution between complex data and complex Green's function and multiply by (i)
 * 
 */
void FFTW_Solver::dothemagic_rhs_complex_nmult3()
{
    opt_complex_ptr mydata = (fftw_complex*) _data;
    const opt_complex_ptr mygreen = (fftw_complex*) _green;

    for(int iz=0; iz<_size_hat[2]; iz++){
        for(int iy=0; iy<_size_hat[1]; iy++){
            for(int ix=0; ix<_size_hat[0]; ix++){   

                const size_t id       = ix + _size_hat[0]       * (iy + _size_hat[1]       * iz);
                const size_t id_green = ix + _size_hat_green[0] * (iy + _size_hat_green[1] * iz);

                const double temp_real = mydata[id][0]* mygreen[id_green][0] - mydata[id][1] * mygreen[id_green][1];
                const double temp_imag = mydata[id][1]* mygreen[id_green][0] + mydata[id][0] * mygreen[id_green][1];

                mydata[id][0] =          _normfact * ( temp_real );
                mydata[id][1] = (-1.0) * _normfact * ( temp_imag );
            }
        }
    }
}