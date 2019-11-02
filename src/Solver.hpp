/**
 * @file Solver.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

#ifndef FFTW_SOLVER_HPP
#define FFTW_SOLVER_HPP

#include <cstring>
#include <map>
#include "FFTW_plan_dim.hpp"
#include "defines.hpp"
#include "green_functions_3d.hpp"
#include "hdf5_io.hpp"

#include "SwitchTopo.hpp"
#include "SwitchTopo_a2a.hpp"
#include "SwitchTopo_nb.hpp"

#include "Profiler.hpp"
#include "omp.h"

#ifdef HAVE_METIS
#include "metis.h"
#endif

using namespace std;

/**
 * @brief The Poisson solver
 * 
 * A collection of 3 FFTW_plan_dim for the forward and backward FFT transform of
 * data, plus the transformed required for the Green's function to solve the Poisson equation.
 * The tranforms are done in-place in each direction successively. Between each transform, data
 * are remapped (i.e. transposed) in order to have the memory aligned with the direction of the
 * transform. This is done using SwitchTopo which changes the layout of data between 2 topos.
 * 
 * @warning
 * The memory alignement follows the rules explained on the mainpage.
 * Yet, a transposition of the data is required to perfom the transfroms in the correct order.
 * For an element (ix,iy,iz) its tranposed location is computed as\code{.cpp}
 * const size_t id_transposed = iz*_dim_multfact[2] + iy*_dim_multfact[1] + ix*_dim_multfact[0] + _offset;
 * \endcode
 *  
 * 
 */
class Solver {
    // the memory allocation is assumed to be data[iz][iy][ix]
    // so the fastest running index is n[0] then n[1] then n[2]
    // even is the dimension is 2, we allocate arrays of dimension 3

   protected:
    int _fftwalignment = 0; /**< @brief alignement assumed by the FFTW Solver  */
    int _orderdiff     = 0; /**< @brief the order of derivative (spectral = 0)  */
    int _nbr_imult     = 0; /**< @brief the number of time we have applied a DST transform */
    int _nbr_spectral  = 0; /** @brief the number of spectral directions involved     */

    double  _normfact = 1.0;   /**< @brief normalization factor so that the forward/backward FFT gives output = input */
    double  _volfact  = 1.0;   /**< @brief volume factor due to the convolution computation */
    double  _hgrid[3] = {0.0}; /**< @brief grid spacing in the tranposed directions */
    double* _data     = NULL;  /**< @brief data pointer to the transposed memory */

    /**
     * @name Forward and backward 
     * transforms related objects
     */
    /**@{ */
    FFTW_plan_dim* _plan_forward[3];  /**< @brief map containing the plans for the forward fft transforms */
    FFTW_plan_dim* _plan_backward[3]; /**< @brief map containing the plans for the backward fft transforms */
    Topology*      _topo_phys     = NULL;
    Topology*      _topo_hat[3]   = {NULL, NULL, NULL}; /**< @brief map containing the topologies (i.e. data memory layout) corresponding to each transform */
    SwitchTopo*    _switchtopo[3] = {NULL, NULL, NULL}; /**< @brief switcher of topologies for the forward transform (phys->topo[0], topo[0]->topo[1], topo[1]->topo[2]).*/
    opt_double_ptr _sendBuf       = NULL;               /**<@brief The send buffer for _switchtopo */
    opt_double_ptr _recvBuf       = NULL;               /**<@brief The recv buffer for _switchtopo */
    /**@} */

    /**
     * @name Green's function (and corresponding forward transform) related vars and objects
     * 
     */
    /**@{ */
    // int       _shiftgreen[3]   = {0, 0, 0}; /**< @brief the shift in the Green's function which chose to take the flip-flop mode or not */
    double    _alphaGreen      = 2.0;       /**< @brief regularization parameter for HEJ_* Green's functions */
    double*   _green           = NULL;      /**< @brief data pointer to the transposed memory for Green */
    GreenType _typeGreen       = CHAT_2;    /**< @brief the type of Green's function */

    FFTW_plan_dim* _plan_green[3]; /**< @brief map containing the plan for the Green's function */
    Topology*   _topo_green[3]       = {NULL, NULL, NULL}; /**< @brief list of topos dedicated to Green's function */
    SwitchTopo* _switchtopo_green[3] = {NULL, NULL, NULL}; /**< @brief switcher of topos for the Green's forward transform*/
    /**@} */

    // time the solve
    Profiler* _prof = NULL;

   protected:
    /**
     * @name Data management
     * 
     * @{
     */
    void _allocate_data(const Topology *const topo[3], const Topology *topo_phys, double **data);
    void _delete_switchtopos(SwitchTopo* switchtopo[3]);
    void _delete_topologies(Topology* topo[3]);
    /**@}  */

    /**
     * @name Plan management
     * 
     * @{
     */
    void _sort_plans(FFTW_plan_dim* plan[3]);
    void _init_plansAndTopos(const Topology* topo, Topology* topomap[3], SwitchTopo* switchtopo[3], FFTW_plan_dim* planmap[3], bool isGreen);
    void _allocate_plans(const Topology* const topo[3], FFTW_plan_dim* planmap[3], double* data);
    void _delete_plans(FFTW_plan_dim* planmap[3]);
    /**@} */

    /**
     * @name SwitchTopo management
     * 
     * @{
     */
    void _allocate_switchTopo(const int ntopo, SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    void _deallocate_switchTopo(SwitchTopo** switchtopo, opt_double_ptr* send_buff, opt_double_ptr* recv_buff);
    /**@} */

    /**
     * @name Do the magic
     * 
     * @{
     */
    void dothemagic_rhs_real(double *data);
    void dothemagic_rhs_complex_nmult0(double *data);
    void dothemagic_rhs_complex_nmult1(double *data);
    void dothemagic_rhs_complex_nmult2(double *data);
    void dothemagic_rhs_complex_nmult3(double *data);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void _cmptGreenFunction(Topology* topo[3], double* green, FFTW_plan_dim* planmap[3]);
    void _cmptGreenSymmetry(const Topology* topo, const int sym_idx, double* data, const bool isComplex);
    void _scaleGreenFunction(const Topology* topo, double* data, bool killModeZero);
    void _finalizeGreenFunction(Topology* topo_field, double* green, const Topology* topo, FFTW_plan_dim* planmap[3]);
    /**@} */

   public:
    Solver(Topology* topo, const BoundaryType mybc[3][2], const double h[3], const double L[3],Profiler* prof);
    ~Solver();

    double* setup(const bool changeTopoComm);
    void set_OrderDiff(const int order) { _orderdiff = order; }
    const Topology* get_innerTopo_physical() ;
    const Topology* get_innerTopo_spectral() ;

    /**
     * @brief Get the total allocated size of the pointer data (returned by setup)
     * 
     * @return size_t 
     */
    size_t get_allocSize() {
        size_t size_tot = 1;
        for (int id = 0; id < 3; id++) {
            size_tot = std::max(_topo_hat[id]->memsize(), size_tot);
        }
        return size_tot;
    };

    /**
     * @brief Get the spectral information to compute the modes k in full spectral space
     * 
     * @param kfact  multiply the index by this factor to obtain the wave number (1/2/3 corresponds to x/y/z )
     * @param koffset  add this to the index to obtain the wave number (1/2/3 corresponds to x/y/z )
     */
    void get_spectralInfo(double kfact[3], double koffset[3], double symstart[3]) {
        for (int ip = 0; ip < 3; ip++) {
            const int dimID = _plan_forward[ip]->dimID();
            kfact[dimID]    = _plan_forward[ip]->kfact();
            koffset[dimID]  = _plan_forward[ip]->koffset();
            symstart[dimID] = _plan_forward[ip]->symstart();
        }
    }

    /**
     * @name Solver use 
     * 
     * @{
     */
    void solve(double* field, double* rhs, const SolverType type);
    /**@} */

    /**
     * @name Solver use (advanced)
     * 
     * @{
     */
    void do_copy(const Topology *topo, double *data, const int sign );
    void do_FFT(double *data, const int sign);
    void do_mult(double *data, const SolverType type);
    /**@} */

    /**
     * @name Green's function
     * 
     * @{
     */
    void set_GreenType(const GreenType type) { _typeGreen = type; }
    void set_alpha(const double alpha) { _alphaGreen = alpha; }
    /**@} */
};

/**
 * @brief compute the pencil layout given the pencil direction
 * 
 * The pencil layout is computed so as to obtain pencils with an aspect
 * ratio close to 1, i.e. the same number points per proc in the the 2 other directions than id.
 * 
 * @param id the pencil direction
 * @param nproc the number of proc in each direction
 * @param comm_size the total communicator size
 * @param nglob the domain size in each direction
 */
static inline void pencil_nproc(const int id, int nproc[3], const int comm_size, const int nglob[3]) {
    int id1 = (id + 1) % 3;
    int id2 = (id + 2) % 3;

    nproc[id] = 1;

    double       n1       = 1;
    double       n2       = (double) comm_size;
    //invert indexes so that id1 is the dimension where nglob is the smallest
    if( nglob[id1] > nglob[id2]){
        const int tmp = id2;
        id2 = id1;
        id1 = tmp;
    }
    double       np1      = (double) nglob[id1];
    double       np2      = (double) nglob[id2]/ comm_size;
    const double npsquare = sqrt((double)(nglob[id1] * nglob[id2]) / comm_size);  //target number of points per dimension

    //keep on deviding as long as ncurr/2>nsquare
    //we want to leave n1=1, and we do not want to reach n2=1
    while ( (np1 > npsquare) && std::floor(n2*.5) == n2*.5) {
        n1  *= 2.0;
        np1 *= 0.5;
        n2  *= 0.5;
        np2 *= 2.0;
    }
    nproc[id1] = (int)n1;
    nproc[id2] = (int)n2;

    FLUPS_INFO("my proc repartition is %d %d %d",nproc[0],nproc[1],nproc[2]);
    if(nproc[0] * nproc[1] * nproc[2] != comm_size){
        FLUPS_ERROR("the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size, LOCATION);
    }
    if(comm_size>8 && (n1==1||n2==1)){
        FLUPS_WARNING("A slab decomposition was used instead of a pencil decomposition in direction %d. This may increase communication time.",id, LOCATION);
        //Loss of performance may originate in slab decompositions, as an actual All2All communication is required, whereas with the pencils,
        // we manage to do All2All communications in subcoms of size sqrt(comm_size).
        //We could prevent this to happen by doing something like:
        // if(n2==1){
        //     n2*=2;
        //     n1*=0.5;
        // }
    }
}

/**
 * @brief compute the pencil layout given the pencil direction, compatible with another pencil decoposition given as a hint
 * 
 * @param id the pencil direction
 * @param nproc the number of proc in each direction
 * @param comm_size the total communicator size
 * @param id_hint the axis of the pencils in another decomposition, which we want this decomposition to be compatible with
 * @param nproc_hint the number of procs in the other decomposition we want to be compatible with
 */
static inline void pencil_nproc_hint(const int id, int nproc[3], const int comm_size, const int id_hint, const int nproc_hint[3]) {
    // get the id shared between the hint topo
    int sharedID = 0;
    for (int i = 0; i < 3; i++) {
        if (i != id && i != id_hint) {
            sharedID = i;
            break;
        }
    }
    nproc[id]       = 1;
    nproc[sharedID] = nproc_hint[sharedID];
    nproc[id_hint]  = comm_size / nproc[sharedID];

    FLUPS_INFO("My proc repartition in this topo is %d %d %d",nproc[0],nproc[1],nproc[2]);
    FLUPS_CHECK(nproc[0] * nproc[1] * nproc[2] == comm_size, "the number of proc %d %d %d does not match the comm size %d", nproc[0], nproc[1], nproc[2], comm_size, LOCATION);

    if(comm_size>8 && (nproc[sharedID]==1||nproc[id_hint]==1)){
        FLUPS_WARNING("A slab decomposition was used instead of a pencil decomposition in direction %d. This may increase communication time.",id, LOCATION);
    }
}

/**
 * @brief reorder the MPI-ranks using metis
 * 
 * @warning this functions assume an evenly distributed amount of procs on the nodes
 * 
 * @param comm 
 * @param sources 
 * @param sourcesW 
 * @param dests 
 * @param destsW 
 * @param n_nodes 
 * @param order 
 */
static void reorder_metis(MPI_Comm comm, int *sources, int *sourcesW, int *dests, int *destsW, int *order) {
    int comm_size;
    int comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

#ifdef HAVE_METIS

    //-------------------------------------------------------------------------
    /** - get the total number of nodes */
    //-------------------------------------------------------------------------
    // create a group where everybody can create a shared memory region
    MPI_Comm nodecomm;
    MPI_Info mpinfo;
    MPI_Info_create(&mpinfo);
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, comm_rank, mpinfo, &nodecomm);
    // we store the comm size
    int local_nodesize;
    MPI_Comm_size(nodecomm, &local_nodesize);

// #define PART_OF_EQUAL_SIZE
#ifdef PART_OF_EQUAL_SIZE
    //_______ OPTION 1 with gcd (suboptimal)________
    // gather on each proc the gcd
    int *vec_nodesize = (int  *)flups_malloc(sizeof(int) * comm_size);
    MPI_Allgather(&local_nodesize, 1, MPI_INT, vec_nodesize, 1, MPI_INT, comm);
    // get the Greatest Common Divider among every process
    int nodesize = comm_size;
    for (int ip = 0; ip < comm_size; ip++) {
        nodesize = gcd(nodesize, vec_nodesize[ip]);
    }
    // store the number of nodes
    int n_nodes = comm_size / nodesize;
    double* tpwgts = NULL;
#else
    //_______ OPTION 2 with various size partitions________
    // gather on proc 1 the number of proc per node
    int *vec_nodesize = (int *)flups_malloc(sizeof(int) * comm_size);
    MPI_Allgather(&local_nodesize, 1, MPI_INT, vec_nodesize, 1, MPI_INT, comm);
    
    // count the number of partitions we'll need:
    int n_nodes = 0;
    int id = 0;
    while( id < comm_size){
        id += vec_nodesize[id];
        n_nodes++;
    }
    
#ifdef DEV_SIMULATE_GRAPHCOMM
    //CHEATING: imposing that there will be 2 groups (there needs to be at least 4 procs)
    n_nodes = 2;
    for (int ip = 0; ip<comm_size; ip++ ){
        vec_nodesize[ip]=comm_size/3;
    }
    vec_nodesize[0]=comm_size-vec_nodesize[1];
#endif

    real_t* tpwgts = (real_t*) flups_malloc(sizeof(real_t)*n_nodes);
    // deduce the size of each partition:
    id = 0;
    for (int ip = 0; ip<n_nodes; ip++ ){
        tpwgts[ip] = ((real_t)vec_nodesize[id])/((real_t) comm_size);
        id += vec_nodesize[id];
    }
    //______________________________________________
#endif

    // free stuffs
    flups_free(vec_nodesize);
    MPI_Comm_free(&nodecomm);

    //-------------------------------------------------------------------------
    /** - get the neighbour list and the associated weights */
    //-------------------------------------------------------------------------
    // we count the number of neighbours
    int n_neighbours = 0;
    // if we have either one way in or one way out, we have a neighbour
    for (int i = 0; i < comm_size; ++i) {
        if ((sourcesW[i] + destsW[i]) > 0 && i != comm_rank) n_neighbours++;
    }
    // allocate the number of neighbours and their weights
    int *neighbours = (int *)flups_malloc(sizeof(int) * n_neighbours);
    int *weights    = (int *)flups_malloc(sizeof(int) * n_neighbours);
    n_neighbours    = 0;
    for (int i = 0; i < comm_size; ++i) {
        if (sourcesW[i] + destsW[i] > 0 && i != comm_rank) {
            neighbours[n_neighbours] = i;
            weights[n_neighbours]    = sourcesW[i] + destsW[i];
            n_neighbours++;
        }
    }

    //-------------------------------------------------------------------------
    /** - build the graph on proc 0 and ask for partioning
     * The graph structure follows metis rules:
     * the edges (= id of the destination of the edges) starting from proc k are located
     * from adj[xadj[k]] to adj[xadj[k+1]-1]
     * Same structure is used for the weights with the ajdw
     * */
    //-------------------------------------------------------------------------
    if (comm_rank == 0) {
        int *xadj = (int *)flups_malloc((comm_size + 1) * sizeof(int));
        int *nadj = (int *)flups_malloc((comm_size) * sizeof(int));

        // get the number of neighbours from everybody
        MPI_Gather(&n_neighbours, 1, MPI_INT, nadj, 1, MPI_INT, 0, comm);
        // get the starting indexes of the neighbour description for everybody
        xadj[0] = 0;
        for (int i = 0; i < comm_size; ++i) {
            xadj[i + 1] = xadj[i] + nadj[i];
        }

        // allocate the adjency list + weights and fill it with the neighbour list from everybody
        int *adj  = (int *)flups_malloc(xadj[comm_size] * sizeof(int));
        int *adjw = (int *)flups_malloc(xadj[comm_size] * sizeof(int));
        MPI_Gatherv(neighbours, n_neighbours, MPI_INT, adj, nadj, xadj, MPI_INT, 0, comm);
        MPI_Gatherv(weights, n_neighbours, MPI_INT, adjw, nadj, xadj, MPI_INT, 0, comm);
#ifdef PROF
        {
            //writing graph to file, CSR format
            string filename = "prof/graph.csr";
            FILE* file      = fopen(filename.c_str(), "w+");
            if(file==NULL){FLUPS_ERROR("Could not create file in ./prof. Did you create the folder?",LOCATION);}
            for(int i=0; i<=comm_size; i++){
                fprintf(file, "%d ",xadj[i]);
            }
            fprintf(file,"\n");
            for(int i=0; i<xadj[comm_size]; i++){
                fprintf(file, "%d (%d), ",adj[i],adjw[i]);
            }
            fprintf(file,"\n");
            fclose(file);

            //writing graph to file, per node
            filename = "prof/graph.txt";
            file     = fopen(filename.c_str(), "w+");
            for(int i=0; i<comm_size; i++){
                fprintf(file, "%d: ",i);
                for(int j = xadj[i]; j<xadj[i+1]; j++){
                        fprintf(file, "%d (%d), ",adj[j],adjw[j]);
                }
                fprintf(file,"\n");
            }
            fclose(file);
        }
#endif

        //prepare vall to metis
        int  ncon = 1;  // the number of balancing constraints
        real_t tol = 1.0001; //tolerance on the constraint 
        int  objval;
        int *part = (int *)flups_malloc(comm_size * sizeof(int));
        int *rids = (int *)flups_malloc(n_nodes * sizeof(int));
        std::memset(rids,0,sizeof(int)*n_nodes);

        // ask of the partitioning. call metis several times in case the tolerance on the partition size is not exactly respected 
        int max_iter = 10;
        if(n_nodes==1){
            max_iter = 0;
            FLUPS_WARNING("METIS: you asked only 1 node. I can't do the partitaioning.",LOCATION);
        }
        int iter;
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);

        //METIS options: ncuts and niter seems to have the most effect
        options[METIS_OPTION_SEED] = 1;
        options[METIS_OPTION_NCUTS] = 50;
        options[METIS_OPTION_NITER] = 50;
        options[METIS_OPTION_UFACTOR] = 1.;
        // options[METIS_OPTION_IPTYPE] = 
        //     METIS_IPTYPE_GROW,
        //     METIS_IPTYPE_RANDOM,
        //     METIS_IPTYPE_EDGE,
        //     METIS_IPTYPE_NODE,
        //     METIS_IPTYPE_METISRB
        // options[METIS_OPTION_RTYPE] = 
        //     METIS_RTYPE_FM,
        //     METIS_RTYPE_GREEDY,
        //     METIS_RTYPE_SEP2SIDED,
        //     METIS_RTYPE_SEP1SIDED

        for(iter = 0; iter<max_iter; iter++){
            FLUPS_INFO("METIS: graph partitioning attempt %d",iter+1);
            METIS_PartGraphRecursive(&comm_size, &ncon, xadj, adj, NULL, NULL, adjw, &n_nodes, tpwgts, &tol, options, &objval, part);
            tol=((tol-1.)/2.)+1.;
            options[METIS_OPTION_NCUTS] +=10;
            options[METIS_OPTION_NITER] +=10;
            options[METIS_OPTION_SEED] += 1; 
            // options[METIS_OPTION_UFACTOR] /= 1.;
            
            // compute how many proc in each group, resulting from metis partitioning
            for (int i = 0; i < comm_size; ++i) {
                rids[part[i]]++;
            }
            // check that we did respect the constraint on the size of the partitions
            bool succeed = (rids[0] == (int)(tpwgts[0] * comm_size));
            FLUPS_INFO("METIS:   part %d: size %d (should be %d)",0, rids[0], (int)(tpwgts[0] * comm_size));
            for (int ip = 1; ip < n_nodes; ++ip) {
                succeed &= (rids[ip] == (int)(tpwgts[ip] * comm_size));
                FLUPS_INFO("METIS:   part %d: size %d (should be %d)",ip, rids[ip], (int)(tpwgts[ip] * comm_size));
                rids[ip] += rids[ip-1]; //switch to cumulative numbering
            }
            for (int ip = n_nodes-1; ip > 0; --ip) {
                rids[ip] = rids[ip-1]; //offset by 1
            }
            rids[0] = 0;
            if(!succeed){
                FLUPS_INFO("METIS:   attempt failed.");
            }else{
                // assign the rank value and redistribute
                for (int i = 0; i < comm_size; ++i) {
                    order[i] = rids[part[i]]++ ;
                }
                break;
            }
        }
        // check that we did not reach max_iter
        if(iter>=max_iter){
            FLUPS_WARNING("Failed to find a graph partitioning with the current allocation. I will not change the rank orderegin in the graph_comm!",LOCATION);
            for (int i = 0; i < comm_size; ++i) {
                order[i] = i;
            }
        }

        // result of the partitioning
    #ifdef PART_OF_EQUAL_SIZE   
        FLUPS_INFO("I have partitioned the graph in %d chunks of size %d\n",n_nodes,comm_size/n_nodes);
    #else            
        FLUPS_INFO("I have partitioned the graph in %d chunks.",n_nodes);
    #endif
#ifdef PROF
        //writing graph to file, CSR format
        string filename = "prof/partitions.txt";
        FILE* file      = fopen(filename.c_str(), "w+");
    #ifdef PART_OF_EQUAL_SIZE   
        fprintf(file,"%d partitions of size %d\n",n_nodes,comm_size/n_nodes);
    #else
        fprintf(file,"%d partitions of size:\n",n_nodes);
        for(int i=0; i<n_nodes; i++){
            fprintf(file, "part %d with %d elems\n",i,(int)(comm_size*tpwgts[i]));
        }
    #endif
        fclose(file);
#endif        

        flups_free(xadj);
        flups_free(nadj);
        flups_free(adj);
        flups_free(adjw);

        flups_free(part);
        flups_free(rids);
    } else {
        MPI_Gather(&n_neighbours, 1, MPI_INT, NULL, 1, MPI_INT, 0, comm);
        MPI_Gatherv(neighbours, n_neighbours, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm);
        MPI_Gatherv(weights, n_neighbours, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, comm);
    }
    flups_free(neighbours);
    flups_free(weights);
#ifndef PART_OF_EQUAL_SIZE
    flups_free(tpwgts);
#endif

    //-------------------------------------------------------------------------
    /** - give the rank info to everybody */
    //-------------------------------------------------------------------------
    MPI_Bcast(order, comm_size, MPI_INT, 0, comm);
#ifdef PROF        
    if (comm_rank == 3) {
        //writing reordering to file
        string filename = "prof/order.txt";
        FILE* file      = fopen(filename.c_str(), "w+");
        for (int i = 0; i < comm_size; ++i) {
            FLUPS_INFO("METIS ORDER %i : %i \n", i, order[i]);
            printf("%i : %i \n", i, order[i]);
            fprintf(file,"%i : %i \n", i, order[i]);
        }
        fclose(file);
    }
#endif
#else
    for (int i = 0; i < comm_size; ++i) {
        order[i] = i;
    }
#endif
}


#endif
