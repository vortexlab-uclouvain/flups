/**
 * @file flups.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * 
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

#include "defines.hpp"
#include "h3lpr/macros.hpp"
#include "h3lpr/profiler.hpp"
#include "toolsinterface.hpp"
#include "Topology.hpp"
#include "Solver.hpp"



extern "C" {

void * flups_malloc(size_t size){
    return m_calloc(size);
}

void flups_free(void* data){
    m_free(data);
}


/**
 * @brief writes the file murphy.info used for tracking of the results, bookkeeping etc
 */
void flups_info(int argc, char** argv) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::string filename(argv[0]);
        // remove the "./" that might be in it
        filename.erase(std::remove(filename.begin(), filename.end(), '.'), filename.end());
        filename.erase(std::remove(filename.begin(), filename.end(), '/'), filename.end());
        filename += ".info";
        FILE* file = fopen(filename.c_str(), "w+");
        fprintf(file, "FLUPS \n");
        fprintf(file, "- commit: %s\n", FLUPS_GIT_COMMIT);
        fprintf(file, "- defines:\n");
        fprintf(file, "\tFLUPS_FFTW_FLAG = %d\n", FLUPS_FFTW_FLAG);
        fprintf(file, "\tFLUPS_ALIGNMENT = %d\n", FLUPS_ALIGNMENT);
#if (FLUPS_MPI_AGGRESSIVE)
        fprintf(file, "\tFLUPS_MPI_AGGRESSIVE ? yes\n");
#else
        fprintf(file, "\tFLUPS_MPI_AGGRESSIVE ? no\n");
#endif
#ifndef NDEBUG
        fprintf(file, "\tNDEBUG ? no\n");
#else
        fprintf(file, "\tNDEBUG ? yes\n");
#endif
#if (FLUPS_NEW_BALANCE)
        fprintf(file, "\tFLUPS_NEW_BALANCE ? yes\n");
#else
        fprintf(file, "\tFLUPS_NEW_BALANCE ? no\n");
#endif
#ifdef FLUPS_WISDOM_PATH
        fprintf(file, "\tFLUPS_WISDOM_PATH = %s\n", FLUPS_WISDOM_PATH);
#else
        fprintf(file, "\tFLUPS_WISDOM_PATH = none\n");
#endif
#if (FLUPS_OLD_MPI)
        fprintf(file, "\tMPI_40 ? no\n");
#else
        fprintf(file, "\tMPI_40 ? yes\n");
#endif
#ifdef COMM_NONBLOCK
        fprintf(file, "\tFLUPS_MPI_BATCH_SEND = %d\n", FLUPS_MPI_BATCH_SEND);
#endif
        fprintf(file, "- argument list:\n");
        for (int i = 1; i < argc; ++i) {
            fprintf(file, "\t%s\n", argv[i]);
        }
        fclose(file);
    }
    FLUPS_INFO("-------------------------------------------------------------------");
    FLUPS_INFO("FLUPS - (c) UCLOUVAIN, MIT");
    FLUPS_INFO("commit = %s", FLUPS_GIT_COMMIT);
    FLUPS_INFO("-------------------------------------------------------------------");
}

//***********************************************************************
// * TOPOLOGIES
// **********************************************************************/
Topology* flups_topo_new(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm){
    Topology* t = new Topology(axis, lda, nglob, nproc, isComplex, axproc, alignment, comm);
    return t;
}

void flups_topo_free(const Topology* t) {
    delete t;
}

bool flups_topo_get_isComplex(const Topology* t) {
    return t->isComplex();
}

int flups_topo_get_axis(const Topology* t) {
    return t->axis();
}

int flups_topo_get_nglob(const Topology* t, const int dim) {
    return t->nglob(dim);
}

int flups_topo_get_nloc(const Topology* t, const int dim) {
    return t->nloc(dim);
}

int flups_topo_get_nmem(const Topology* t, const int dim) {
    return t->nmem(dim);
}

int flups_topo_get_nproc(const Topology* t, const int dim) {
    return t->nproc(dim);
}

void flups_topo_get_istartGlob(const Topology* t, int istart[3]) {
    t->get_istart_glob(istart);
}

size_t flups_topo_get_locsize(const Topology* t) {
    return (size_t)t->locsize();
}

size_t flups_topo_get_memsize(const Topology* t) {
    return (size_t)t->memsize();
}

int flups_topo_cmpt_rank_fromid(const Topology* t, const int global_id, const int id){
    return t->cmpt_rank_fromid(global_id, id);
}

int flups_topo_cmpt_start_id_from_rank(const Topology* t, const int rank_id, const int id){
    return t->cmpt_start_id_from_rank(rank_id, id);
}

MPI_Comm flups_topo_get_comm(Topology* t){
    return t->get_comm();
}

void flups_topo_ranksplit(const Topology* t, const int rank, int rankd[3]) {
    int axproc[3] = {t->axproc(0),t->axproc(1),t->axproc(2)};
    int nproc[3]  = {t->nproc(0),t->nproc(1),t->nproc(2)};
    ranksplit(rank, axproc, nproc, t->get_comm(), rankd);
}

int flups_topo_rankindex(const Topology *topo, const int rankd[3]) {
    return rankindex(rankd, topo);
}

//***********************************************************************
//*  SOLVER
//********************************************************************* */

// get a new solver
Solver* flups_init(Topology* t, BoundaryType* bc[3][2], const double h[3], const double L[3], DiffType orderDiff, const CenterType center_type[3]) {
    Solver* s = new Solver(t, bc, h, L, orderDiff, center_type, NULL);
    return s;
}
Solver* flups_init_timed(Topology* t, BoundaryType* bc[3][2], const double h[3], const double L[3], const DiffType orderDiff, const CenterType center_type[3], H3LPR::Profiler* prof) {
    Solver* s = new Solver(t, bc, h, L, orderDiff, center_type, prof);
    return s;
}

// destroy the solver
void flups_cleanup(Solver* s){
    delete s;
}

// setup the solver
void flups_set_greenType(Solver* s, const GreenType type){
    s->set_GreenType(type);
}

double* flups_setup(Solver* s,const bool changeComm){
    return s->setup(changeComm);
}

// solve
void flups_solve(Solver* s, double* field, double* rhs, const SolverType type) {
    s->solve(field, rhs, type);
}


// -- ADVANCED FEATURES --

size_t flups_get_allocSize(Solver* s){
    return(size_t) s->get_allocSize();
}

void flups_get_spectralInfo(Solver* s, double kfact[3], double koffset[3], double symstart[3]){
    s->get_spectralInfo(kfact,koffset,symstart);
}

void flups_set_alpha(Solver* s, const double alpha){
    s->set_alpha(alpha);   
}

const Topology* flups_get_innerTopo_physical(Solver* s){
    return s->get_innerTopo_physical();
}

const Topology* flups_get_innerTopo_spectral(Solver* s){
    return s->get_innerTopo_spectral();
}

void flups_do_copy(Solver* s, const Topology* topo, double* data, const int sign){
    s->do_copy(topo,data,sign);
}

void flups_do_FFT(Solver* s, double* data, const int sign){
    s->do_FFT(data,sign);
}

void flups_do_mult(Solver* s, double* data,const SolverType type){
    s->do_mult(data,type);
}

int flups_hint_proc_repartition(const int lda, const double h[3], const double L[3], BoundaryType* bc[3][2], const CenterType center_type[3]){
    return hint_proc_repartition(lda, h, L, bc, center_type);
}

void flups_switchtopo_info(Solver* s){
    s->get_switchtopo_info();
}

//**********************************************************************
//  PROFILER - TIMERS
//**********************************************************************

H3LPR::Profiler* flups_profiler_new() {
    H3LPR::Profiler* p = new H3LPR::Profiler();
    // return reinterpret_cast<void*>(p);
    return p;
}

H3LPR::Profiler* flups_profiler_new_n(const char name[]){
    H3LPR::Profiler* p = new H3LPR::Profiler(name);
    // return reinterpret_cast<void*>(p);
    return p; 
}

void flups_profiler_free(H3LPR::Profiler* p) {
    // delete reinterpret_cast<H3LPR::Profiler*>(p);
    delete p;
}

void flups_profiler_disp(H3LPR::Profiler* p) {
    // m_profDisp(reinterpret_cast<H3LPR::Profiler*>(p));
    m_profDisp(p);
}


//**********************************************************************
//  HDF5
//**********************************************************************

void flups_hdf5_dump(const Topology *topo, const char filename[], const double *data){
    const std::string fn(filename);
    hdf5_dump(topo,fn, data);
}

void flups_print_data(const Topology *topo, double* data){
    FLUPS_print_data(topo, data);    
}

}
