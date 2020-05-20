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
#include "Topology.hpp"
#include "Solver.hpp"
#include "Profiler.hpp"


extern "C" {

void * flups_malloc(size_t size){
    return flups_mem_malloc(size);
}

void flups_free(void* data){
    flups_mem_free(data);
}

//***********************************************************************
// * TOPOLOGIES
// **********************************************************************/
FLUPS_Topology* flups_topo_new(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm){
    Topology* t = new Topology(axis, lda, nglob, nproc, isComplex, axproc, alignment, comm);
    return t;
}

void flups_topo_free(const FLUPS_Topology* t) {
    delete t;
}

bool flups_topo_get_isComplex(const FLUPS_Topology* t) {
    return t->isComplex();
}

int flups_topo_get_axis(const FLUPS_Topology* t) {
    return t->axis();
}

int flups_topo_get_nglob(const FLUPS_Topology* t, const int dim) {
    return t->nglob(dim);
}

int flups_topo_get_nloc(const FLUPS_Topology* t, const int dim) {
    return t->nloc(dim);
}

int flups_topo_get_nmem(const FLUPS_Topology* t, const int dim) {
    return t->nmem(dim);
}

int flups_topo_get_nproc(const FLUPS_Topology* t, const int dim) {
    return t->nproc(dim);
}

void flups_topo_get_istartGlob(const FLUPS_Topology* t, int istart[3]) {
    t->get_istart_glob(istart);
}

size_t flups_topo_get_locsize(const FLUPS_Topology* t) {
    return (size_t)t->locsize();
}

size_t flups_topo_get_memsize(const FLUPS_Topology* t) {
    return (size_t)t->memsize();
}

MPI_Comm flups_topo_get_comm(FLUPS_Topology* t){
    return t->get_comm();
}

//***********************************************************************
//*  SOLVER
//********************************************************************* */

// get a new solver
FLUPS_Solver* flups_init(FLUPS_Topology* t, FLUPS_BoundaryType* bc[3][2], const double h[3], const double L[3], const int orderDiff) {
    Solver* s = new Solver(t, bc, h, L, orderDiff, NULL);
    return s;
}
FLUPS_Solver* flups_init_timed(FLUPS_Topology* t, FLUPS_BoundaryType* bc[3][2], const double h[3], const double L[3], const int orderDiff, Profiler* prof) {
#ifndef PROF
    Solver* s = new Solver(t, bc, h, L, orderDiff, NULL);
#else
    Solver* s = new Solver(t, bc, h, L, orderDiff, prof);
#endif
    return s;
}

// destroy the solver
void flups_cleanup(FLUPS_Solver* s){
    delete s;
}

// setup the solver
void flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type){
    s->set_GreenType(type);
}

double* flups_setup(FLUPS_Solver* s,const bool changeComm){
    return s->setup(changeComm);
}

// solve
void flups_solve(FLUPS_Solver* s, double* field, double* rhs, const FLUPS_SolverType type) {
    s->solve(field, rhs, type);
}


// -- ADVANCED FEATURES --

size_t flups_get_allocSize(FLUPS_Solver* s){
    return(size_t) s->get_allocSize();
}

void flups_get_spectralInfo(FLUPS_Solver* s, double kfact[3], double koffset[3], double symstart[3]){
    s->get_spectralInfo(kfact,koffset,symstart);
}

void flups_set_alpha(FLUPS_Solver* s, const double alpha){
    s->set_alpha(alpha);   
}

// void flups_set_OrderDiff(FLUPS_Solver* s, const int order){
//     s->set_OrderDiff(order);
// }

const FLUPS_Topology* flups_get_innerTopo_physical(FLUPS_Solver* s){
    return s->get_innerTopo_physical();
}

const FLUPS_Topology* flups_get_innerTopo_spectral(FLUPS_Solver* s){
    return s->get_innerTopo_spectral();
}

void flups_do_copy(FLUPS_Solver* s, const FLUPS_Topology* topo, double* data, const int sign){
    s->do_copy(topo,data,sign);
}

void flups_do_FFT(FLUPS_Solver* s, double* data, const int sign){
    s->do_FFT(data,sign);
}

void flups_do_mult(FLUPS_Solver* s, double* data,const FLUPS_SolverType type){
    s->do_mult(data,type);
}

//**********************************************************************
//  PROFILER - TIMERS
//**********************************************************************

FLUPS_Profiler* flups_profiler_new() {
    Profiler* p = new Profiler();
    return p;
}

FLUPS_Profiler* flups_profiler_new_n(const char name[]){
    Profiler* p = new Profiler(name);
    return p;
}

void flups_profiler_free(FLUPS_Profiler* p) {
    delete p;
}

void flups_profiler_disp_root(FLUPS_Profiler* p) {
    p->disp();
}

void flups_profiler_disp(FLUPS_Profiler* p, const char* name) {
    const std::string myname(name);
    p->disp(myname);
}

//**********************************************************************
//  HDF5
//**********************************************************************

void flups_hdf5_dump(const FLUPS_Topology *topo, const char filename[], const double *data){
    const std::string fn(filename);
    hdf5_dump(topo,fn, data);
}


}
