/**
 * @file flups.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-10-09
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
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

#include "defines.hpp"
#include "Topology.hpp"
#include "Solver.hpp"
#include "Profiler.hpp"


extern "C" {

void * flups_malloc(FLUPS_SIZE size){
    return flups_mem_malloc(size);
}

void flups_free(void* data){
    flups_mem_free(data);
}


//***********************************************************************
// * TOPOLOGIES
// **********************************************************************/
const FLUPS_Topology* flups_topo_new(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment){
    Topology* t = new Topology(axis, nglob, nproc, isComplex, axproc, alignment);
    return t;
}

void flups_topo_free(const FLUPS_Topology* t) {
    t->~Topology();
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

FLUPS_SIZE flups_topo_get_locsize(const FLUPS_Topology* t) {
    return (FLUPS_SIZE)t->locsize();
}

FLUPS_SIZE flups_topo_get_memsize(const FLUPS_Topology* t) {
    return (FLUPS_SIZE)t->memsize();
}

//***********************************************************************
//*  SOLVER
//********************************************************************* */

// get a new solver

FLUPS_Solver* flups_init(const FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]){
    Solver* s = new Solver(t, bc, h, L, NULL);
    return s;
}
FLUPS_Solver* flups_init_timed(const FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3],Profiler* prof){
    Solver* s = new Solver(t, bc, h, L, prof);
    return s;
}

// destroy the solver
void flups_cleanup(FLUPS_Solver* s){
    s->~Solver();
}

// setup the solver
void flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type){
    s->set_GreenType(type);
}

double* flups_setup(FLUPS_Solver* s){
    return s->setup();
}

// solve
void flups_solve(FLUPS_Solver* s, double* field, double* rhs, const FLUPS_SolverType type){
    s->solve(field, rhs, type);
}


// -- ADVANCED FEATURES --

FLUPS_SIZE flups_get_allocSize(FLUPS_Solver* s){
    return(FLUPS_SIZE) s->get_allocSize();
}

void flups_set_alpha(FLUPS_Solver* s, const double alpha){
    s->set_alpha(alpha);   
}

void flups_set_OrderDiff(FLUPS_Solver* s, const int order){
    s->set_OrderDiff(order);
}

const FLUPS_Topology* flups_get_topo_physical(FLUPS_Solver* s){
    return s->get_innerTopo_physical();
}

const FLUPS_Topology* flups_get_topo_spectral(FLUPS_Solver* s){
    return s->get_innerTopo_spectral();
}

void flups_do_copy(FLUPS_Solver* s, const FLUPS_Topology* topo, double* data, const int sign){
    s->do_copy(topo,data,sign);
}

void flups_do_FFT(FLUPS_Solver* s, double* data, const int sign){
    s->do_FFT(data,sign);
}

void flups_do_mult(FLUPS_Solver* s, double* data, const FLUPS_SolverType type){
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
    p->~Profiler();
}

void flups_profiler_disp(FLUPS_Profiler* p) {
    p->disp();
}

void flups_profiler_disp_root(FLUPS_Profiler* p, const char* name) {
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
