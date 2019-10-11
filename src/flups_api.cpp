/**
 * @file flups_api.cpp
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

/** *********************************************************************
 * TOPOLOGIES
 ********************************************************************* */

FLUPS_Topology* flups_new_topo_(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3]){
    Topology* t = new Topology(axis, nglob, nproc, isComplex, axproc, 1);
    return t;
}

FLUPS_Topology* flups_new_topo_(const int axis, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment){
    Topology* t = new Topology(axis, nglob, nproc, isComplex, axproc, alignment);
    return t;
}

void        flups_free_topo(FLUPS_Topology* t){
    t->~Topology();
}


bool flups_topo_get_isComplex(FLUPS_Topology* t){
    return t->isComplex();
}

int flups_topo_get_axis(FLUPS_Topology* t){
    return t->axis();
}

int flups_topo_get_nglob(FLUPS_Topology* t, const int dim){
    return t->nglob(dim);
}

int flups_topo_get_nloc(FLUPS_Topology* t, const int dim){
    return t->nloc(dim);
}

int flups_topo_get_nmem(FLUPS_Topology* t, const int dim){
    return t->nmem(dim);
}

int flups_topo_get_nproc(FLUPS_Topology* t, const int dim){
    return t->nproc(dim);
}

void flups_topo_get_istartGlob(FLUPS_Topology* t, int istart[3]){
    t->get_istart_glob(istart);
}

unsigned long long flups_topo_get_locsize(FLUPS_Topology* t){
    return (unsigned long long) t->locsize();
}

unsigned long long flups_topo_get_memsize(FLUPS_Topology* t){
    return (unsigned long long) t->memsize();
}


/** *********************************************************************
 *  SOLVER
 ********************************************************************* */

// get a new solver
#ifndef PROF
FLUPS_Solver* flups_new_solver(FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3]){
    Solver* s = new Solver(t, bc, h, L, NULL);
    return s;
}
#else
FLUPS_Solver* flups_new_solver_timed(FLUPS_Topology* t, const FLUPS_BoundaryType bc[3][2], const double h[3], const double L[3],Profiler* prof){
    Solver* s = new Solver(t, bc, h, L, prof);
    return s;
}
#endif

// destroy the solver
void flups_free_solver(FLUPS_Solver* s){
    s->~Solver();
}

// setup the solver
void flups_set_greenType(FLUPS_Solver* s, const FLUPS_GreenType type){
    s->set_GreenType(type);
}

void flups_setup(FLUPS_Solver* s){
    s->setup();
}

// solve
void flups_solve(FLUPS_Solver* s, const FLUPS_Topology* t, double* field, double* rhs, const FLUPS_SolverType type){
    s->solve(t, field, rhs, type);
}


// -- ADVANCED FEATURES --

void flups_set_alpha(FLUPS_Solver* s, const double alpha){
    s->set_alpha(alpha);   
}

void flups_set_OrderDiff(FLUPS_Solver* s, const int order){
    s->set_OrderDiff(order);
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

// //**********************************************************************
// //  PROFILER - TIMERS
// //**********************************************************************

FLUPS_Profiler* profiler_new() {
    Profiler* p = new Profiler();
    return p;
}

FLUPS_Profiler* profiler_new_n(char name[]){
    Profiler* p = new Profiler(name);
    return p;
}

void profiler_free(FLUPS_Profiler* p) {
    p->~Profiler();
}

void profiler_disp(FLUPS_Profiler* p) {
    p->disp();
}
