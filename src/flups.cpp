/**
 * @file flups.cpp
 * @copyright Copyright © UCLouvain 2020
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "flups.h"

#include <algorithm>  // std::remove

#include "Solver.hpp"
#include "Topology.hpp"
#include "FFTW_plan_dim.hpp"
#include "defines.hpp"
#include "h3lpr/profiler.hpp"

extern "C" {

void* flups_malloc(size_t size) {
    return m_calloc(size); 
}

void flups_free(void* data) {
    m_free(data);
}

//***********************************************************************
// * TOPOLOGIES
// **********************************************************************/
Topology* flups_topo_new(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm) {
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

int flups_topo_cmpt_rank_fromid(const Topology* t, const int global_id, const int id) {
    return t->cmpt_rank_fromid(global_id, id);
}

int flups_topo_cmpt_start_id_from_rank(const Topology* t, const int rank_id, const int id) {
    return t->cmpt_start_id_from_rank(rank_id, id);
}

MPI_Comm flups_topo_get_comm(Topology* t) {
    return t->get_comm();
}

void flups_topo_ranksplit(const Topology* t, const int rank, int rankd[3]) {
    int axproc[3] = {t->axproc(0), t->axproc(1), t->axproc(2)};
    int nproc[3]  = {t->nproc(0), t->nproc(1), t->nproc(2)};
    ranksplit(rank, axproc, nproc, t->get_comm(), rankd);
}

int flups_topo_rankindex(const Topology* topo, const int rankd[3]) {
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
Solver* flups_init_timed(Topology* t, BoundaryType* bc[3][2], const double h[3], const double L[3], const DiffType orderDiff, const CenterType center_type[3], FLUPS_Profiler* prof) {
    auto    h_prof = reinterpret_cast<H3LPR::Profiler*>(prof);
    Solver* s      = new Solver(t, bc, h, L, orderDiff, center_type, h_prof);
    return s;
}

// destroy the solver
void flups_cleanup(Solver* s) {
    delete s;
}

// setup the solver
void flups_set_greenType(Solver* s, const GreenType type) {
    s->set_GreenType(type);
}

void flups_setup(Solver* s, const bool changeComm) {
    s->setup(changeComm);
}

// solve
void flups_solve(Solver* s, double* field, double* rhs, const SolverType type) {
    s->solve(field, rhs, type);
}

// -- ADVANCED FEATURES --

size_t flups_get_allocSize(Solver* s) {
    return (size_t)s->get_allocSize();
}

void flups_get_spectralInfo(Solver* s, double kfact[3], double koffset[3], double symstart[3]) {
    s->get_spectralInfo(kfact, koffset, symstart);
}

void flups_set_alpha(Solver* s, const double alpha) {
    s->set_alpha(alpha);
}

double* flups_get_innerBuffer(FLUPS_Solver* s){
    return s->get_innerBuffer();
}

Topology* flups_get_innerTopo_physical(Solver* s) {
    return s->get_innerTopo_physical();
}

Topology* flups_get_innerTopo_spectral(Solver* s) {
    return s->get_innerTopo_spectral();
}

void flups_skip_firstSwitchtopo(Solver* s){
    s->skip_firstSwitchtopo();
}

void flups_do_copy(Solver* s, const Topology* topo, double* data, const int sign) {
    s->do_copy(topo, data, sign);
}

void flups_do_FFT(Solver* s, double* data, const int sign) {
    s->do_FFT(data, sign);
}

void flups_do_mult(Solver* s, double* data, const SolverType type) {
    s->do_mult(data, type);
}

void flups_pencilDirs(const FLUPS_BoundaryType* bc[3][2], int dirs[3]) {
    // get the priorities from the bcs
    std::array<std::tuple<int, int>, 3> priority = {std::make_tuple(bc_to_types(bc[0]), 0),
                                                    std::make_tuple(bc_to_types(bc[1]), 1),
                                                    std::make_tuple(bc_to_types(bc[2]), 2)};
    // sort the priority and get the direction of the pencils
    sort_priority(&priority);
    for (int id = 0; id < 3; ++id) {
        dirs[id] = std::get<1>(priority[id]);
    }
}

void flups_switchtopo_info(Solver* s) {
    s->get_switchtopo_info();
}

//**********************************************************************
//  PROFILER - TIMERS
//**********************************************************************
FLUPS_Profiler* flups_profiler_new() {
    H3LPR::Profiler* p = new H3LPR::Profiler();
    return reinterpret_cast<FLUPS_Profiler*>(p);
}

FLUPS_Profiler* flups_profiler_new_n(const char name[]) {
    H3LPR::Profiler* p = new H3LPR::Profiler(name);
    return reinterpret_cast<FLUPS_Profiler*>(p);
}

void flups_profiler_free(FLUPS_Profiler* p) {
    H3LPR::Profiler* h3lpr_p = reinterpret_cast<H3LPR::Profiler*>(p);
    delete h3lpr_p;
}

void flups_profiler_disp(FLUPS_Profiler* p) {
    H3LPR::Profiler* h3lpr_p = reinterpret_cast<H3LPR::Profiler*>(p);
    m_profDisp(h3lpr_p);
}

//**********************************************************************
//  HDF5
//**********************************************************************
void flups_hdf5_dump(const Topology* topo, const char filename[], const double* data) {
    const std::string fn(filename);
    hdf5_dump(topo, fn, data);
}

void flups_print_data(const Topology* topo, double* data) {
    FLUPS_print_data(topo, data);
}

//**********************************************************************
//  MISC
//**********************************************************************
/**
 * @brief writes the file murphy.info used for tracking of the results, bookkeeping etc
 */
void flups_info(int argc, char** argv) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    if (rank == 0) {
        std::string filename(argv[0]);
        // remove the "./" that might be in it
        filename.erase(std::remove(filename.begin(), filename.end(), '.'), filename.end());
        filename.erase(std::remove(filename.begin(), filename.end(), '/'), filename.end());
        filename += ".flups.info";
        FILE* file = fopen(filename.c_str(), "w+");
        fprintf(file, "FLUPS \n");
        fprintf(file, "- commit: %s\n", FLUPS_GIT_COMMIT);
        fprintf(file, "- comm_world size: %d\n", comm_size);
        fprintf(file, "- defines:\n");
        fprintf(file, "\tFLUPS_FFTW_FLAG = %d (estimate: %d, measure: %d, patient: %d)\n", FLUPS_FFTW_FLAG, FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT);
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
        fprintf(file, "\tPersistent, non blocking implementation \n");
        fprintf(file, "\tFLUPS_MPI_BATCH_SEND = %d\n", FLUPS_MPI_BATCH_SEND);
        fprintf(file, "\tFLUPS_MPI_MAX_NBSEND = %d\n", FLUPS_MPI_MAX_NBSEND);
#endif
#ifdef COMM_ISR
        fprintf(file, "\tNon blocking implementation -- MPI data type \n");
        fprintf(file, "\tFLUPS_MPI_BATCH_SEND = %d\n", FLUPS_MPI_BATCH_SEND);
        fprintf(file, "\tFLUPS_MPI_MAX_NBSEND = %d\n", FLUPS_MPI_MAX_NBSEND);
#endif
#if (FLUPS_HDF5)
        fprintf(file, "\tHDF5 ? yes\n");
#else
        fprintf(file, "\tHDF5 ? no\n");
#endif
#if (FLUPS_PRIORITYLIST)
        fprintf(file, "\tPriority list order ? yes\n");
#else
        fprintf(file, "\tPriority list order ? no\n");
#endif

#if (FLUPS_ROLLING_RANK)
        fprintf(file, "\tRolling rank ? yes\n");
#else
        fprintf(file, "\tRolling rank  ? no\n");
#endif

#if (FLUPS_MPI_ALLOC)
        fprintf(file, "\tMPI ALLOC ? yes\n");
#else
        fprintf(file, "\tMPI ALLOC ? no\n");
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

}  // end extern
