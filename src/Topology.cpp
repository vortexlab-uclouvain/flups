/**
 * @file Topology.cpp
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

#include "Topology.hpp"

/**
 * @brief Construct a new Topology
 * 
 * @param nf the number of field (eg: real = 1, complex = 2)
 * @param lda leading dimension of array=the number of components (eg scalar=1, vector=3)
 * @param axis the dimension of the fastest rotating index (eg: 0, 1 or 2)
 * @param nglob the global size per dim
 * @param nproc the number of proc per dim
 * @param isComplex indicate if the Topology uses complex indexing or not
 * @param axproc gives the order of the rank decomposition (eg. (0,2,1) to start decomposing in X then Z then Y). If NULL is passed, use by default (0,1,2).
 * @param alignment the number of bytes on which we want the topology to be aligned along the #axis only
 * @param comm the communicator associated to the topology.
 * 
 * If the MPI comm is associated with a MPI_CART topology, axproc is ignored and we use the MPI routines to determine the 3D rank from the global rank (and vice versa).
 * 
 */
Topology::Topology(const int axis, const int lda, const int nglob[3], const int nproc[3], const bool isComplex, const int axproc[3], const int alignment, MPI_Comm comm):_alignment(alignment) {
    BEGIN_FUNC;

    _comm = comm;
    _lda  = lda;

    int comm_size, rank;
    MPI_Comm_size(_comm,&comm_size);
    MPI_Comm_rank(_comm,&rank);

    FLUPS_CHECK(nproc[0]*nproc[1]*nproc[2] == comm_size,"the total number of procs (=%d) have to be = to the comm size (=%d)",nproc[0]*nproc[1]*nproc[2], comm_size, LOCATION);

    //-------------------------------------------------------------------------
    /** - get memory axis and complex information  */
    //-------------------------------------------------------------------------
    _axis = axis;
    if (!isComplex) {
        _nf = 1;
    } else {
        _nf = 2;
    }

    //-------------------------------------------------------------------------
    /** - get proc information  */
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        // total dimension
        _nglob[id] = nglob[id];
        // number of proc in each direction
        _nproc[id] = nproc[id];
        // store the proc axis, used to split the rank
        _axproc[id] = (axproc == NULL) ? id : axproc[id];
    }

    //-------------------------------------------------------------------------
    /** - split the rank and get rankd  */
    //-------------------------------------------------------------------------
    ranksplit(rank, _axproc, _nproc, _comm, _rankd);

    //-------------------------------------------------------------------------
    /** - Get the new sizes  */
    //-------------------------------------------------------------------------
    cmpt_sizes();

    FLUPS_INFO("New topo created: nf = %d, axis = %d, local sizes = %d %d %d vs mem size = %d %d %d",_nf,_axis,_nloc[0],_nloc[1],_nloc[2],_nmem[0],_nmem[1],_nmem[2]);
    END_FUNC;
}

/**
 * @brief compute the nloc and nmem sizes using _rankd, _nglob, _nproc, _nloc
 * 
 * This function padds the size of the domain if needed
 * 
 */
void Topology::cmpt_sizes() {
    BEGIN_FUNC;
    for (int id = 0; id < 3; id++) {
        // we get the max between the nglob and
        _nloc[id] = cmpt_nbyproc(id);
        _nmem[id] = _nloc[id];
        // if we are in the axis and the last proc, we pad to ensure that every pencil is ok with alignment
        // if (id == _axis && _rankd[id] == (_nproc[id] - 1)) {
        if (id == _axis) {
            // compute by how many we are not aligned: the global size in double = nglob * nf
            const int modulo = (_nloc[id] * _nf * sizeof(double)) % _alignment;
            // compute the number of points to add (in double indexing)
            const int delta = (_alignment - modulo) / sizeof(double);
            _nmem[id] += (modulo == 0) ? 0 : delta / _nf;
        }
    }
    END_FUNC;
}

/**
 * @brief Set a new communicator for the topology
 * 
 * @param comm 
 */
void Topology::change_comm(MPI_Comm comm) {
    BEGIN_FUNC;

    //-------------------------------------------------------------------------
    /** - Store the new communicator  */
    //-------------------------------------------------------------------------
    _comm = comm;
    //-------------------------------------------------------------------------
    /** - Recompute the new rank and rankd info  */
    //-------------------------------------------------------------------------
    int newrank;
    MPI_Comm_rank(comm, &newrank);
    // split the rank to get the new rankd
    ranksplit(newrank, _axproc, _nproc, _comm, _rankd);
    //-------------------------------------------------------------------------
    /** - Get the sizes  */
    //-------------------------------------------------------------------------
    cmpt_sizes();

    END_FUNC;
}

/**
 * @brief Destroy the Topology
 * 
 */
Topology::~Topology() {}

/**
 * @brief compute the sarting and ending ids for the current topo in order to be inside the other topology's bounds
 * 
 * @param shift the shift between the 2 topos: current topo in (0,0,0) = other topo in (shift)
 * @param other the other topology
 * @param start the start index to use in the current topo
 * @param end the end index to use in the current topo
 */
void Topology::cmpt_intersect_id(const int shift[3], const Topology* other, int start[3], int end[3]) const {
    BEGIN_FUNC;
    FLUPS_CHECK(this->isComplex() == other->isComplex(), "The two topo have to be both complex or real", LOCATION);

    for (int id = 0; id < 3; id++) {
        const int onglob = other->nglob(id);

        start[id] = 0;
        end[id]   = 0;
        // for the input configuration
        for (int i = 0; i < _nloc[id]; ++i) {
            // get the global id in the other topology
            int oid_global = cmpt_start_id(id) + i + shift[id];
            if (oid_global <= 0) start[id] = i;
            if (oid_global < onglob) end[id] = i + 1;
        }
    }
    END_FUNC;
}

/**
 * @brief display topology most important info
 * 
 */
void Topology::disp() const {
    BEGIN_FUNC;

    int comm_size, rank;
    MPI_Comm_size(_comm,&comm_size);
    MPI_Comm_rank(_comm,&rank);

    FLUPS_INFO("------------------------------------------");
    FLUPS_INFO("## Topology created on proc %d/%d", rank, comm_size);
    FLUPS_INFO(" - axis = %d", _axis);
    FLUPS_INFO(" - nglob = %d %d %d", _nglob[0], _nglob[1], _nglob[2]);
    FLUPS_INFO(" - nloc = %d %d %d", _nloc[0], _nloc[1], _nloc[2]);
    FLUPS_INFO(" - nmem = %d %d %d", _nmem[0], _nmem[1], _nmem[2]);
    FLUPS_INFO(" - nproc = %d %d %d", _nproc[0], _nproc[1], _nproc[2]);
    FLUPS_INFO(" - rankd = %d %d %d", _rankd[0], _rankd[1], _rankd[2]);
    // FLUPS_INFO(" - nbyproc = %d %d %d", _nbyproc[0], _nbyproc[1], _nbyproc[2]);
    FLUPS_INFO(" - axproc = %d %d %d", _axproc[0], _axproc[1], _axproc[2]);
    FLUPS_INFO(" - isComplex = %d", _nf == 2);
    // FLUPS_INFO(" - h = %f %f %f",_h[0],_h[1],_h[2]);
    // FLUPS_INFO(" - L = %f %f %f",_L[0],_L[1],_L[2]);
    FLUPS_INFO("------------------------------------------");
}

void Topology::disp_rank() {
    BEGIN_FUNC;
#ifdef DUMP_DBG
    // we only focus on the real size = local size
    double* rankdata = (double*)flups_malloc(sizeof(double) * this->locsize() * 2);
    int     rank, rank_new;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(_comm, &rank_new);

    for (int i = 0; i < this->locsize(); i++) {
        rankdata[2 * i]     = rank + rank_new / 100.;
        rankdata[2 * i + 1] = _rankd[0] + _rankd[1] / 10. + _rankd[2] / 100.;
    }

    int  rlen;
    char commname[MPI_MAX_OBJECT_NAME];
    MPI_Comm_get_name(_comm, commname, &rlen);
    std::string cn(commname, rlen);
    std::string name = "rank_topo_axis" + std::to_string(this->axis()) + "_procs" + std::to_string(this->nproc(0)) + std::to_string(this->nproc(1)) + std::to_string(this->nproc(2)) + "_" + cn;
    if (this->isComplex()) {
        hdf5_dump(this, name, rankdata);
    } else {
        this->switch2complex();
        hdf5_dump(this, name, rankdata);
        this->switch2real();
    }
    flups_free(rankdata);
#endif
    END_FUNC;
}

/**
 * @brief this function shift the current axis by 1 element in the specified component
 * 
 * FLUPS_FORWARD = (-1) shift:
 *      +---+---+---+---+---+---+
 *      | 0 | 1 | 2 | 3 | 4 | 5 |
 *      +---+---+---+---+---+---+
 *  +---+---+---+---+---+---+
 *  | 0 | 1 | 2 | 3 | 4 | 5 |
 *  +---+---+---+---+---+---+
 * 
 * FLUPS_BACKWARD = (+1) shift:
 * +---+---+---+---+---+---+
 * | 0 | 1 | 2 | 3 | 4 | 5 |
 * +---+---+---+---+---+---+
 *     +---+---+---+---+---+---+
 *     | 0 | 1 | 2 | 3 | 4 | 5 |
 *     +---+---+---+---+---+---+
 * 
 * @param sign 
 * @param lia 
 * @param data 
 */
void Topology::memshift(const int sign,const int lia, double* data){
    BEGIN_FUNC;
    FLUPS_CHECK(_nproc[_axis] == 1, "You shouldn't shift the memory in a non continuous direction", LOCATION);

    // get the current axis
    const int nf  = _nf;
    const int lda = _lda;
    const int ax0 = _axis;
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    // compute some of the indexes
    const int    nmem[3] = {_nmem[0], _nmem[1], _nmem[2]};
    const size_t sonmax   = _nloc[ax1] * _nloc[ax2] * lia;
    const size_t eonmax   = _nloc[ax1] * _nloc[ax2] * (lia+1);
    const size_t inmax   = _nloc[ax0] * nf;

    // put the memory in the correct place
    if (sign == FLUPS_FORWARD) {
        FLUPS_INFO("runing memshift with FLUPS_FORWARD and nf = %d",nf);

        if (nf == 1) {
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(sonmax,eonmax, inmax, data, nmem, ax0)
            for (int io = sonmax; io < eonmax; io++) {
                opt_double_ptr dataloc = data + collapsedIndex(ax0, 0, io, nmem, 1);
                // set the alignment
                FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                // if forward we go from nloc-2 to 0 included (sign = -1)
                for (int ii = (inmax - 2); ii >= 0; ii--) {
                    // the source id is always the current one
                    dataloc[ii + 1] = dataloc[ii];
                }
                dataloc[0] = 0.0;
            }
        } else if (nf == 2) {
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, inmax, data, nmem, ax0)
            for (int io = sonmax; io < eonmax; io++) {
                opt_double_ptr dataloc = data + collapsedIndex(ax0, 0, io, nmem, 2);
                FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                // if forward we go from nloc-2 to 0 included (sign = -1)
                for (int ii = (inmax - 3); ii >= 0; ii--) {
                    // the source id is always the current one
                    dataloc[ii + 2] = dataloc[ii];
                }
                dataloc[0] = 0.0;
                dataloc[1] = 0.0;
            }
        }
    } else if (sign == FLUPS_BACKWARD) {
        FLUPS_INFO("runing memshift with FLUPS_FORWARD and nf = %d",nf);
        FLUPS_INFO("nloc = %d %d %d, onmax = %d, inmax = %d",_nloc[ax0],_nloc[ax1],_nloc[ax2],onmax,inmax);
        if (nf == 1) {
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, inmax, data, nmem, ax0)
            for (int io = sonmax; io < eonmax; io++) {
                opt_double_ptr dataloc = data + collapsedIndex(ax0, 0, io, nmem, 1);
                FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                // if forward we go from nloc-2 to 0 included (sign = -1)
                for (int ii = 1; ii < inmax; ii++) {
                    // the source id is always the current one
                    dataloc[ii - 1] = dataloc[ii];
                }
                dataloc[inmax-1] = 0.0;
            }
        } else if (nf == 2) {
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, inmax, data, nmem, ax0)
            for (int io = sonmax; io < eonmax; io++) {
                opt_double_ptr dataloc = data + collapsedIndex(ax0, 0, io, nmem, 2);
                FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
                // if forward we go from nloc-2 to 0 included (sign = -1)
                for (int ii = 2; ii < inmax; ii++) {
                    // the source id is always the current one
                    dataloc[ii - 2] = dataloc[ii];
                }
                dataloc[inmax-1] = 0.0;
                dataloc[inmax-2] = 0.0;
            }
        }
    }
    END_FUNC;
}