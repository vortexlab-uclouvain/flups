/**
 * @file Reorder_MPI.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-25
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#ifndef REODER_MPI_HPP
#define REODER_MPI_HPP

#include "defines.hpp"
#include "mpi.h"
#include <cstring>
#include "topology.hpp"


class Reorder_MPI
{

protected:
    // int _nf   ;

    int _ishift[3];
    int _oshift[3];

    int _istart[3];
    int _iend[3];
    int _ostart[3];
    int _oend[3];

    int *_nsend = NULL; // number of unknowns send to each proc
    int *_nrecv = NULL; // number of unknowns received from each proc
    int *_ssend = NULL; // start index in my memory to send to each proc
    int *_srecv = NULL; // start index in my memory to receive from each proc
    int *_count = NULL;

    double *_bufsend = NULL;
    double *_bufrecv = NULL;

    const Topology* _topo_in;
    const Topology* _topo_out;

public:
    // Reorder_MPI(const int mynf, const int nglob[3], const int inproc[3], int myaxis0, int onproc[3], int myaxis1);
    Reorder_MPI(const Topology* topo_input,const Topology* topo_output, const int shift[3]);
    ~Reorder_MPI();

    void execute(opt_double_ptr v,const int sign);
};

#endif