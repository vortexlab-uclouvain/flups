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

class Reorder_MPI
{

protected:
    int _nf    = 1;
    int _axis0 = 1;
    int _axis1 = 1;

    int _inloc[3];
    int _onloc[3];
    int _inproc[3];
    int _onproc[3];
    int _irankd[3];
    int _orankd[3];

    int *_nsend = NULL; // number of unknowns send to each proc
    int *_nrecv = NULL; // number of unknowns received from each proc
    int *_ssend = NULL; // start index in my memory to send to each proc
    int *_srecv = NULL; // start index in my memory to receive from each proc
    int *_count = NULL;

    double *_bufsend = NULL;
    double *_bufrecv = NULL;

public:
    Reorder_MPI(const int nglob[3], const int nf, const int inproc[3], int axis0, int onproc[3], int axis1);
    ~Reorder_MPI();

    void execute(double *v);
};

#endif