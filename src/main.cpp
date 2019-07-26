/**
 * @file main.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-16
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include <cmath>
#include <iostream>

#include "expint.hpp"
#include "tools.hpp"

#include "FFTW_Solver.hpp"
#include "validation_2d.hpp"
#include "validation_3d.hpp"
#include "Reorder_MPI.hpp"

#include "mpi.h"

typedef double pos[4];

int main(int argc, char *argv[])
{

#if (DIM == 2)
    // int nsample = 2; int size[2] = {64,128};
    int nsample = 3;
    int size[3] = {64, 128, 256};
    // int nsample = 4; int size[4] = {64,128,256,512};

    validation_2d_UU_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_UU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UU_UE(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_UE(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_UU_EU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_EU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UE_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UE_UU(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_EU_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_EU_UU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UU_UO(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_UO(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_UU_OU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UU_OU(nsample, size, UP_SRHS, HEJ_4);

    validation_2d_UO_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_UO_UU(nsample, size, UP_SRHS, HEJ_4);
    validation_2d_OU_UU(nsample, size, UP_SRHS, HEJ_2);
    validation_2d_OU_UU(nsample, size, UP_SRHS, HEJ_4);
#elif (DIM == 3)
    MPI_Init(&argc, &argv);


    int n[3] = {16,16,16};
    
    int nproc[3] = {2, 2, 2};
    size_t nloc[3];
    for (int i = 0; i < 3; ++i)
    {
        nloc[i] = n[i] / nproc[i];
    }
    int rank;

    int proc_px[3] = {1,2,4};
    int proc_py[3] = {2,1,4};
    int proc_pz[3] = {4,2,1};

    Reorder_MPI* reorder_g2x = new Reorder_MPI(n, 1, nproc  ,0,proc_px,0);
    Reorder_MPI* reorder_x2y = new Reorder_MPI(n, 1, proc_px,0,proc_py,1);
    Reorder_MPI* reorder_y2z = new Reorder_MPI(n, 1, proc_py,1,proc_pz,2);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rankd[3] = {
        rank / (nproc[1] * nproc[2]),
        (rank % (nproc[1] * nproc[2])) / nproc[2],
        rank % nproc[2]};

    double *v =(double*) fftw_malloc(sizeof(double) * nloc[0] * nloc[1] * nloc[2]);
    for (size_t i = 0; i < nloc[0]; i++)
    {
        for (size_t j = 0; j < nloc[1]; j++)
        {
            for (size_t k = 0; k < nloc[2]; k++)
            {
                double x = 1. / n[0] * (i + rankd[0] * nloc[0]);
                double y = 1. / n[1] * (j + rankd[1] * nloc[1]);
                double z = 1. / n[2] * (k + rankd[2] * nloc[2]);
                v[(i * nloc[1] + j) * nloc[2] + k] = z;
            }
        }
    }
    
    for (size_t i = 0; i < nloc[0]; i++)
    {
        for (size_t j = 0; j < nloc[1]; j++)
        {
            for (size_t k = 0; k < nloc[2]; k++)
            {
                printf("%lu %lu %lu : %g\n", i, j, k, v[(i * nloc[1] + j) * nloc[2] + k]);
            }
        }
    }
    reorder_g2x->execute(v);
    // printf("x to y\n");
    reorder_x2y->execute( v);
    
    // for (size_t i = 0; i < nloc[0]; i++)
    // {
    //     for (size_t j = 0; j < nloc[1]; j++)
    //     {
    //         for (size_t k = 0; k < nloc[2]; k++)
    //         {
    //             printf("%lu %lu %lu : %g\n", i, j, k, v[(i * nloc[1] + j) * nloc[2] + k]);
    //         }
    //     }
    // }
    // printf("y to z\n");
    reorder_y2z->execute( v);

    for (size_t i = 0; i < nloc[0]; i++)
    {
        for (size_t j = 0; j < nloc[1]; j++)
        {
            for (size_t k = 0; k < nloc[2]; k++)
            {
                printf("%lu %lu %lu : %g\n", i, j, k, v[(i * nloc[1] + j) * nloc[2] + k]);
            }
        }
    }
    fftw_free(v);

    delete reorder_g2x;
    delete reorder_x2y;
    delete reorder_y2z;
    
    MPI_Finalize();
#endif
}