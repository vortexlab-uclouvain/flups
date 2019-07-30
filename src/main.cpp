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

    MPI_Init(&argc,&argv);
    int nsample = 3;
    int size[3] = {16, 32, 64};
    validation_3d_UU_UU(nsample,size,UP_SRHS,HEJ_2);
    MPI_Finalize();


    // MPI_Init(&argc, &argv);
    // int n[3] = {4,3,5};
    // int nproc[3] = {2, 1, 2};
    // int rank;
    // int proc_px[3] = {1,2,2};
    // int proc_py[3] = {2,1,2};
    // int proc_pz[3] = {2,2,1};
    // int ishift[3] = {0,0,0};

    // Topology* topo_glob = new Topology(0,n,nproc);
    // Topology* topo_penx = new Topology(0,n,proc_px);
    // Topology* topo_peny = new Topology(1,n,proc_py);
    // Topology* topo_penz = new Topology(2,n,proc_pz);

    // Reorder_MPI* reorder_g2x = new Reorder_MPI(1,topo_glob,topo_penx,ishift);
    // Reorder_MPI* reorder_x2y = new Reorder_MPI(1,topo_penx,topo_peny,ishift);
    // Reorder_MPI* reorder_y2z = new Reorder_MPI(1,topo_peny,topo_penz,ishift);

    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // int alloc_size = std::max(std::max(topo_glob->locsize(),topo_penx->locsize()),std::max(topo_peny->locsize(),topo_penz->locsize()));
    // double *v =(double*) fftw_malloc(sizeof(double) * alloc_size);

    // int istart[3];
    // topo_glob->get_idstart_glob(istart);

    // for (size_t i = 0; i < topo_glob->nloc(2); i++)
    // {
    //     for (size_t j = 0; j < topo_glob->nloc(1); j++)
    //     {
    //         for (size_t k = 0; k < topo_glob->nloc(0); k++)
    //         {
    //             double x = 1. / n[0] * (i + istart[0]);
    //             double y = 1. / n[1] * (j + istart[1]);
    //             double z = 1. / n[2] * (k + istart[2]);
    //             v[(i * topo_glob->nloc(1) + j) * topo_glob->nloc(0)+ k] = x;
    //         }
    //     }
    // }
    // for (size_t i = 0; i < topo_glob->nloc(2); i++)
    // {
    //     for (size_t j = 0; j < topo_glob->nloc(1); j++)
    //     {
    //         for (size_t k = 0; k < topo_glob->nloc(0); k++)
    //         {
    //             double x = 1. / n[0] * (i + istart[0]);
    //             double y = 1. / n[1] * (j + istart[1]);
    //             double z = 1. / n[2] * (k + istart[2]);
    //             printf("%lu %lu %lu : %g\n", i, j, k, v[(i * topo_glob->nloc(1) + j) * topo_glob->nloc(0)+ k]);
    //         }
    //     }
    // }
    
    // reorder_g2x->execute(v,FFTW_FORWARD);
    // // printf("x to y\n");
    // reorder_x2y->execute(v,FFTW_FORWARD);
    
    // // for (size_t i = 0; i < nloc[0]; i++)
    // // {
    // //     for (size_t j = 0; j < nloc[1]; j++)
    // //     {
    // //         for (size_t k = 0; k < nloc[2]; k++)
    // //         {
    // //             printf("%lu %lu %lu : %g\n", i, j, k, v[(i * nloc[1] + j) * nloc[2] + k]);
    // //         }
    // //     }
    // // }
    // // printf("y to z\n");
    // reorder_y2z->execute(v,FFTW_FORWARD);

    // // for (size_t i2 = 0; i2 < topo_penz->nloc(2); i2++)
    // // {
    // //     for (size_t i1 = 0; i1 < topo_penz->nloc(1); i1++)
    // //     {
    // //         for (size_t i0 = 0; i0 < topo_penz->nloc(0); i0++)
    // //         {
    // //             double x = 1. / n[0] * (i0 + istart[0]);
    // //             double y = 1. / n[1] * (i1 + istart[1]);
    // //             double z = 1. / n[2] * (i2 + istart[2]);
    // //             const int id[3] = {i0,i1,i2};
    // //             printf("%lu %lu %lu : %g\n", i0, i1, i2, v[localindex(id,topo_penz)]);
    // //         }
    // //     }
    // // }

    // reorder_y2z->execute(v,FFTW_BACKWARD);
    // reorder_x2y->execute(v,FFTW_BACKWARD);
    // reorder_g2x->execute(v,FFTW_BACKWARD);


    // for (size_t i2 = 0; i2 < topo_glob->nloc(2); i2++)
    // {
    //     for (size_t i1 = 0; i1 < topo_glob->nloc(1); i1++)
    //     {
    //         for (size_t i0 = 0; i0 < topo_glob->nloc(0); i0++)
    //         {
    //             double x = 1. / n[0] * (i0 + istart[0]);
    //             double y = 1. / n[1] * (i1 + istart[1]);
    //             double z = 1. / n[2] * (i2 + istart[2]);
    //             const int id[3] = {i0,i1,i2};
    //             printf("%lu %lu %lu : %g\n", i0, i1, i2, v[localindex(id,topo_glob)]);
    //         }
    //     }
    // }
    // fftw_free(v);

    // delete reorder_g2x;
    // delete reorder_x2y;
    // delete reorder_y2z;
    // delete topo_glob;
    // delete topo_penx;
    // delete topo_peny;
    // delete topo_penz;
    
    // MPI_Finalize();
#endif
}