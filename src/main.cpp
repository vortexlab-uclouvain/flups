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

#include "mpi.h"

typedef double pos[4];


int main(int argc, char* argv[])
{
    
    
    #if (DIM==2)
        // int nsample = 2; int size[2] = {64,128}; 
        int nsample = 3; int size[3] = {64,128,256};
       // int nsample = 4; int size[4] = {64,128,256,512};

        validation_2d_UU_UU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UU_UU(nsample,size,UP_SRHS,HEJ_4);

        validation_2d_UU_UE(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UU_UE(nsample,size,UP_SRHS,HEJ_4);
        validation_2d_UU_EU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UU_EU(nsample,size,UP_SRHS,HEJ_4);

        validation_2d_UE_UU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UE_UU(nsample,size,UP_SRHS,HEJ_4);
        validation_2d_EU_UU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_EU_UU(nsample,size,UP_SRHS,HEJ_4);

        validation_2d_UU_UO(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UU_UO(nsample,size,UP_SRHS,HEJ_4);
        validation_2d_UU_OU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UU_OU(nsample,size,UP_SRHS,HEJ_4);

        validation_2d_UO_UU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_UO_UU(nsample,size,UP_SRHS,HEJ_4);
        validation_2d_OU_UU(nsample,size,UP_SRHS,HEJ_2);
        validation_2d_OU_UU(nsample,size,UP_SRHS,HEJ_4);
    #elif (DIM == 3)
        MPI_Init(&argc, &argv);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        

        int nproc;
        MPI_Comm_size(MPI_COMM_WORLD,&nproc);

        int n = sqrt(nproc);


        int global_size[3] = {2,2,4};
        int local_size[3] = {global_size[0],global_size[1]/n,global_size[2]/n};
        pos* array =(pos*) fftw_malloc(local_size[0]*local_size[1]*local_size[2]*sizeof(pos));

        for(int i2=0; i2<local_size[2]; i2++){
            for(int i1=0; i1<local_size[1]; i1++){
                for(int i0=0; i0<local_size[0]; i0++){
                    array[i0 + local_size[0]*i1][0] = i0;
                    array[i0 + local_size[0]*i1][1] = i1;
                    array[i0 + local_size[0]*i1][2] = i2;
                    array[i0 + local_size[0]*i1][3] = rank;
                }
            }
        }

        ptrdiff_t local_n[3];
        ptrdiff_t start_n[3];
        fftw_mpi_local_size_3d(global_size[1],global_size[1],global_size[0],MPI_COMM_WORLD,local_n,start_n);
        printf("rank = %d - advised size = %d\n",rank,local_n[0]);

        // fftw_mpi_local_size_3d()

        // printf("Hello from rank %d/%d, create a grid of size %d %d %d\n",rank,nproc,local_size[0],local_size[1],local_size[2]);

        // printf("creating the transpose of a %d by %d matrix\n",global_size[0],global_size[1]*global_size[2]);
        // printf("before the size of %d is split in blocs of %d\n",global_size[1]*global_size[2],local_size[1]*local_size[2]);

        // fftw_plan transplan = fftw_mpi_plan_many_transpose(nproc,nproc,4,1,1,(double*) array,(double*) array,MPI_COMM_WORLD,FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_OUT);

        // for(int i2=0; i2<local_size[2]; i2++){
        //     for(int i1=0; i1<local_size[1]; i1++){
        //         for(int i0=0; i0<local_size[0]; i0++){
        //             printf("rank = %d, value = (%f, %f, %f)\n",rank,array[i0 + local_size[0]*i1][0],array[i0 + local_size[0]*i1][1],array[i0 + local_size[0]*i1][2]);
        //         }
        //     }
        // }

        // // fftw_execute(transplan);

        // MPI_Barrier(MPI_COMM_WORLD);

        

        // printf("\n\nTRANSPOSE\n\n");


        // for(int i1=0; i1<local_size[0]; i1++){
        //     for(int i0=0; i0<local_size[1]; i0++){
        //         for(int i0=0; i0<local_size[1]; i0++){
        //             printf("rank = %d, value = (%f, %f, %f)\n",rank,array[i0 + local_size[1]*i1][0],array[i0 + local_size[1]*i1][1],array[i0 + local_size[1]*i1][2]);
        //         }
        //     }
        // }

        MPI_Finalize();
    #endif
}