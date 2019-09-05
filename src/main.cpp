/**
 * @file main.cpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
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
#include "SwitchTopo.hpp"
#include "Validation_2d.hpp"
#include "Validation_3d.hpp"

#include "mpi.h"


int main(int argc, char *argv[]) {

    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // int nsample = 3; int size[3] = {16, 32, 64};
    int nsample = 4; int size[4] = {16, 32, 64,128};
    // int nsample = 3; int size[3] = {64,128,256};
    // int nsample = 1; int size[1] = {128};
    // int nsample = 4; int size[4] = {128,128,128,128};

    // loop over the resolution for convergence study
    for (int is = 0; is < nsample; is++) {

        struct DomainDescr valCase ;
        valCase.nglob[0] = size[is] * valCase.L[0];
        valCase.nglob[1] = size[is] * valCase.L[1];
        valCase.nglob[2] = size[is] * valCase.L[2];

        if(rank==0) printf("\n==============================  FULLY UNBOUNDED TEST ============================================\n");
        valCase.mybc[0][0] = PER; valCase.mybc[0][1] = PER;
        valCase.mybc[1][0] = UNB; valCase.mybc[1][1] = UNB;
        valCase.mybc[2][0] = UNB; valCase.mybc[2][1] = UNB;

        // if(rank==0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_2);
        // if(rank==0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_4);
        // if(rank==0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
        // validation_3d(valCase,UP_SRHS,HEJ_6);
        if(rank==0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
        validation_3d(valCase,UP_SRHS,CHAT_2);

        // if(rank==0) printf("\n==============================  MIX UNBOUNDED TEST ============================================\n");
        // for (int id = 0; id < 3; id++) {
        //     for (int lr = 0; lr < 2; lr++) {
        //         valCase.mybc[0][0] = UNB;
        //         valCase.mybc[0][1] = UNB;
        //         valCase.mybc[1][0] = UNB;
        //         valCase.mybc[1][1] = UNB;
        //         valCase.mybc[2][0] = UNB;
        //         valCase.mybc[2][1] = UNB;
        //         // set to even
        //         valCase.mybc[id][lr] = EVEN;
        //         if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, HEJ_2);
        //         if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, HEJ_4);
        //         if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, HEJ_6);
        //         if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, CHAT_2);
        //         // set to odd
        //         valCase.mybc[id][lr] = ODD;
        //         if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, HEJ_2);
        //         if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, HEJ_4);
        //         if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, HEJ_6);
        //         if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
        //         validation_3d(valCase, UP_SRHS, CHAT_2);
        //         // reset to unbounded
        //         valCase.mybc[id][0] = UNB;
        //     }
        // }
    }
    MPI_Finalize();
}