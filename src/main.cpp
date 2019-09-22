/**
 * @file main.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief 
 * @version
 * @date 2019-09-22
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

#include <cmath>
#include <iostream>

#include "expint.hpp"
#include "tools.hpp"

#include "Solver.hpp"
#include "SwitchTopo.hpp"
#include "validation_3d.hpp"

#include "mpi.h"


int main(int argc, char *argv[]) {

    int rank;

    // MPI_Init(&argc, &argv);
    int provided;
    // set MPI_THREAD_MULTIPLE or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_SERIALIZED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested){
        FLUPS_ERROR("The MPI-provided thread behavior does not match");
    }

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
        valCase.mybc[0][0] = FLUPS::UNB; valCase.mybc[0][1] = FLUPS::UNB;
        valCase.mybc[1][0] = FLUPS::UNB; valCase.mybc[1][1] = FLUPS::UNB;
        valCase.mybc[2][0] = FLUPS::UNB; valCase.mybc[2][1] = FLUPS::UNB;

        if(rank==0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
        validation_3d(valCase,FLUPS::SRHS,FLUPS::HEJ_2);
        if(rank==0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
        validation_3d(valCase,FLUPS::SRHS,FLUPS::HEJ_4);
        if(rank==0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
        validation_3d(valCase,FLUPS::SRHS,FLUPS::HEJ_6);
        if(rank==0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
        validation_3d(valCase,FLUPS::SRHS,FLUPS::CHAT_2);

        if(rank==0) printf("\n==============================  MIX UNBOUNDED TEST ============================================\n");
        for (int id = 0; id < 3; id++) {
            for (int lr = 0; lr < 2; lr++) {
                valCase.mybc[0][0] = FLUPS::UNB;
                valCase.mybc[0][1] = FLUPS::UNB;
                valCase.mybc[1][0] = FLUPS::UNB;
                valCase.mybc[1][1] = FLUPS::UNB;
                valCase.mybc[2][0] = FLUPS::UNB;
                valCase.mybc[2][1] = FLUPS::UNB;
                // set to even
                valCase.mybc[id][lr] = FLUPS::EVEN;
                if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_2);
                if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_4);
                if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_6);
                if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::CHAT_2);
                // set to odd
                valCase.mybc[id][lr] = FLUPS::ODD;
                if (rank == 0) printf("------------------------------  HEJ_2  -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_2);
                if (rank == 0) printf("------------------------------  HEJ_4  -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_4);
                if (rank == 0) printf("------------------------------  HEJ_6  -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::HEJ_6);
                if (rank == 0) printf("------------------------------  CHAT_2 -----------------------------------------------------------------\n");
                validation_3d(valCase, FLUPS::SRHS, FLUPS::CHAT_2);
                // reset to unbounded
                valCase.mybc[id][0] = FLUPS::UNB;
            }
        }
    }
    MPI_Finalize();
}