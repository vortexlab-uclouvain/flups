/**
 * @file validation_3d.cpp
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

#include "validation_3d.hpp"
#include "omp.h"


using namespace std;

void validation_3d(const DomainDescr myCase, const FLUPS_SolverType type, const FLUPS_GreenType typeGreen) {
    validation_3d(myCase, type, typeGreen, 1);
}

/**
 * @brief computes the reference solution and the numerical one, outputs errors in a file
 * 
 * @param myCase description of the domain and initial condition
 * @param type type of solver
 * @param typeGreen type of Green function
 * @param nSolve number of times we call the same solver (for timing)
 */
void validation_3d(const DomainDescr myCase, const FLUPS_SolverType type, const FLUPS_GreenType typeGreen, const int nSolve) {
// void validation_3d(const DomainDescr myCase, const SolverType type, const GreenType typeGreen, const int nSolve) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int *   nglob  = myCase.nglob;
    const int *   nproc  = myCase.nproc;
    const double *L      = myCase.L;

    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};

    const FLUPS_BoundaryType mybc[3][2] = {myCase.mybc[0][0], myCase.mybc[0][1],
                                     myCase.mybc[1][0], myCase.mybc[1][1],
                                     myCase.mybc[2][0], myCase.mybc[2][1]};

    // create a real topology
    const FLUPS_Topology *topo    = flups_topo_new(0, nglob, nproc, false, NULL, FLUPS_ALIGNMENT);
    // const Topology *topo    = new Topology(0, nglob, nproc, false, NULL,FLUPS_ALIGNMENT);

    //-------------------------------------------------------------------------
    /** - Initialize the solver */
    //-------------------------------------------------------------------------
    std::string name = "validation_res" + std::to_string((int)(nglob[0]/L[0])) + "_nrank" + std::to_string(comm_size)+"_nthread" + std::to_string(omp_get_max_threads());
    FLUPS_Profiler* prof = flups_profiler_new_n(name.c_str());
    FLUPS_Solver *mysolver = flups_init_timed(topo, mybc, h, L,prof);
    flups_set_greenType(mysolver,typeGreen);
    flups_setup(mysolver);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    // FLUPS_INFO("topo memsize = %d vs %d %d %d",flups_topo_get_memsize(topo),topo->nmem(0),topo->nmem(1),topo->nmem(2));

    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *sol   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(rhs, 0, sizeof(double ) * flups_topo_get_memsize(topo));
    std::memset(sol, 0, sizeof(double ) * flups_topo_get_memsize(topo));
    std::memset(field, 0, sizeof(double ) * flups_topo_get_memsize(topo));

#ifndef MANUFACTURED_SOLUTION
    //-------------------------------------------------------------------------
    /** - fill the rhs and the solution */
    //-------------------------------------------------------------------------
    const double sigma     = 0.1;
    const double center[3] = {0.5, 0.5, 0.5};
    const double oosigma   = 1.0 / (sigma);
    const double oosigma2  = 1.0 / (sigma * sigma);
    const double oosigma3  = 1.0 / (sigma * sigma * sigma);

    int istart[3];
    flups_topo_get_istartGlob(topo,istart);

    /**
     * also accounting for various symmetry conditions. CAUTION: the solution for the Gaussian blob does not go to 0 fast enough
     * for `anal` to be used as a reference solution for cases where there is at least 1 symmetric (left AND right) or periodic direction
     */
    for (int j2 = -1; j2 < 2; j2++) {
        if (j2 != 0 && mybc[2][(j2 + 1) / 2] == UNB) continue;  //skip unbounded dirs
        for (int j1 = -1; j1 < 2; j1++) {
            if (j1 != 0 && mybc[1][(j1 + 1) / 2] == UNB) continue;  //skip unbounded dirs
            for (int j0 = -1; j0 < 2; j0++) {
                if (j0 != 0 && mybc[0][(j0 + 1) / 2] == UNB) continue;  //skip unbounded dirs

                double sign = 1.0;
                double centerPos[3];
                double orig[3] = {j0 * L[0], j1 * L[1], j2 * L[2]};  //inner left corner of the current block i'm looking at

                sign *= j0 == 0 ? 1.0 : 1 - 2 * (mybc[0][(j0 + 1) / 2] == ODD);  //multiply by -1 if the symm is odd
                sign *= j1 == 0 ? 1.0 : 1 - 2 * (mybc[1][(j1 + 1) / 2] == ODD);  //multiply by -1 if the symm is odd
                sign *= j2 == 0 ? 1.0 : 1 - 2 * (mybc[2][(j2 + 1) / 2] == ODD);  //multiply by -1 if the symm is odd

                centerPos[0] = orig[0] + ((j0 != 0) && (mybc[0][(j0 + 1) / 2] != PER) ? (1.0 - center[0]) * L[0] : (center[0] * L[0]));
                centerPos[1] = orig[1] + ((j1 != 0) && (mybc[1][(j1 + 1) / 2] != PER) ? (1.0 - center[1]) * L[1] : (center[1] * L[1]));
                centerPos[2] = orig[2] + ((j2 != 0) && (mybc[2][(j2 + 1) / 2] != PER) ? (1.0 - center[2]) * L[2] : (center[2] * L[2]));

                // printf("CENTER HERE IS: %d,%d,%d -- %lf,%lf,%lf -- %lf,%lf,%lf ++ %lf,%lf,%lf ** %lf\n",j0,j1,j2,orig[0],orig[1],orig[2],centerPos[0],centerPos[1],centerPos[2],\
                ( (j0!=0)&&(mybc[0][(j0+1)/2]!=PER )) ? (1.0-center[0])*L[0] : (center[0]*L[0]),\
                ( (j1!=0)&&(mybc[1][(j1+1)/2]!=PER )) ? (1.0-center[1])*L[1] : (center[1]*L[1]),\
                ( (j2!=0)&&(mybc[2][(j2+1)/2]!=PER )) ? (1.0-center[2])*L[2] : (center[2]*L[2]), sign);

                for (int i2 = 0; i2 < topo->nloc(2); i2++) {
                    for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                        for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                            double       x    = (istart[0] + i0 + 0.5) * h[0] - centerPos[0];
                            double       y    = (istart[1] + i1 + 0.5) * h[1] - centerPos[1];
                            double       z    = (istart[2] + i2 + 0.5) * h[2] - centerPos[2];
                            double       rho2 = (x * x + y * y + z * z) * oosigma2;
                            double       rho  = sqrt(rho2);
                            const size_t id    = flups_locID(0, i0, i1, i2, 0, nmem, 2);

                            // Gaussian
                            rhs[id] -= sign * c_1o4pi * oosigma3 * sqrt(2.0 / M_PI) * exp(-rho2 * 0.5);
                            sol[id] += sign * c_1o4pi * oosigma * 1.0 / rho * erf(rho * c_1osqrt2);
                        }
                    }
                }
            }
        }
    }

    // double lIs = 1.e10, 
    double gIs = 0.0;
    double lIs = 0.0;

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                const size_t id    = flups_locID(0, i0, i1, i2, 0, nmem, 2);
                lIs += sol[id];
                // lIs = min(sol[id],lIs);
            }
        }
    }
    // MPI_Allreduce(&lIs, &gIs, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lIs, &gIs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    gIs *= (h[0]*h[1]*h[2]);
    // const int nmem[3] = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};
    // for (int i2 = 0; i2 < topo->nloc(2); i2++) {
    //     for (int i1 = 0; i1 < topo->nloc(1); i1++) {
    //         for (int i0 = 0; i0 < topo->nloc(0); i0++) {
    //             const size_t id    = localIndex(0, i0, i1, i2, 0, nmem, 2);
    //             sol[id] -= gIs;
    //         }
    //     }
    // }
    printf("Integral sol : %lf\n",gIs);
#else
    //-------------------------------------------------------------------------
    /** - fill the rhs and the solution */
    //-------------------------------------------------------------------------

    int istart[3];
    flups_topo_get_istartGlob(topo,istart);

    {
        const int ax0     = flups_topo_get_axis(topo);
        const int ax1     = (ax0 + 1) % 3;
        const int ax2     = (ax0 + 2) % 3;
        const int nmem[3] = {flups_topo_get_nmem(topo,0),flups_topo_get_nmem(topo,1), flups_topo_get_nmem(topo,2)};
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    const size_t id = flups_locID(ax0, i0, i1, i2, ax0, nmem, 1);
                    sol[id]         = 1.0;
                }
            }
        }
    }

    manuF manuRHS[3] ;
    manuF manuSol[3] ;
    
    struct manuParams params[3]; 
    params[0].freq = 1;
    params[1].freq = 2;
    params[2].freq = 4;

    // Selecting manufactured solution compatible with the BCs
    for (int dir = 0; dir < 3; dir++) {
        if (mybc[dir][0] == PER && mybc[dir][1] == PER) {
            manuRHS[dir] = &d2dx2_fOddOdd;
            manuSol[dir] = &fOddOdd;
            if (params[dir].freq < 1) params[dir].freq = 1;
        } else if (mybc[dir][0] == ODD && mybc[dir][1] == ODD) {
            manuRHS[dir] = &d2dx2_fOddOdd;
            manuSol[dir] = &fOddOdd;
        } else if (mybc[dir][0] == EVEN && mybc[dir][1] == EVEN) {
            manuRHS[dir] = &d2dx2_fEvenEven;
            manuSol[dir] = &fEvenEven;
        }  else if (mybc[dir][0] == ODD && mybc[dir][1] == EVEN) {
            manuRHS[dir] = &d2dx2_fOddEven;
            manuSol[dir] = &fOddEven;
            if (params[dir].freq < 1) params[dir].freq = 1;
        }  else if (mybc[dir][0] == EVEN && mybc[dir][1] == ODD) {
            manuRHS[dir] = &d2dx2_fEvenOdd;
            manuSol[dir] = &fEvenOdd;
            if (params[dir].freq < 1) params[dir].freq = 1;            
        } else if (mybc[dir][0] == UNB) {
            params[dir].center  = .5;
            if (mybc[dir][1] == ODD) {
                params[dir].center  = .7;
                params[dir].sign[1] = -1.;
            } else if (mybc[dir][1] == EVEN) {
                params[dir].center  = .7;
                params[dir].sign[1] = +1.;
            }
            manuRHS[dir] = &d2dx2_fUnb;
            manuSol[dir] = &fUnb;
//still blobs missing if multiple mix
        } else if (mybc[dir][1] == UNB) {
            params[dir].center  = .5;
            if (mybc[dir][0] == ODD) {
                params[dir].center  = .3;
                params[dir].sign[0] = -1.;
            } else if (mybc[dir][0] == EVEN) {
                params[dir].center  = .3;
                params[dir].sign[0] = +1.;
            }
            manuRHS[dir] = &d2dx2_fUnb;
            manuSol[dir] = &fUnb;
        } else {
            // FLUPS_ERROR("I don''t know how to generate an analytical solution for this combination of BC.", LOCATION);
        }
    }

    {
        const int ax0     = flups_topo_get_axis(topo);
        const int ax1     = (ax0 + 1) % 3;
        const int ax2     = (ax0 + 2) % 3;
        const int nmem[3] = {flups_topo_get_nmem(topo,0),flups_topo_get_nmem(topo,1), flups_topo_get_nmem(topo,2)};
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    const size_t id = flups_locID(ax0, i0, i1, i2, ax0, nmem, 1);
                    const double x[3] = {(istart[ax0] + i0 + 0.5) * h[ax0],
                                         (istart[ax1] + i1 + 0.5) * h[ax1],
                                         (istart[ax2] + i2 + 0.5) * h[ax2]};

                    for (int dir = 0; dir < 3; dir++) {
                        const int dir2 = (dir + 1) % 3;
                        const int dir3 = (dir + 2) % 3;
                        sol[id] *= manuSol[dir](x[dir], L[dir], params[dir]);
                        rhs[id] += manuRHS[dir](x[dir], L[dir], params[dir]) * manuSol[dir2](x[dir2], L[dir2], params[dir2]) * manuSol[dir3](x[dir3], L[dir3], params[dir3]);
                    }
                }
            }
        }
    }

#endif



#ifdef DUMP_H5
    char msg[512];
    // write the source term and the solution
    sprintf(msg, "rhs_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0], mybc[0][1], mybc[1][0], mybc[1][1], mybc[2][0], mybc[2][1], nglob[0], nglob[1], nglob[2]);
    hdf5_dump(topo, msg, rhs);
    sprintf(msg, "anal_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0], mybc[0][1], mybc[1][0], mybc[1][1], mybc[2][0], mybc[2][1], nglob[0], nglob[1], nglob[2]);
    hdf5_dump(topo, msg, sol);
#endif

    //-------------------------------------------------------------------------
    /** - solve the equations */
    //-------------------------------------------------------------------------
    for(int is=0; is<nSolve; is++){
        flups_solve(mysolver,field, rhs, type);
    }

    flups_profiler_disp_root(prof,"solve");
    flups_profiler_free(prof);

    // lIs = 1.e10, gIs = 0.0;
    // for (int i2 = 0; i2 < topo->nloc(2); i2++) {
    //     for (int i1 = 0; i1 < topo->nloc(1); i1++) {
    //         for (int i0 = 0; i0 < topo->nloc(0); i0++) {
    //             const size_t id   = localindex_xyz(i0, i1, i2, topo);
    //             // lIs += rhs[id];
    //             lIs = min(rhs[id],lIs);
    //         }
    //     }
    // }
    // MPI_Allreduce(&lIs, &gIs, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    // // gIs *= (h[0]*h[1]*h[2]);
    // for (int i2 = 0; i2 < topo->nloc(2); i2++) {
    //     for (int i1 = 0; i1 < topo->nloc(1); i1++) {
    //         for (int i0 = 0; i0 < topo->nloc(0); i0++) {
    //             const size_t id   = localindex_xyz(i0, i1, i2, topo);
    //             rhs[id] -= gIs;
    //         }
    //     }
    // }

#ifdef DUMP_H5
    // write the source term and the solution
    sprintf(msg, "sol_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0], mybc[0][1], mybc[1][0], mybc[1][1], mybc[2][0], mybc[2][1], nglob[0], nglob[1], nglob[2]);
    hdf5_dump(topo, msg, rhs);
#endif    

    //-------------------------------------------------------------------------
    /** - compute the error */
    //-------------------------------------------------------------------------
    double lerr2 = 0.0;
    double lerri = 0.0;

    {
        const int ax0     = flups_topo_get_axis(topo);
        const int ax1     = (ax0 + 1) % 3;
        const int ax2     = (ax0 + 2) % 3;
        const int nmem[3] = {flups_topo_get_nmem(topo,0),flups_topo_get_nmem(topo,1), flups_topo_get_nmem(topo,2)};
        for (int i2 = 0; i2 < flups_topo_get_nloc(topo,ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo,ax1); i1++) {
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo,ax0); i0++) {
                    const size_t id = flups_locID(ax0, i0, i1, i2, ax0, nmem, 1);
                    const double err = sol[id] - field[id];

                    lerri = max(lerri, fabs(err));
                    lerr2 += (err * err) * h[0] * h[1] * h[2];
                }
            }
        }
    }
    double erri = 0.0;
    double err2 = 0.0;
    MPI_Allreduce(&lerr2, &err2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&lerri, &erri, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    err2 = sqrt(err2);

    char filename[512];
    string folder = "./data";

    sprintf(filename, "%s/%s_%d%d%d%d%d%d_typeGreen=%d.txt",folder.c_str(),__func__, mybc[0][0], mybc[0][1], mybc[1][0], mybc[1][1], mybc[2][0], mybc[2][1],typeGreen);

    if (rank == 0) {
        struct stat st = {0};
        if (stat(folder.c_str(), &st) == -1) {
                mkdir(folder.c_str(), 0770);
        }

        FILE *myfile = fopen(filename, "a+");
        if (myfile != NULL) {
            fprintf(myfile, "%d %12.12e %12.12e\n", nglob[0], err2, erri);
            fclose(myfile);
        } else {
            // FLUPS_CHECK(false, "unable to open file %s", filename, LOCATION);
        }
    }

    flups_free(sol);
    flups_free(rhs);
    flups_free(field);
    flups_cleanup(mysolver);
    flups_topo_free(topo);
}

