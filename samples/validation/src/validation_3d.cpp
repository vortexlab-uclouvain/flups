/**
 * @file validation_3d.cpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 *
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 *
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 *
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "validation_3d.hpp"

#include <limits.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <vector>

using namespace std;

void validation_3d(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int lda, const std::string output_dir) {
    validation_3d(myCase, typeGreen, 1, output_dir);
}

/**
 * @brief computes the reference solution and the numerical one, outputs errors in a file
 *
 * @param myCase description of the domain and initial condition
 * @param typeGreen type of Green function
 * @param lda leading dimension of array = number of vector components
 * @param nSolve number of times we call the same solver (for timing)
 */
void validation_3d(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int lda, const int nSolve, const std::string output_dir) {
    int      rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    const int    *nproc = myCase.nproc;
    const double *L     = myCase.L;

    const bool   is_cell = myCase.center_type[0] == CELL_CENTER;
    const double fact    = (double)(!is_cell);

    const int    nglob[3] = {myCase.nglob[0] + (int)fact, myCase.nglob[1] + (int)fact,
                          (myCase.nglob[2] == 1) ? 1 : myCase.nglob[2] + (int)fact};
    const double h[3]     = {L[0] / (nglob[0] - fact), L[1] / (nglob[1] - fact),
                         nglob[2] == 1 ? 1 : L[2] / (nglob[2] - fact)};

    FLUPS_CenterType center_type[3];
    for (int i = 0; i < 3; i++) {
        center_type[i] = myCase.center_type[i];
    }

    FLUPS_BoundaryType *mybc[3][2];
    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            mybc[id][is] = (FLUPS_BoundaryType *)flups_malloc(sizeof(int) * lda);
        }
    }

    if (myCase.dovectorbc) {
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                for (int lia = 0; lia < lda; lia++) {
                    mybc[id][is][lia] = myCase.mybcv[id][is][lia];
                }
            }
        }
    } else {
        for (int id = 0; id < 3; id++) {
            for (int is = 0; is < 2; is++) {
                for (int lia = 0; lia < lda; lia++) {
                    mybc[id][is][lia] = myCase.mybc[id][is];
                }
            }
        }
    }

    // create a real topology
    FLUPS_Topology *topo = flups_topo_new(0, lda, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);
    // const Topology *topo    = new Topology(0, 1, nglob, nproc, false, NULL,FLUPS_ALIGNMENT,comm);

    //-------------------------------------------------------------------------
    /** - Initialize the solver */
    //-------------------------------------------------------------------------
    std::string     name = "validation_nx" + std::to_string((int)(nglob[0])) + "_ny" + std::to_string((int)(nglob[1])) + "_nz" + std::to_string((int)(nglob[2])) + "_nrank" + std::to_string(comm_size);
    FLUPS_Profiler *prof = flups_profiler_new_n(name.c_str());

    m_profStart(prof, "Validation--Init-Flups");
    FLUPS_Solver *mysolver;
    mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, prof);

    flups_set_greenType(mysolver, typeGreen);
    flups_setup(mysolver, true);

    // update the comm and the rank
    comm = flups_topo_get_comm(topo);
    MPI_Comm_rank(comm, &rank);
    m_profStop(prof, "Validation--Init-Flups");

    // flups_switchtopo_info(mysolver);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    m_profStart(prof, "Validation--Init-RHS");
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *sol   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(sol, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

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
    flups_topo_get_istartGlob(topo, istart);

    /**
     * also accounting for various symmetry conditions. CAUTION: the solution for the Gaussian blob does not go to 0 fast enough
     * for `anal` to be used as a reference solution for cases where there is at least 1 symmetric (left AND right) or periodic direction
     */
    for (int lia = 0; lia < lda; lia++) {
        for (int j2 = -1; j2 < 2; j2++) {
            if (j2 != 0 && mybc[2][(j2 + 1) / 2] == UNB) continue;  // skip unbounded dirs
            for (int j1 = -1; j1 < 2; j1++) {
                if (j1 != 0 && mybc[1][(j1 + 1) / 2] == UNB) continue;  // skip unbounded dirs
                for (int j0 = -1; j0 < 2; j0++) {
                    if (j0 != 0 && mybc[0][(j0 + 1) / 2] == UNB) continue;  // skip unbounded dirs

                    double sign = 1.0;
                    double centerPos[3];
                    double orig[3] = {j0 * L[0], j1 * L[1], j2 * L[2]};  // inner left corner of the current block i'm looking at

                    sign *= j0 == 0 ? 1.0 : 1 - 2 * (mybc[0][(j0 + 1) / 2] == ODD);  // multiply by -1 if the symm is odd
                    sign *= j1 == 0 ? 1.0 : 1 - 2 * (mybc[1][(j1 + 1) / 2] == ODD);  // multiply by -1 if the symm is odd
                    sign *= j2 == 0 ? 1.0 : 1 - 2 * (mybc[2][(j2 + 1) / 2] == ODD);  // multiply by -1 if the symm is odd

                    centerPos[0] = orig[0] + ((j0 != 0) && (mybc[0][(j0 + 1) / 2] != PER) ? (1.0 - center[0]) * L[0] : (center[0] * L[0]));
                    centerPos[1] = orig[1] + ((j1 != 0) && (mybc[1][(j1 + 1) / 2] != PER) ? (1.0 - center[1]) * L[1] : (center[1] * L[1]));
                    centerPos[2] = orig[2] + ((j2 != 0) && (mybc[2][(j2 + 1) / 2] != PER) ? (1.0 - center[2]) * L[2] : (center[2] * L[2]));

                    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
                        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                                const double shift = 0.5;
                                double       x     = (istart[0] + i0 + shift) * h[0] - centerPos[0];
                                double       y     = (istart[1] + i1 + shift) * h[1] - centerPos[1];
                                double       z     = (istart[2] + i2 + shift) * h[2] - centerPos[2];
                                double       rho2  = (x * x + y * y + z * z) * oosigma2;
                                double       rho   = sqrt(rho2);
                                const size_t id    = flups_locID(0, i0, i1, i2, lia, 0, nmem, 2);

                                // Gaussian
                                rhs[id] -= sign * c_1o4pi * oosigma3 * sqrt(2.0 / M_PI) * exp(-rho2 * 0.5);
                                sol[id] += sign * c_1o4pi * oosigma * 1.0 / rho * erf(rho * c_1osqrt2);
                            }
                        }
                    }
                }
            }
        }
    }

    double *gIs = (double *)malloc(lda * sizeof(double));
    double *lIs = (double *)malloc(lda * sizeof(double));

    std::memset(gIs, 0, sizeof(double) * lda);
    std::memset(lIs, 0, sizeof(double) * lda);

    for (int lia = 0; lia < lda; lia++) {
        for (int i2 = 0; i2 < topo->nloc(2); i2++) {
            for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                    const size_t id = flups_locID(0, i0, i1, i2, lia, 0, nmem, 2);
                    lIs[lia] += sol[id];
                    // lIs = min(sol[id],lIs);
                }
            }
        }
    }
    // MPI_Allreduce(&lIs, &gIs, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(lIs, gIs, lda, MPI_DOUBLE, MPI_SUM, comm);

    printf("Integral sol : ");
    for (int lia = 0; lia < lda; lia++) {
        gIs[lia] *= (h[0] * h[1] * h[2]);
        printf("%lf ", gIs[lia]);
    }
    printf("\n");

#else
    //-------------------------------------------------------------------------
    /** - fill the rhs and the solution */
    //-------------------------------------------------------------------------

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);

    {
        const int ax0     = flups_topo_get_axis(topo);
        const int ax1     = (ax0 + 1) % 3;
        const int ax2     = (ax0 + 2) % 3;
        const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};
        for (int lia = 0; lia < lda; lia++) {
            for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
                for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                    for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                        const size_t id = flups_locID(ax0, i0, i1, i2, lia, ax0, nmem, 1);
                        sol[id]         = 1.0;
                    }
                }
            }
        }
    }

    manuF manuRHS[3];
    manuF manuSol[3];
    manuF manuDer[3];

    for (int lia = 0; lia < lda; lia++) {
        // Selecting manufactured solution compatible with the BCs
        struct manuParams params[3];
        params[0].freq = 1;
        params[1].freq = 2;
        params[2].freq = 4;
        for (int dir = 0; dir < 3; dir++) {
            if (mybc[dir][0][lia] == PER && mybc[dir][1][lia] == PER) {
                manuRHS[dir] = &d2dx2_fOddOdd;
                manuSol[dir] = &fOddOdd;
                if (params[dir].freq < 1) params[dir].freq = 1;
            } else if (mybc[dir][0][lia] == ODD && mybc[dir][1][lia] == ODD) {
                manuRHS[dir] = &d2dx2_fOddOdd;
                manuSol[dir] = &fOddOdd;
            } else if (mybc[dir][0][lia] == EVEN && mybc[dir][1][lia] == EVEN) {
                manuRHS[dir] = &d2dx2_fEvenEven;
                manuSol[dir] = &fEvenEven;
            } else if (mybc[dir][0][lia] == ODD && mybc[dir][1][lia] == EVEN) {
                manuRHS[dir] = &d2dx2_fOddEven;
                manuSol[dir] = &fOddEven;
                if (params[dir].freq < 1) params[dir].freq = 1;
            } else if (mybc[dir][0][lia] == EVEN && mybc[dir][1][lia] == ODD) {
                manuRHS[dir] = &d2dx2_fEvenOdd;
                manuSol[dir] = &fEvenOdd;
                if (params[dir].freq < 1) params[dir].freq = 1;
            } else if (mybc[dir][0][lia] == UNB) {
                if (mybc[dir][1][lia] == ODD) {
                    params[dir].center  = .7;
                    params[dir].sign[1] = -1.;
                } else if (mybc[dir][1][lia] == EVEN) {
                    params[dir].center  = .7;
                    params[dir].sign[1] = +1.;
                }
                // manuRHS[dir] = &d2dx2_fUnb;
                // manuSol[dir] = &fUnb;
                manuRHS[dir] = &d2dx2_fUnbSpietz;
                manuSol[dir] = &fUnbSpietz;
            } else if (mybc[dir][1][lia] == UNB) {
                if (mybc[dir][0][lia] == ODD) {
                    params[dir].center  = .3;
                    params[dir].sign[0] = -1.;
                } else if (mybc[dir][0][lia] == EVEN) {
                    params[dir].center  = .3;
                    params[dir].sign[0] = +1.;
                }
                // manuRHS[dir] = &d2dx2_fUnb;
                // manuSol[dir] = &fUnb;
                manuRHS[dir] = &d2dx2_fUnbSpietz;
                manuSol[dir] = &fUnbSpietz;
            } else {
                manuRHS[dir] = &fZero;
                manuSol[dir] = &fCst;
                // FLUPS_ERROR("I don''t know how to generate an analytical solution for this combination of BC.");
            }
        }

        // USE THE FOLLOWING TO TEST THE K=0 PART OF THE 1DIRUNBOUNDED KERNEL
        //  manuRHS[0] = &fZero;
        //  manuSol[0] = &fCst;
        //  manuRHS[1] = &fZero;
        //  manuSol[1] = &fCst;
        //  manuRHS[2] = &d2dx2_fUnbSpietz;
        //  manuSol[2] = &fUnbSpietz;

        {
            // Obtaining the reference sol and rhs
            const int ax0     = flups_topo_get_axis(topo);
            const int ax1     = (ax0 + 1) % 3;
            const int ax2     = (ax0 + 2) % 3;
            const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};

            // printf("for dim %d, we use the sign = %f %f %f\n",lia,params[0].sign[0],params[1].sign[0],params[2].sign[0]);

            for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
                for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                    for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                        const size_t id    = flups_locID(ax0, i0, i1, i2, lia, ax0, nmem, 1);
                        const double shift = is_cell ? 0.5 : 0.0;
                        const double x[3]  = {(istart[ax0] + i0 + shift) * h[ax0],
                                              (istart[ax1] + i1 + shift) * h[ax1],
                                              (istart[ax2] + i2 + shift) * h[ax2]};

                        // const double y[3] = {0.,(x[1]-.5)/params.sigma[1],(x[2]-.5)/params.sigma[2] };
                        // const double rsq = y[1]*y[1] + y[2]*y[2] ; //((x[1]-.5)*(x[1]-.5)+(x[2]-.5)*(x[2]-.5)) / .25;
                        // const double r = sqrt(rsq);

                        // CST(z):
                        //  sol[id] = fabs(rsq)>=1. ?  0.0 : exp(c_C * (1. - 1. / (1. - rsq))) ;
                        //  rhs[id] = fab s(rsq) >= 1. ? 0.0 : 4.*c_C* exp(c_C * (1. - 1. / (1. - rsq))) / (pow(rsq - 1., 4) * params.sigma[id] * params.sigma[id]) * \
                        //              (c_C * rsq - 1. + pow(y[1],4) + pow(y[2],4) + 2.* y[1]*y[1]*y[2]*y[2]) ;

                        // SIN(z):
                        //  sol[id] = fabs(rsq) >= 1. ? 0.0 : exp(c_C * (1. - 1. / (1. - rsq))) * (sin(2. * M_PI * x[0] / L[0]));
                        //  rhs[id] = fab s(rsq) >= 1. ? 0.0 : 4.*c_C* exp(c_C * (1. - 1. / (1. - rsq))) / (pow(rsq - 1., 4) * params.sigma[id] * params.sigma[id]) * \
                        //              (c_C * rsq - 1. + pow(y[1],4) + pow(y[2],4) + 2.* y[1]*y[1]*y[2]*y[2]) * (sin(2.*M_PI*x[0] /L[0])) ;
                        //  rhs[id] += fabs(rsq) >= 1. ? 0.0 : -sin(2*M_PI*x[0] /L[0]) * (2. * M_PI / L[0])* (2. * M_PI / L[0])  * exp(c_C * (1. - 1. / (1. - rsq))) ;

                        // CST(z)+sin(z):
                        //  sol[id] = fabs(rsq) >= 1. ? 0.0 : exp(c_C * (1. - 1. / (1. - rsq))) * (1. + sin(2. * M_PI * x[0] / L[0]));
                        //  rhs[id] = fab s(rsq) >= 1. ? 0.0 : 4.*c_C* exp(c_C * (1. - 1. / (1. - rsq))) / (pow(rsq - 1., 4) * params.sigma[id] * params.sigma[id]) * \
                        //              (c_C * rsq - 1. + pow(y[1],4) + pow(y[2],4) + 2.* y[1]*y[1]*y[2]*y[2]) * (1.+sin(2.*M_PI*x[0] /L[0])) ;
                        //  rhs[id] += fabs(rsq) >= 1. ? 0.0 : -sin(2*M_PI*x[0] /L[0]) * (2. * M_PI / L[0])* (2. * M_PI / L[0])  * exp(c_C * (1. - 1. / (1. - rsq))) ;

                        for (int dir = 0; dir < 3; dir++) {
                            const int dir2 = (dir + 1) % 3;
                            const int dir3 = (dir + 2) % 3;
                            // if we do the RHS type, fill the solution now
                            sol[id] *= manuSol[dir](x[dir], L[dir], params[dir]);
                            rhs[id] += manuRHS[dir](x[dir], L[dir], params[dir]) * manuSol[dir2](x[dir2], L[dir2], params[dir2]) * manuSol[dir3](x[dir3], L[dir3], params[dir3]);
                        }
                    }
                }
            }
        }
    }

#endif
    m_profStop(prof, "Validation--Init-RHS");

#ifdef DUMP_DBG
    char msg[512];
    // write the source term and the solution
    sprintf(msg, "rhs_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0][0], mybc[0][1][0], mybc[1][0][0], mybc[1][1][0], mybc[2][0][0], mybc[2][1][0], nglob[0], nglob[1], nglob[2]);
    flups_hdf5_dump(topo, msg, rhs);
    sprintf(msg, "anal_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0][0], mybc[0][1][0], mybc[1][0][0], mybc[1][1][0], mybc[2][0][0], mybc[2][1][0], nglob[0], nglob[1], nglob[2]);
    flups_hdf5_dump(topo, msg, sol);
#endif

    //-------------------------------------------------------------------------
    /** - Warm up */
    //-------------------------------------------------------------------------
    for (int is = 0; is < 5; is++) {
        MPI_Barrier(MPI_COMM_WORLD);
        m_profStart(prof, "Validation--Warm-up");
        flups_solve(mysolver, field, rhs, STD);
        m_profStop(prof, "Validation--Warm-up");
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //-------------------------------------------------------------------------
    /** - solve the equations */
    //-------------------------------------------------------------------------
    for (int is = 0; is < nSolve; is++) {
        MPI_Barrier(MPI_COMM_WORLD);
        m_profStart(prof, "Validation--Solve");
        flups_solve(mysolver, field, rhs, STD);
        m_profStop(prof, "Validation--Solve");
        MPI_Barrier(MPI_COMM_WORLD);
    }

    flups_profiler_disp(prof);

#ifdef DUMP_DBG
    // write the source term and the solution
    sprintf(msg, "sol_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0][0], mybc[0][1][0], mybc[1][0][0], mybc[1][1][0], mybc[2][0][0], mybc[2][1][0], nglob[0], nglob[1], nglob[2]);
    flups_hdf5_dump(topo, msg, field);
#endif

    //-------------------------------------------------------------------------
    /** - compute the error */
    //-------------------------------------------------------------------------
    double *lerr2 = (double *)malloc(lda * sizeof(double));
    double *lerri = (double *)malloc(lda * sizeof(double));
    double *err2  = (double *)malloc(lda * sizeof(double));
    double *erri  = (double *)malloc(lda * sizeof(double));

    std::memset(lerr2, 0, sizeof(double) * lda);
    std::memset(lerri, 0, sizeof(double) * lda);
    std::memset(err2, 0, sizeof(double) * lda);
    std::memset(erri, 0, sizeof(double) * lda);

    // determine the volume associated to a mesh
    double vol = 1.0;
    for (int id = 0; id < 3; id++) {
        if (mybc[id][0][0] != NONE && mybc[id][1][0] != NONE) {
            vol *= h[id];
        }
    }

    {
        const int ax0     = flups_topo_get_axis(topo);
        const int ax1     = (ax0 + 1) % 3;
        const int ax2     = (ax0 + 2) % 3;
        const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};
        for (int lia = 0; lia < lda; lia++) {
            for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
                for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                    for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                        const size_t id  = flups_locID(ax0, i0, i1, i2, lia, ax0, nmem, 1);
                        const double err = sol[id] - field[id];

                        lerri[lia] = max(lerri[lia], fabs(err));
                        lerr2[lia] += (err * err) * vol;
                    }
                }
            }
        }
    }
    MPI_Allreduce(lerr2, err2, lda, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(lerri, erri, lda, MPI_DOUBLE, MPI_MAX, comm);

    for (int i = 0; i < lda; i++) {
        err2[i] = sqrt(err2[i]);
    }

    string folder  = output_dir + "/data";
    string ct_name = is_cell ? "CellCenter" : "NodeCenter";
    char   filename[PATH_MAX];
    sprintf(filename, "%s/%s_%s_%d%d%d%d%d%d_typeGreen=%d.txt", folder.c_str(), __func__, ct_name.c_str(), mybc[0][0][0], mybc[0][1][0], mybc[1][0][0], mybc[1][1][0], mybc[2][0][0], mybc[2][1][0], typeGreen);

    if (rank == 0) {
        struct stat st_output = {0};
        struct stat st_folder = {0};
        if (stat(output_dir.c_str(), &st_output) == -1) {
            mkdir(output_dir.c_str(), 0770);
        }

        if (stat(folder.c_str(), &st_folder) == -1) {
            mkdir(folder.c_str(), 0770);
        }

        FILE *myfile = fopen(filename, "a+");
        if (myfile != NULL) {
            printf("Opening file %s !\n", filename);
            fprintf(myfile, "%d ", nglob[0]);
            for (int i = 0; i < lda; i++) {
                fprintf(myfile, "%12.12e %12.12e ", err2[i], erri[i]);
            }
            fprintf(myfile, "\n");

            fclose(myfile);
        } else {
            printf("unable to open file %s ! Here is what I would have written:", filename);
            printf("%d ", nglob[0]);
            for (int i = 0; i < lda; i++) {
                printf("%12.12e %12.12e ", err2[i], erri[i]);
            }
            printf("\n");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(lerr2);
    free(lerri);
    free(err2);
    free(erri);

    flups_free(sol);
    flups_free(rhs);
    flups_free(field);
    flups_cleanup_solver(mysolver);
    flups_cleanup_backend();
    flups_profiler_free(prof);
    flups_topo_free(topo);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }
}

//==============================================================================
// compact-source pad check (--source=1 --padcheck=P)
//
// With 2 unbounded + 1 spectral directions, every k_z != 0 mode uses the Green's function of a
// *screened* two-dimensional lattice operator, built by quadrature (lgf_2unb1spe.hpp). An inexact
// kernel, e.g. its symbol sampled on the doubled grid, is a 2N-periodised kernel: each k_z != 0 mode
// sees images of the source at a distance d ~ L in the unbounded directions, a relative error
// ~ exp(-|k_z| d), largest for the smallest k_z, i.e. for a long spectral direction. That error does
// not decrease with h, and the manufactured-solution test of validation_3d is blind to it (see
// validation_3d.hpp). With 1 or 3 unbounded directions the kernel is a closed form or a table in
// physical space and has no such error.
//
// The source is white noise, uniform in [-1, 1] on the reference box and zero outside: it is compact
// and excites every mode of the spectral direction, hence every slice of the kernel. Its free-space
// solution does not depend on how far the box extends, so a solve on a box P times larger in the two
// unbounded directions (same h, same origin, same source) gives the same values on the common box.
// Both solves share the discretisation error, which cancels in the difference: max|dU| is the image
// error alone, round-off for an exact kernel. The difference is also reported mode by mode along the
// spectral direction (padnoise_mode lines). The images are strongest for L_z > L_x.
//==============================================================================

/**
 * @brief white noise, uniform in [-1, 1], as a hash (splitmix64) of the point index @p i
 *
 * A hash rather than a random sequence: the value at a point depends neither on the decomposition
 * nor on the size of the box it is solved in.
 */
static inline double padcheck_noise_(const uint64_t i) {
    uint64_t z = i + 0x9e3779b97f4a7c15ULL;
    z          = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z          = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z          = z ^ (z >> 31);
    return 2.0 * ((double)(z >> 11) / 9007199254740992.0) - 1.0;
}

/**
 * @brief solve the compact source on a box padded by @p pad in the two unbounded directions
 *
 * The origin, the grid spacing and the physical location of the source are the same for every
 * @p pad, so the first `nc[d]` points of each direction of the padded run coincide with the points
 * of the reference run. The solution restricted to that common box is accumulated (globally, on
 * every rank) into @p ucommon, in the layout `i0 + nc0 * (i1 + nc1 * i2)` with the indices given in
 * *dimension* order.
 */
static void padcheck_solve_(const DomainDescr &myCase, const FLUPS_GreenType typeGreen, const int pad,
                            const int dir_u[2], const int dir_s, const int nc[3], double *ucommon, MPI_Comm comm) {
    const bool   is_cell = myCase.center_type[0] == CELL_CENTER;
    const double fact    = (double)(!is_cell);

    // the padded domain: same h, same origin, pad times larger in the two unbounded directions
    int    ncell[3];
    double L[3];
    for (int d = 0; d < 3; ++d) {
        const bool isunb = (d == dir_u[0] || d == dir_u[1]);
        ncell[d]         = isunb ? myCase.nglob[d] * pad : myCase.nglob[d];
        L[d]             = isunb ? myCase.L[d] * pad : myCase.L[d];
    }
    const int    nglob[3] = {ncell[0] + (int)fact, ncell[1] + (int)fact, ncell[2] + (int)fact};
    const double h[3]     = {L[0] / ncell[0], L[1] / ncell[1], L[2] / ncell[2]};

    // the source: noise on the *reference* box, indexed in its layout; the last node of a periodic
    // direction is the first one
    const bool zper  = myCase.mybc[dir_s][0] == PER;
    const int  nzsrc = (zper && !is_cell) ? nc[dir_s] - 1 : nc[dir_s];

    FLUPS_CenterType center_type[3];
    for (int d = 0; d < 3; ++d) {
        center_type[d] = myCase.center_type[d];
    }
    FLUPS_BoundaryType *mybc[3][2];
    for (int d = 0; d < 3; ++d) {
        for (int s = 0; s < 2; ++s) {
            mybc[d][s]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int));
            mybc[d][s][0] = myCase.mybc[d][s];
        }
    }

    FLUPS_Topology *topo     = flups_topo_new(0, 1, nglob, myCase.nproc, false, NULL, FLUPS_ALIGNMENT, comm);
    FLUPS_Profiler *prof     = flups_profiler_new_n("padcheck");
    FLUPS_Solver   *mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, prof);
    flups_set_greenType(mysolver, typeGreen);
    flups_setup(mysolver, true);

    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                int          g[3];
                g[ax0]            = istart[ax0] + i0;
                g[ax1]            = istart[ax1] + i1;
                g[ax2]            = istart[ax2] + i2;
                if (g[dir_u[0]] < nc[dir_u[0]] && g[dir_u[1]] < nc[dir_u[1]]) {
                    int gr[3] = {g[0], g[1], g[2]};
                    gr[dir_s] = g[dir_s] % nzsrc;
                    rhs[id]   = padcheck_noise_(gr[0] + nc[0] * (gr[1] + ((uint64_t)nc[1]) * gr[2]));
                }
            }
        }
    }

    flups_solve(mysolver, field, rhs, STD);

    // gather the common box
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                int          g[3];
                g[ax0] = istart[ax0] + i0;
                g[ax1] = istart[ax1] + i1;
                g[ax2] = istart[ax2] + i2;
                if (g[0] < nc[0] && g[1] < nc[1] && g[2] < nc[2]) {
                    ucommon[g[0] + nc[0] * (g[1] + ((size_t)nc[1]) * g[2])] = field[id];
                }
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, ucommon, nc[0] * nc[1] * nc[2], MPI_DOUBLE, MPI_SUM, comm);

    flups_free(rhs);
    flups_free(field);
    flups_cleanup_solver(mysolver);
    flups_cleanup_backend();
    flups_profiler_free(prof);
    flups_topo_free(topo);
    for (int d = 0; d < 3; ++d) {
        for (int s = 0; s < 2; ++s) {
            flups_free(mybc[d][s]);
        }
    }
}

void padcheck_3d(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int pad, const std::string output_dir) {
    int      rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    //-------------------------------------------------------------------------
    /** - identify the two unbounded directions and the spectral one */
    //-------------------------------------------------------------------------
    int dir_u[2] = {-1, -1};
    int dir_s    = -1;
    int n_unb    = 0;
    for (int d = 0; d < 3; ++d) {
        if (myCase.mybc[d][0] == UNB && myCase.mybc[d][1] == UNB) {
            m_assert(n_unb < 2, "the pad check needs exactly 2 unbounded directions");
            dir_u[n_unb] = d;
            n_unb += 1;
        } else {
            m_assert(dir_s < 0, "the pad check needs exactly 1 non-unbounded direction");
            m_assert((myCase.mybc[d][0] == PER && myCase.mybc[d][1] == PER) ||
                         ((myCase.mybc[d][0] == EVEN || myCase.mybc[d][0] == ODD) && (myCase.mybc[d][1] == EVEN || myCase.mybc[d][1] == ODD)),
                     "the non-unbounded direction of the pad check must be PER, or EVEN/ODD on both sides (got %d %d)", myCase.mybc[d][0], myCase.mybc[d][1]);
            dir_s = d;
        }
    }
    m_assert(n_unb == 2 && dir_s >= 0, "the pad check needs 2 UNB-UNB directions and 1 spectral one");

    const bool is_cell = myCase.center_type[0] == CELL_CENTER;
    const int  fact    = (int)(!is_cell);
    const int  nc[3]   = {myCase.nglob[0] + fact, myCase.nglob[1] + fact, myCase.nglob[2] + fact};
    const size_t ncom  = ((size_t)nc[0]) * nc[1] * nc[2];

    double *u1 = (double *)flups_malloc(sizeof(double) * ncom);
    double *up = (double *)flups_malloc(sizeof(double) * ncom);
    std::memset(u1, 0, sizeof(double) * ncom);
    std::memset(up, 0, sizeof(double) * ncom);

    padcheck_solve_(myCase, typeGreen, 1, dir_u, dir_s, nc, u1, comm);
    if (pad > 1) {
        padcheck_solve_(myCase, typeGreen, pad, dir_u, dir_s, nc, up, comm);
    } else {
        std::memcpy(up, u1, sizeof(double) * ncom);
    }

    //-------------------------------------------------------------------------
    /** - the modes of the spectral direction: the transform FLUPS uses there, on the points it takes */
    //-------------------------------------------------------------------------
    // mode j at the point k: w_k trig(dk (j + oj) (k + ok)) for k in [klo, khi), w_k = wlo / whi on the
    // first / last point, 1 elsewhere
    const FLUPS_BoundaryType bl     = myCase.mybc[dir_s][0];
    const FLUPS_BoundaryType br     = myCase.mybc[dir_s][1];
    const int                ns     = nc[dir_s];
    const bool               is_sin = (bl == ODD);

    int    klo = 0, khi = ns, nj;
    double dk, oj = 0.0, ok = is_cell ? 0.5 : 0.0, wlo = 1.0, whi = 1.0;
    bool   is_dft = false;
    if (bl == PER) {
        // DFT, the last node is the first one
        khi    = is_cell ? ns : ns - 1;
        nj     = khi / 2 + 1;
        dk     = c_2pi / khi;
        ok     = 0.0;
        is_dft = true;
    } else if (is_cell) {
        // DCT-II (EVEN-EVEN), DST-II (ODD-ODD), DCT-IV (EVEN-ODD), DST-IV (ODD-EVEN)
        nj = ns;
        dk = M_PI / ns;
        oj = (bl != br) ? 0.5 : ((bl == ODD) ? 1.0 : 0.0);
    } else {
        // on n = ns - 1 cells, without the zero ODD end points: DCT-I (EVEN-EVEN), DST-I (ODD-ODD),
        // DCT-III (EVEN-ODD), DST-III (ODD-EVEN)
        const int n = ns - 1;
        dk          = M_PI / n;
        klo         = (bl == ODD) ? 1 : 0;
        khi         = (br == ODD) ? n : ns;
        if (bl == br) {
            nj  = (bl == EVEN) ? n + 1 : n - 1;
            oj  = (bl == ODD) ? 1.0 : 0.0;
            wlo = (bl == EVEN) ? 0.5 : 1.0;
            whi = wlo;
        } else {
            nj  = n;
            oj  = 0.5;
            wlo = (bl == EVEN) ? 0.5 : 1.0;
            whi = (br == EVEN) ? 0.5 : 1.0;
        }
    }
    const int      nk = khi - klo;
    vector<double> bre((size_t)nj * nk), bim(is_dft ? (size_t)nj * nk : 0);
    for (int j = 0; j < nj; ++j) {
        for (int kk = 0; kk < nk; ++kk) {
            const double w           = ((kk == 0) ? wlo : 1.0) * ((kk == nk - 1) ? whi : 1.0);
            const double ph          = dk * (j + oj) * (klo + kk + ok);
            bre[(size_t)j * nk + kk] = w * (is_sin ? sin(ph) : cos(ph));
            if (is_dft) {
                bim[(size_t)j * nk + kk] = -w * sin(ph);
            }
        }
    }

    //-------------------------------------------------------------------------
    /** - the difference and the solution: max over the box, and of each mode over the (unbounded) plane */
    //-------------------------------------------------------------------------
    // red = [max_p |dU_j(p)| for j < nj, max_p |U_j(p)| for j < nj, max|dU|, max|U|], the planar points
    // p split over the ranks
    const int      nu0 = nc[dir_u[0]];
    const int      nu1 = nc[dir_u[1]];
    vector<double> red(2 * nj + 2, 0.0);
    vector<double> dz(nk), uz(nk);
    for (int p = rank; p < nu0 * nu1; p += comm_size) {
        int g[3];
        g[dir_u[0]] = p % nu0;
        g[dir_u[1]] = p / nu0;
        for (int kk = 0; kk < nk; ++kk) {
            g[dir_s]         = klo + kk;
            const size_t idx = g[0] + nc[0] * (g[1] + ((size_t)nc[1]) * g[2]);
            uz[kk]           = u1[idx];
            dz[kk]           = u1[idx] - up[idx];
            red[2 * nj]      = max(red[2 * nj], fabs(dz[kk]));
            red[2 * nj + 1]  = max(red[2 * nj + 1], fabs(uz[kk]));
        }
        for (int j = 0; j < nj; ++j) {
            const double *br_j = bre.data() + (size_t)j * nk;
            double        dre = 0.0, ure = 0.0, dim = 0.0, uim = 0.0;
            for (int kk = 0; kk < nk; ++kk) {
                dre += br_j[kk] * dz[kk];
                ure += br_j[kk] * uz[kk];
            }
            if (is_dft) {
                const double *bi_j = bim.data() + (size_t)j * nk;
                for (int kk = 0; kk < nk; ++kk) {
                    dim += bi_j[kk] * dz[kk];
                    uim += bi_j[kk] * uz[kk];
                }
            }
            red[j]      = max(red[j], sqrt(dre * dre + dim * dim));
            red[nj + j] = max(red[nj + j], sqrt(ure * ure + uim * uim));
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, red.data(), 2 * nj + 2, MPI_DOUBLE, MPI_MAX, comm);
    const double dmax   = red[2 * nj];
    const double umax   = red[2 * nj + 1];
    double       dUjmax = 0.0, Ujmax = 0.0, dUorel = 0.0;
    for (int j = 0; j < nj; ++j) {
        dUjmax = max(dUjmax, red[j]);
        Ujmax  = max(Ujmax, red[nj + j]);
        if (red[nj + j] > 0.0) {
            dUorel = max(dUorel, red[j] / red[nj + j]);
        }
    }

    const double LzoLx = myCase.L[dir_s] / myCase.L[dir_u[0]];
    if (rank == 0) {
        printf("padcheck: N=%d,%d,%d Lz/Lx=%g pad=%d kernel=%d max|dU|=%12.12e max|U|=%12.12e maxj(dU_j)/maxj(U_j)=%12.12e maxj(dU_j/U_j)=%12.12e\n",
               nc[0], nc[1], nc[2], LzoLx, pad, (int)typeGreen, dmax, umax, dUjmax / Ujmax, dUorel);
        // res and h2Rf, the lattice residual L u - h^2 R f and its scale, are not computed (nan)
        for (int j = 0; j < nj; ++j) {
            printf("padnoise_mode: kernel=%d Lz/Lx=%g pad=%d j=%d dU=%12.12e U=%12.12e res=nan h2Rf=nan\n",
                   (int)typeGreen, LzoLx, pad, j, red[j], red[nj + j]);
        }

        string folder  = output_dir + "/data";
        string ct_name = is_cell ? "CellCenter" : "NodeCenter";
        char   filename[PATH_MAX];
        // NOTE: deliberately NOT the file name used by validation_3d, so that the data_ref
        // regression of scripts/check_res.py is never confused by these lines
        sprintf(filename, "%s/padcheck_%s_%d%d%d%d%d%d_typeGreen=%d.txt", folder.c_str(), ct_name.c_str(),
                myCase.mybc[0][0], myCase.mybc[0][1], myCase.mybc[1][0], myCase.mybc[1][1], myCase.mybc[2][0], myCase.mybc[2][1], typeGreen);

        struct stat st_output = {0};
        struct stat st_folder = {0};
        if (stat(output_dir.c_str(), &st_output) == -1) {
            mkdir(output_dir.c_str(), 0770);
        }
        if (stat(folder.c_str(), &st_folder) == -1) {
            mkdir(folder.c_str(), 0770);
        }
        FILE *myfile = fopen(filename, "a+");
        if (myfile != NULL) {
            fprintf(myfile, "padcheck %d %d %d %g %d %12.12e %12.12e %12.12e %12.12e\n", nc[0], nc[1], nc[2], LzoLx, pad, dmax, umax,
                    dUjmax / Ujmax, dUorel);
            fclose(myfile);
        }
    }

    flups_free(u1);
    flups_free(up);
    MPI_Barrier(comm);
}

//==============================================================================
// axially modulated ("wavy") vortex tube (--source=2)
//==============================================================================
/**
 * @brief Gauss-Legendre nodes and weights on [-1,1] (Newton iteration on the Legendre polynomial)
 */
static void wavy_gauss_legendre_(const int n, std::vector<double> &x, std::vector<double> &w) {
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z  = cos(M_PI * (i + 0.75) / (n + 0.5));
        double pp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = 0.0;
            for (int j = 0; j < n; ++j) {
                const double p2 = p1;
                p1              = p0;
                p0              = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1.0);
            }
            pp             = n * (z * p0 - p1) / (z * z - 1.0);
            const double dz = -p0 / pp;
            z += dz;
            if (fabs(dz) < 1.0e-16) break;
        }
        x[i]         = -z;
        x[n - 1 - i] = z;
        w[i]         = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

/**
 * @brief the radial part of the wavy-tube solution for one axial harmonic
 *
 * Everything is evaluated in the variable @f$ v = t/(1-t) @f$ with @f$ t = (r/R)^2 @f$, in which
 * the compactly supported bump becomes a plain exponential:
 * @f[ \omega(s) = A e^{-1/(1-t)} = A e^{-1} e^{-v}, \qquad s\,ds = \frac{R^2}{2 (1+v)^2} dv @f]
 * so that any radial moment is a Gauss-Laguerre-like integral of an *analytic* integrand:
 * @f[ \int_0^r F(s)\,\omega(s)\, s\, ds = A e^{-1} R^2 \int_0^{v(r)} \frac{F(s(v))}{2(1+v)^2} e^{-v} dv @f]
 * The C-infinity (but non-analytic) flat point at @f$ s = R @f$ is pushed to @f$ v = \infty @f$
 * where the integrand is exponentially small, and the only remaining singularity is the
 * logarithm of @f$ K_0 @f$ at @f$ v = 0 @f$, handled by geometric panels: on a panel [a, 2a] the
 * logarithm is analytic and Gauss-Legendre converges geometrically. The cumulative integrals at
 * the panel boundaries are tabulated once, so evaluating psi at one radius costs a single
 * partial-panel quadrature.
 */
struct WavyRadial {
    double R    = 0.0;  //!< the tube radius
    double k    = 0.0;  //!< the axial wavenumber
    double amp  = 0.0;  //!< the vorticity amplitude A
    double pref = 0.0;  //!< A exp(-1) R^2

    std::vector<double> vb;      //!< the panel boundaries in v
    std::vector<double> ca, cb;  //!< the cumulative integrals at the panel boundaries
    std::vector<double> gx, gw;  //!< the Gauss-Legendre nodes/weights on [-1,1]
    double              atot = 0.0;  //!< int_0^R I_0(ks) omega(s) s ds
    double              btot = 0.0;  //!< int_0^R K_0(ks) omega(s) s ds

    double s_of_v(const double v) const { return R * sqrt(v / (1.0 + v)); }

    void panel_(const double a, const double b, double &ia, double &ib) const {
        const double c = 0.5 * (a + b);
        const double hh = 0.5 * (b - a);
        ia = 0.0;
        ib = 0.0;
        for (size_t q = 0; q < gx.size(); ++q) {
            const double v = c + hh * gx[q];
            const double s = s_of_v(v);
            const double e = exp(-v) * gw[q] * hh / (2.0 * (1.0 + v) * (1.0 + v));
            ia += e * std::cyl_bessel_i(0.0, k * s);
            ib += e * std::cyl_bessel_k(0.0, k * s);
        }
    }

    void build(const double aR, const double ak, const double aamp, const int nq, const double vmin, const double vmax) {
        R    = aR;
        k    = ak;
        amp  = aamp;
        pref = amp * exp(-1.0) * R * R;
        wavy_gauss_legendre_(nq, gx, gw);
        vb.clear();
        vb.push_back(0.0);
        for (double v = vmin; v < vmax; v *= 2.0) vb.push_back(v);
        vb.push_back(vmax);
        const int np = (int)vb.size() - 1;
        ca.assign(np + 1, 0.0);
        cb.assign(np + 1, 0.0);
        for (int p = 0; p < np; ++p) {
            double ia, ib;
            panel_(vb[p], vb[p + 1], ia, ib);
            ca[p + 1] = ca[p] + ia;
            cb[p + 1] = cb[p] + ib;
        }
        atot = pref * ca[np];
        btot = pref * cb[np];
    }

    //! @brief `a = int_0^r I_0(ks) omega s ds`, `bp = int_0^r K_0(ks) omega s ds`
    void partial_(const double r, double &a, double &bp) const {
        const int    np = (int)vb.size() - 1;
        const double t  = (r / R) * (r / R);
        if (t >= 1.0) {
            a  = atot;
            bp = btot;
            return;
        }
        const double v = t / (1.0 - t);
        if (v >= vb[np]) {
            a  = atot;
            bp = btot;
            return;
        }
        int p = (int)(std::upper_bound(vb.begin(), vb.end(), v) - vb.begin()) - 1;
        p     = std::max(p, 0);
        double ia, ib;
        panel_(vb[p], v, ia, ib);
        a  = pref * (ca[p] + ia);
        bp = pref * (cb[p] + ib);
    }

    //! @brief the solution of `-(lap_2D - k^2) psi = omega`
    double psi(const double r) const {
        double a, bp;
        partial_(r, a, bp);
        const double b = btot - bp;
        if (k * r < 1.0e-13) return b;  // psi(0) = B(0) = btot, the K_0 A product vanishes
        return std::cyl_bessel_k(0.0, k * r) * a + std::cyl_bessel_i(0.0, k * r) * b;
    }
};

//! @brief the compact C-infinity vorticity profile
static inline double wavy_omega_(const double r, const double R, const double amp) {
    const double t = (r / R) * (r / R);
    return (t >= 1.0) ? 0.0 : amp * exp(-1.0 / (1.0 - t));
}

/**
 * @brief prints the self-check of the analytical solution on rank 0
 *
 * Three independent checks:
 *  -# the circulation `2 pi int_0^R omega s ds` must be 1,
 *  -# the radial quadrature must be converged: psi is recomputed with a finer rule,
 *  -# `(lap_2D - k^2) psi + omega` must vanish, by a 4th-order centred finite difference. The
 *     finite difference has its own round-off floor of ~1e-9 in absolute value (~1e-10 relative
 *     to max|omega|), so check 2 is the sharp statement about the accuracy of psi.
 */
static void wavy_selfcheck_(const double R, const double amp, const double *kk, const int nmode) {
    // ..... 1. the circulation .................................................
    // 2 pi int_0^R omega s ds, in the same v variable (the integrand is then just e^-v)
    {
        std::vector<double> gx, gw;
        wavy_gauss_legendre_(32, gx, gw);
        std::vector<double> vb(1, 0.0);
        for (double v = 1.0e-18; v < 40.0; v *= 2.0) vb.push_back(v);
        vb.push_back(40.0);
        double tot = 0.0;
        for (size_t p = 0; p + 1 < vb.size(); ++p) {
            const double c = 0.5 * (vb[p] + vb[p + 1]);
            const double hh = 0.5 * (vb[p + 1] - vb[p]);
            for (size_t q = 0; q < gx.size(); ++q) {
                const double v = c + hh * gx[q];
                tot += exp(-v) * gw[q] * hh / (2.0 * (1.0 + v) * (1.0 + v));
            }
        }
        tot *= amp * exp(-1.0) * R * R;
        m_log("wavytube-selfcheck: circulation = %.16e (err %.3e)", 2.0 * M_PI * tot, 2.0 * M_PI * tot - 1.0);
    }

    for (int j = 0; j < nmode; ++j) {
        const double k = kk[j];
        WavyRadial   ref, fine;
        ref.build(R, k, amp, 24, 1.0e-18, 40.0);
        fine.build(R, k, amp, 48, 1.0e-24, 60.0);

        // ..... 2. the quadrature convergence ..................................
        double dev = 0.0, scale = 0.0;
        for (int i = 0; i < 200; ++i) {
            const double r = (i + 0.5) * (4.0 * R / 200.0);
            const double p = ref.psi(r);
            scale          = std::max(scale, fabs(p));
            dev            = std::max(dev, fabs(p - fine.psi(r)));
        }

        // ..... 3. the finite-difference residual of the radial equation .......
        const double d      = 1.0e-3 * R;
        double       resmax = 0.0;
        for (int i = 1; i <= 40; ++i) {
            const double r = i * (3.0 * R / 40.0);
            if (r - 2.0 * d <= 0.0) continue;
            const double pm2 = ref.psi(r - 2.0 * d), pm1 = ref.psi(r - d), p0 = ref.psi(r);
            const double pp1 = ref.psi(r + d), pp2 = ref.psi(r + 2.0 * d);
            const double d2  = (-pm2 + 16.0 * pm1 - 30.0 * p0 + 16.0 * pp1 - pp2) / (12.0 * d * d);
            const double d1  = (pm2 - 8.0 * pm1 + 8.0 * pp1 - pp2) / (12.0 * d);
            resmax           = std::max(resmax, fabs(d2 + d1 / r - k * k * p0 + wavy_omega_(r, R, amp)));
        }
        m_log("wavytube-selfcheck: mode %d k=%.6e |psi|max=%.6e | quadrature dev=%.3e (rel %.3e) | FD residual of (lap_2D-k^2)psi+omega = %.3e (rel to max|omega| %.3e)",
              j + 1, k, scale, dev, dev / scale, resmax, resmax / (amp * exp(-1.0)));
    }
}

/**
 * @brief solves the axially modulated vortex tube and reports the error against the exact solution
 *
 * @param myCase description of the domain, must be UNB-UNB in two directions and PER (or
 *               EVEN-EVEN / ODD-ODD) in the third one
 * @param typeGreen the Green kernel
 * @param nSolve the number of solves (for timing)
 * @param output_dir the output directory
 */
void wavytube_3d(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int nSolve, const std::string output_dir) {
    int      rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    //-------------------------------------------------------------------------
    /** - identify the two unbounded directions and the spectral one */
    //-------------------------------------------------------------------------
    int dir_u[2] = {-1, -1};
    int dir_s    = -1;
    int n_unb    = 0;
    for (int d = 0; d < 3; ++d) {
        if (myCase.mybc[d][0] == UNB && myCase.mybc[d][1] == UNB) {
            m_assert(n_unb < 2, "the wavy tube needs exactly 2 unbounded directions");
            dir_u[n_unb] = d;
            n_unb += 1;
        } else {
            m_assert(dir_s < 0, "the wavy tube needs exactly 1 non-unbounded direction");
            dir_s = d;
        }
    }
    m_assert(n_unb == 2 && dir_s >= 0, "the wavy tube needs 2 UNB-UNB directions and 1 spectral one");

    // the axial profile, chosen so that c_j'' = -k_j^2 c_j *and* the symmetry of the axial BC is
    // satisfied.  0 = cos(k z - phi) [PER], 1 = cos(k z) [EVEN-EVEN], 2 = sin(k z) [ODD-ODD]
    int axial_kind = -1;
    if (myCase.mybc[dir_s][0] == PER && myCase.mybc[dir_s][1] == PER) {
        axial_kind = 0;
    } else if (myCase.mybc[dir_s][0] == EVEN && myCase.mybc[dir_s][1] == EVEN) {
        axial_kind = 1;
    } else if (myCase.mybc[dir_s][0] == ODD && myCase.mybc[dir_s][1] == ODD) {
        axial_kind = 2;
    }
    m_assert(axial_kind >= 0, "the axial direction of the wavy tube must be PER-PER, EVEN-EVEN or ODD-ODD (got %d %d)",
             myCase.mybc[dir_s][0], myCase.mybc[dir_s][1]);

    //-------------------------------------------------------------------------
    /** - the grid */
    //-------------------------------------------------------------------------
    const bool   is_cell = myCase.center_type[0] == CELL_CENTER;
    const double fact    = (double)(!is_cell);
    const double shift   = is_cell ? 0.5 : 0.0;

    const int    *nproc   = myCase.nproc;
    const double *L       = myCase.L;
    const int     nglob[3] = {myCase.nglob[0] + (int)fact, myCase.nglob[1] + (int)fact, myCase.nglob[2] + (int)fact};
    const double  h[3]     = {L[0] / (nglob[0] - fact), L[1] / (nglob[1] - fact), L[2] / (nglob[2] - fact)};

    FLUPS_CenterType center_type[3];
    for (int d = 0; d < 3; ++d) {
        center_type[d] = myCase.center_type[d];
    }
    FLUPS_BoundaryType *mybc[3][2];
    for (int d = 0; d < 3; ++d) {
        for (int s = 0; s < 2; ++s) {
            mybc[d][s]    = (FLUPS_BoundaryType *)flups_malloc(sizeof(int));
            mybc[d][s][0] = myCase.mybc[d][s];
        }
    }

    FLUPS_Topology *topo = flups_topo_new(0, 1, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);

    std::string     name = "wavytube_nx" + std::to_string(nglob[0]) + "_ny" + std::to_string(nglob[1]) + "_nz" + std::to_string(nglob[2]) + "_nrank" + std::to_string(comm_size);
    FLUPS_Profiler *prof = flups_profiler_new_n(name.c_str());

    m_profStart(prof, "Validation--Init-Flups");
    FLUPS_Solver *mysolver = flups_init_timed(topo, mybc, h, L, NOD, center_type, prof);
    flups_set_greenType(mysolver, typeGreen);
    flups_setup(mysolver, true);
    comm = flups_topo_get_comm(topo);
    MPI_Comm_rank(comm, &rank);
    m_profStop(prof, "Validation--Init-Flups");

    //-------------------------------------------------------------------------
    /** - the analytical source and solution */
    //-------------------------------------------------------------------------
    m_profStart(prof, "Validation--Init-RHS");
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *sol   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(sol, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));

    // the tube
    const int    nmode      = 2;
    const double aj[nmode]  = {1.0, 0.5};
    const double phi[nmode] = {0.0, 0.5 * M_PI};
    const double Lz         = L[dir_s];
    const double Rtube      = 0.2 * std::min(L[dir_u[0]], L[dir_u[1]]);
    // A = 1 / (pi R^2 E_2(1)) so that int omega dA = 1, with E_2(1) = e^-1 - E_1(1) = e^-1 + Ei(-1)
    const double E2one = exp(-1.0) + std::expint(-1.0);
    const double amp   = 1.0 / (M_PI * Rtube * Rtube * E2one);
    double       kk[nmode];
    for (int j = 0; j < nmode; ++j) {
        kk[j] = c_2pi * (j + 1) / Lz;
    }
    const double ctr[2] = {0.5 * L[dir_u[0]], 0.5 * L[dir_u[1]]};

    if (rank == 0) {
        m_log("wavytube: R=%.6e A=%.6e Lz=%.6e k=(%.6e,%.6e) axial_kind=%d (0=PER,1=EVEN,2=ODD)", Rtube, amp, Lz, kk[0], kk[1], axial_kind);
        wavy_selfcheck_(Rtube, amp, kk, nmode);
    }

    WavyRadial rad[nmode];
    for (int j = 0; j < nmode; ++j) {
        rad[j].build(Rtube, kk[j], amp, 24, 1.0e-18, 40.0);
    }

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);
    const int ax0     = flups_topo_get_axis(topo);
    const int ax1     = (ax0 + 1) % 3;
    const int ax2     = (ax0 + 2) % 3;
    const int nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};

    // tabulate omega and psi_j on the local slab of the transverse plane -- the radial functions
    // are the expensive part and they only depend on the two unbounded coordinates
    const int    ntr0 = flups_topo_get_nloc(topo, dir_u[0]);
    const int    ntr1 = flups_topo_get_nloc(topo, dir_u[1]);
    const size_t ntr  = ((size_t)ntr0) * ntr1;
    std::vector<double> tab_om(ntr, 0.0);
    std::vector<double> tab_psi(nmode * ntr, 0.0);
    for (int l1 = 0; l1 < ntr1; ++l1) {
        for (int l0 = 0; l0 < ntr0; ++l0) {
            const double x0 = (istart[dir_u[0]] + l0 + shift) * h[dir_u[0]] - ctr[0];
            const double x1 = (istart[dir_u[1]] + l1 + shift) * h[dir_u[1]] - ctr[1];
            const double r  = sqrt(x0 * x0 + x1 * x1);
            const size_t il = l0 + ((size_t)ntr0) * l1;
            tab_om[il]      = wavy_omega_(r, Rtube, amp);
            for (int j = 0; j < nmode; ++j) {
                tab_psi[j * ntr + il] = rad[j].psi(r);
            }
        }
    }

    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                int          g[3];
                g[ax0] = istart[ax0] + i0;
                g[ax1] = istart[ax1] + i1;
                g[ax2] = istart[ax2] + i2;

                const size_t il = (g[dir_u[0]] - istart[dir_u[0]]) + ((size_t)ntr0) * (g[dir_u[1]] - istart[dir_u[1]]);
                const double z  = (g[dir_s] + shift) * h[dir_s];

                double f = 0.0, u = 0.0;
                for (int j = 0; j < nmode; ++j) {
                    const double arg = kk[j] * z;
                    const double cj  = (axial_kind == 2) ? sin(arg) : cos(arg - ((axial_kind == 0) ? phi[j] : 0.0));
                    f -= aj[j] * cj * tab_om[il];
                    u += aj[j] * cj * tab_psi[j * ntr + il];
                }
                rhs[id] = f;
                sol[id] = u;
            }
        }
    }
    m_profStop(prof, "Validation--Init-RHS");

    //-------------------------------------------------------------------------
    /** - solve */
    //-------------------------------------------------------------------------
    for (int is = 0; is < nSolve; is++) {
        MPI_Barrier(comm);
        m_profStart(prof, "Validation--Solve");
        flups_solve(mysolver, field, rhs, STD);
        m_profStop(prof, "Validation--Solve");
        MPI_Barrier(comm);
    }
    flups_profiler_disp(prof);

    //-------------------------------------------------------------------------
    /** - the error */
    //-------------------------------------------------------------------------
    const double vol = h[0] * h[1] * h[2];
    double       lerr2 = 0.0, lerri = 0.0, err2 = 0.0, erri = 0.0;
    for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
        for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
            for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                const size_t id  = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                const double err = sol[id] - field[id];
                lerri            = std::max(lerri, fabs(err));
                lerr2 += err * err * vol;
            }
        }
    }
    MPI_Allreduce(&lerr2, &err2, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&lerri, &erri, 1, MPI_DOUBLE, MPI_MAX, comm);
    err2 = sqrt(err2);

    const int ntrans = nglob[dir_u[0]];
    if (rank == 0) {
        printf("wavytube: N=%d,%d,%d Ntrans=%d Lz/Lx=%g kernel=%d err2=%12.12e errinf=%12.12e\n",
               nglob[0], nglob[1], nglob[2], ntrans, L[dir_s] / L[dir_u[0]], (int)typeGreen, err2, erri);

        string folder  = output_dir + "/data";
        string ct_name = is_cell ? "CellCenter" : "NodeCenter";
        char   filename[PATH_MAX];
        // NOTE: deliberately NOT the file name used by validation_3d, so that the data_ref
        // regression of scripts/check_res.py is never confused by these lines
        sprintf(filename, "%s/wavytube_%s_%d%d%d%d%d%d_typeGreen=%d.txt", folder.c_str(), ct_name.c_str(),
                myCase.mybc[0][0], myCase.mybc[0][1], myCase.mybc[1][0], myCase.mybc[1][1], myCase.mybc[2][0], myCase.mybc[2][1], typeGreen);

        struct stat st_output = {0};
        struct stat st_folder = {0};
        if (stat(output_dir.c_str(), &st_output) == -1) {
            mkdir(output_dir.c_str(), 0770);
        }
        if (stat(folder.c_str(), &st_folder) == -1) {
            mkdir(folder.c_str(), 0770);
        }
        FILE *myfile = fopen(filename, "a+");
        if (myfile != NULL) {
            printf("Opening file %s !\n", filename);
            fprintf(myfile, "%d %12.12e %12.12e\n", ntrans, err2, erri);
            fclose(myfile);
        } else {
            printf("unable to open file %s ! Here is what I would have written: %d %12.12e %12.12e\n", filename, ntrans, err2, erri);
        }
    }
    MPI_Barrier(comm);

    flups_free(sol);
    flups_free(rhs);
    flups_free(field);
    flups_cleanup_solver(mysolver);
    flups_cleanup_backend();
    flups_profiler_free(prof);
    flups_topo_free(topo);
    for (int d = 0; d < 3; ++d) {
        for (int s = 0; s < 2; ++s) {
            flups_free(mybc[d][s]);
        }
    }
}
