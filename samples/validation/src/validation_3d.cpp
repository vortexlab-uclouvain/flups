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
    std::string     name = "validation_nx" + std::to_string((int)(nglob[0])) + "_ny" + std::to_string((int)(nglob[1])) + "_nz" + std::to_string((int)(nglob[2])) + "_nrank" + std::to_string(comm_size) + "_nthread" + std::to_string(1);
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
    flups_cleanup(mysolver);
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
//==============================================================================

/**
 * @brief the compact Spietz bump, support |t| < 1
 */
static inline double padcheck_bump_(const double t) {
    return (fabs(t) >= 1.0) ? 0.0 : exp(c_C * (1.0 - 1.0 / (1.0 - t * t)));
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
    const double shift   = is_cell ? 0.5 : 0.0;

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

    // the source: a bump centred in the *reference* box, of radius 0.35 of the reference box
    double bump_c[3], bump_r[3];
    for (int d = 0; d < 3; ++d) {
        bump_c[d] = 0.50 * myCase.L[d];
        bump_r[d] = 0.35 * myCase.L[d];
    }
    const double Lz = myCase.L[dir_s];

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
                const double xu0  = (g[dir_u[0]] + shift) * h[dir_u[0]];
                const double xu1  = (g[dir_u[1]] + shift) * h[dir_u[1]];
                const double z    = (g[dir_s] + shift) * h[dir_s];
                rhs[id]           = padcheck_bump_((xu0 - bump_c[dir_u[0]]) / bump_r[dir_u[0]]) *
                          padcheck_bump_((xu1 - bump_c[dir_u[1]]) / bump_r[dir_u[1]]) *
                          (1.0 + sin(c_2pi * z / Lz) + 0.5 * cos(2.0 * c_2pi * z / Lz));
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
    flups_cleanup(mysolver);
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
    /** - the difference, and the difference of the spectral-mean-free part */
    //-------------------------------------------------------------------------
    // the mean along the spectral direction, at each point of the two unbounded ones
    const int nsz  = nc[dir_s];
    double    dmax = 0.0, dmax_nomean = 0.0, umax = 0.0;
    {
        // loop over the (dir_u[0], dir_u[1]) plane, the spectral direction inside
        const int nu0 = nc[dir_u[0]];
        const int nu1 = nc[dir_u[1]];
        for (int j1 = 0; j1 < nu1; ++j1) {
            for (int j0 = 0; j0 < nu0; ++j0) {
                double m1 = 0.0, mp = 0.0;
                for (int js = 0; js < nsz; ++js) {
                    int g[3];
                    g[dir_u[0]]      = j0;
                    g[dir_u[1]]      = j1;
                    g[dir_s]         = js;
                    const size_t idx = g[0] + nc[0] * (g[1] + ((size_t)nc[1]) * g[2]);
                    m1 += u1[idx];
                    mp += up[idx];
                }
                m1 /= nsz;
                mp /= nsz;
                for (int js = 0; js < nsz; ++js) {
                    int g[3];
                    g[dir_u[0]]      = j0;
                    g[dir_u[1]]      = j1;
                    g[dir_s]         = js;
                    const size_t idx = g[0] + nc[0] * (g[1] + ((size_t)nc[1]) * g[2]);
                    umax             = max(umax, fabs(u1[idx]));
                    dmax             = max(dmax, fabs(u1[idx] - up[idx]));
                    dmax_nomean      = max(dmax_nomean, fabs((u1[idx] - m1) - (up[idx] - mp)));
                }
            }
        }
    }

    const double LzoLx = myCase.L[dir_s] / myCase.L[dir_u[0]];
    if (rank == 0) {
        printf("padcheck: N=%d,%d,%d Lz/Lx=%g pad=%d kernel=%d max|dU|=%12.12e max|dU_nomean|=%12.12e max|U|=%12.12e\n",
               nc[0], nc[1], nc[2], LzoLx, pad, (int)typeGreen, dmax, dmax_nomean, umax);

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
            fprintf(myfile, "padcheck %d %d %d %g %d %12.12e %12.12e %12.12e\n", nc[0], nc[1], nc[2], LzoLx, pad, dmax, dmax_nomean, umax);
            fclose(myfile);
        }
    }

    flups_free(u1);
    flups_free(up);
    MPI_Barrier(comm);
}
