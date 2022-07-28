/**
 * @file vtube.cpp
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

#include "vtube.hpp"
#include "omp.h"


using namespace std;

void vtube(const DomainDescr myCase, const FLUPS_GreenType typeGreen, const int nSolve, int type, FLUPS_DiffType order, int vdir) {
    int rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    const int lda =3;

    const bool tube = (type == 0);
    const bool ring = (type == 1);

    const int    *nproc = myCase.nproc;
    const double *L     = myCase.L;

    const bool is_cell = myCase.center[0] == CELL_CENTER; 
    const double fact = (double) (!is_cell);

    // Ensure that 64 points in cell-centred mode is equivalent 64 points in node-centred mode
    const int    nglob[3] = {myCase.nglob[0] + (int)fact,
                             myCase.nglob[1] + (int)fact,
                             myCase.nglob[2] + (int)fact};

    // h definition depends on if we are cell- of node-centred
    const double h[3] = {L[0] / (nglob[0] - fact),
                         L[1] / (nglob[1] - fact),
                         L[2] / (nglob[2] - fact)};

    FLUPS_CenterType center_type[3];
    for (int i = 0; i < 3; i++) {
            center_type[i] = myCase.center[i];
    }

    FLUPS_BoundaryType* mybc[3][2];
    for(int id=0; id<3; id++){
        for(int is=0; is<2; is++){
            mybc[id][is] = (FLUPS_BoundaryType*) flups_malloc(sizeof(int)*lda);
        }
    }

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            for (int lia = 0; lia < 3; lia++) {
                mybc[id][is][lia] = myCase.mybcv[id][is][lia];
            }
        }
    }

    // create a real topology
    FLUPS_Topology *topo = flups_topo_new(0, 3, nglob, nproc, false, NULL, FLUPS_ALIGNMENT, comm);

    //-------------------------------------------------------------------------
    /** - Initialize the solver */
    //-------------------------------------------------------------------------
    std::string     name = "tube_" + std::to_string((int)(nglob[0] / L[0])) + "_nrank" + std::to_string(comm_size) + "_nthread" + std::to_string(omp_get_max_threads());
    FLUPS_Profiler *prof = flups_profiler_new_n(name.c_str());
    FLUPS_Solver *  mysolver;
    mysolver = flups_init_timed(topo, mybc, h, L, order, center_type, prof);

    flups_set_greenType(mysolver, typeGreen);
    flups_setup(mysolver, true);

    // update the comm and the rank
    comm = flups_topo_get_comm(topo);
    MPI_Comm_rank(comm, &rank);

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    double *rhs   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *sol   = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *field = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    double *error = (double *)flups_malloc(sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(rhs, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(sol, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(field, 0, sizeof(double) * flups_topo_get_memsize(topo));
    std::memset(error, 0, sizeof(double) * flups_topo_get_memsize(topo));

    //-------------------------------------------------------------------------
    /** - fill the rhs and the solution */
    //-------------------------------------------------------------------------

    int istart[3];
    flups_topo_get_istartGlob(topo, istart);

    {
        const double sigma   = 0.1;
        const double rad     = 0.3;
        const double oo_rad2 = 1.0/(rad*rad);
        const int    ax0     = flups_topo_get_axis(topo);
        const int    ax1     = (ax0 + 1) % 3;
        const int    ax2     = (ax0 + 2) % 3;
        const int    nmem[3] = {flups_topo_get_nmem(topo, 0), flups_topo_get_nmem(topo, 1), flups_topo_get_nmem(topo, 2)};

        double *rhs0 = rhs + flups_locID(ax0, 0, 0, 0, 0, ax0, nmem, 1);
        double *rhs1 = rhs + flups_locID(ax0, 0, 0, 0, 1, ax0, nmem, 1);
        double *rhs2 = rhs + flups_locID(ax0, 0, 0, 0, 2, ax0, nmem, 1);
        double *sol0 = sol + flups_locID(ax0, 0, 0, 0, 0, ax0, nmem, 1);
        double* sol1 = sol + flups_locID(ax0, 0, 0, 0, 1, ax0, nmem, 1);
        double* sol2 = sol + flups_locID(ax0, 0, 0, 0, 2, ax0, nmem, 1);

        int dir0 = (vdir + 1) % 3;
        int dir1 = (vdir + 2) % 3;
        int dir2 = vdir;

        // printf("center = %f %f - signs = %f %f",myCase.xcntr,myCase.ycntr,myCase.xsign,myCase.ysign);

        for (int i2 = 0; i2 < flups_topo_get_nloc(topo, ax2); i2++) {
            for (int i1 = 0; i1 < flups_topo_get_nloc(topo, ax1); i1++) {
                for (int i0 = 0; i0 < flups_topo_get_nloc(topo, ax0); i0++) {
                    const size_t id     = flups_locID(ax0, i0, i1, i2, 0, ax0, nmem, 1);
                    const double shift  = is_cell ? 0.5 : 0.0;
                    const double pos[3] = {(istart[ax0] + i0 + shift) * h[ax0],
                                           (istart[ax1] + i1 + shift) * h[ax1],
                                           (istart[ax2] + i2 + shift) * h[ax2]};
                    if (tube) {
                        //---------------------------------------------------------------
                        // main tube
                        {
                            
                            const double x     = pos[dir0] - (myCase.xcntr * L[dir0]);
                            const double y     = pos[dir1] - (myCase.ycntr * L[dir1]);
                            const double theta = std::atan2(y, x);  // get the angle in the x-y plane
                            const double r     = sqrt(x * x + y * y);
                            const double rho1  = r / sigma;
                            const double rho2  = r / rad;
                            
                            const double E21 = gexpint<2>(1.0);
                            const double vort  = (rho2 >= 1.0) ? (0.0) : (1.0 / c_2pi * (2.0 * oo_rad2) * (1.0 / E21) * exp(-1.0 / (1.0 - rho2 * rho2)));
                            const double tol = 100*std::numeric_limits<double>::epsilon();
                            const double fact =  (rho2 >= 1.0) ? (1.0) : (1.0 - (1.0 / E21) * (1.0 - rho2 * rho2) * gexpint<2>(1.0 / (1.0 - rho2 * rho2)));
                            const double vel   =(r < tol) ? (0.0) : (fact/ (c_2pi * r));

                            // const double vel   = (r < 100*std::numeric_limits<double>::epsilon()) ? 0 : 1.0 / (c_2pi * r) * (1.0 - exp(-rho1 * rho1 * 0.5));
                            // const double vort  = 1.0 / (c_2pi * sigma * sigma) * exp(-rho1 * rho1 * 0.5);

                            if (dir2 == 0) {
                                rhs0[id] = -vort;
                                rhs1[id] = 0.0;
                                rhs2[id] = 0.0;
                                sol0[id] = 0.0;
                                sol1[id] = -sin(theta) * vel;
                                sol2[id] = +cos(theta) * vel;
                            } else if (dir2 == 1) {
                                rhs0[id] = 0.0;
                                rhs1[id] = -vort;
                                rhs2[id] = 0.0;
                                sol0[id] = +cos(theta) * vel;
                                sol1[id] = 0.0;
                                sol2[id] = -sin(theta) * vel;
                            } else if (dir2 == 2) {
                                rhs0[id] = 0.0;
                                rhs1[id] = 0.0;
                                rhs2[id] = -vort;
                                sol0[id] = -sin(theta) * vel;
                                sol1[id] = +cos(theta) * vel;
                                sol2[id] = 0.0;
                            }
                        }

                        //---------------------------------------------------------------
                        // x sym tube
                        {
                            const double x     = pos[dir0] + (myCase.xcntr * L[dir0]);
                            const double y     = pos[dir1] - (myCase.ycntr * L[dir1]);
                            const double theta = std::atan2(y, x);
                            const double r     = sqrt(x * x + y * y);
                            const double rho   = r / sigma;
                            const double vel   = 1.0 / (c_2pi * r) * (1.0 - exp(-rho * rho * 0.5));
                            const double vort  = 1.0 / (c_2pi * sigma * sigma) * exp(-rho * rho * 0.5);

                            if (dir2 == 0) {
                                rhs0[id] += -vort * myCase.xsign;
                                rhs1[id] += 0.0 * myCase.xsign;
                                rhs2[id] += 0.0 * myCase.xsign;
                                sol0[id] += 0.0 * myCase.xsign;
                                sol1[id] += -sin(theta) * vel * myCase.xsign;
                                sol2[id] += +cos(theta) * vel * myCase.xsign;
                            } else if (dir2 == 1) {
                                rhs0[id] += 0.0 * myCase.xsign;
                                rhs1[id] += -vort * myCase.xsign;
                                rhs2[id] += 0.0 * myCase.xsign;
                                sol0[id] += +cos(theta) * vel * myCase.xsign;
                                sol1[id] += 0.0 * myCase.xsign;
                                sol2[id] += -sin(theta) * vel * myCase.xsign;
                            } else if (dir2 == 2) {
                                rhs0[id] += 0.0 * myCase.xsign;
                                rhs1[id] += 0.0 * myCase.xsign;
                                rhs2[id] += -vort * myCase.xsign;
                                sol0[id] += -sin(theta) * vel * myCase.xsign;
                                sol1[id] += +cos(theta) * vel * myCase.xsign;
                                sol2[id] += 0.0 * myCase.xsign;
                            }
                        }

                        //---------------------------------------------------------------
                        // y sym tube
                        {
                            const double x     = pos[dir0] - (myCase.xcntr * L[dir0]);
                            const double y     = pos[dir1] + (myCase.ycntr * L[dir1]);
                            const double theta = std::atan2(y, x);
                            const double r     = sqrt(x * x + y * y);
                            const double rho   = r / sigma;
                            const double vel   = 1.0 / (c_2pi * r) * (1.0 - exp(-rho * rho * 0.5));
                            const double vort  = 1.0 / (c_2pi * sigma * sigma) * exp(-rho * rho * 0.5);

                            if (dir2 == 0) {
                                rhs0[id] += -vort * myCase.ysign;
                                rhs1[id] += 0.0 * myCase.ysign;
                                rhs2[id] += 0.0 * myCase.ysign;
                                sol0[id] += 0.0 * myCase.ysign;
                                sol1[id] += -sin(theta) * vel * myCase.ysign;
                                sol2[id] += +cos(theta) * vel * myCase.ysign;
                            } else if (dir2 == 1) {
                                rhs0[id] += 0.0 * myCase.ysign;
                                rhs1[id] += -vort * myCase.ysign;
                                rhs2[id] += 0.0 * myCase.ysign;
                                sol0[id] += +cos(theta) * vel * myCase.ysign;
                                sol1[id] += 0.0 * myCase.ysign;
                                sol2[id] += -sin(theta) * vel * myCase.ysign;
                            } else if (dir2 == 2) {
                                rhs0[id] += 0.0 * myCase.ysign;
                                rhs1[id] += 0.0 * myCase.ysign;
                                rhs2[id] += -vort * myCase.ysign;
                                sol0[id] += -sin(theta) * vel * myCase.ysign;
                                sol1[id] += +cos(theta) * vel * myCase.ysign;
                                sol2[id] += 0.0 * myCase.ysign;
                            }
                        }
                        //---------------------------------------------------------------
                        // xy sym tube
                        {
                            const double x     = pos[dir0] + (myCase.xcntr * L[dir0]);
                            const double y     = pos[dir1] + (myCase.ycntr * L[dir1]);
                            const double theta = std::atan2(y, x);
                            const double r     = sqrt(x * x + y * y);
                            const double rho   = r / sigma;
                            const double vel   = 1.0 / (c_2pi * r) * (1.0 - exp(-rho * rho * 0.5));
                            const double vort  = 1.0 / (c_2pi * sigma * sigma) * exp(-rho * rho * 0.5);

                            if (dir2 == 0) {
                                rhs0[id] += -vort * myCase.ysign * myCase.xsign;
                                rhs1[id] += 0.0 * myCase.ysign * myCase.xsign;
                                rhs2[id] += 0.0 * myCase.ysign * myCase.xsign;
                                sol0[id] += 0.0 * myCase.ysign * myCase.xsign;
                                sol1[id] += -sin(theta) * vel * myCase.ysign * myCase.xsign;
                                sol2[id] += +cos(theta) * vel * myCase.ysign * myCase.xsign;
                            } else if (dir2 == 1) {
                                rhs0[id] += 0.0 * myCase.ysign * myCase.xsign;
                                rhs1[id] += -vort * myCase.ysign * myCase.xsign;
                                rhs2[id] += 0.0 * myCase.ysign * myCase.xsign;
                                sol0[id] += +cos(theta) * vel * myCase.ysign * myCase.xsign;
                                sol1[id] += 0.0 * myCase.ysign * myCase.xsign;
                                sol2[id] += -sin(theta) * vel * myCase.ysign * myCase.xsign;
                            } else if (dir2 == 2) {
                                rhs0[id] += 0.0 * myCase.ysign * myCase.xsign;
                                rhs1[id] += 0.0 * myCase.ysign * myCase.xsign;
                                rhs2[id] += -vort * myCase.ysign * myCase.xsign;
                                sol0[id] += -sin(theta) * vel * myCase.ysign * myCase.xsign;
                                sol1[id] += +cos(theta) * vel * myCase.ysign * myCase.xsign;
                                sol2[id] += 0.0 * myCase.ysign * myCase.xsign;
                            }
                        }
                        }else if (ring){

                            //---------------------------------------------------------------
                            // main ring
                            {
                                const double x     = pos[0] - (myCase.xcntr * L[0]);
                                const double y     = pos[1] - (myCase.ycntr * L[1]);
                                const double z     = pos[2] - (myCase.zcntr * L[2]);

                                // get the coordinates in the radial plan
                                const double xp = sqrt(x * x + y * y);
                                const double alphap = atan2(y,x);

                                // in the plan, the coordinates are (xp,z)
                                // get the distance to the 1st vortex
                                const double r1 = sqrt((xp-rad)*(xp-rad) + z*z);
                                const double theta1 = std::atan2(z,(xp-rad));  // get the angle in the x-y plane
                                const double rho1   = r1 / sigma;
                                const double r2 = sqrt((xp+rad)*(xp+rad) + z*z);
                                const double theta2 = std::atan2(z,(xp+rad));  // get the angle in the x-y plane
                                const double rho2   = r2 / sigma;

                                const double vel1  = 1.0 / (c_2pi * r1) * (1.0 - exp(-rho1 * rho1 * 0.5));
                                const double vel2 = - 1.0 / (c_2pi * r2) * (1.0 - exp(-rho2 * rho2 * 0.5));
                                const double vn = cos(theta1)*vel1+cos(theta2)*vel2;
                                const double vr = -sin(theta1)*vel1-sin(theta2)*vel2;
                                const double vort = +1.0 / (c_2pi * sigma * sigma) * exp(-rho1 * rho1 * 0.5) - 1.0 / (c_2pi * sigma * sigma) * exp(-rho2 * rho2 * 0.5);

                                rhs0[id] = +sin(alphap) * (-vort);
                                rhs1[id] = -cos(alphap) * (-vort);
                                rhs2[id] = 0.0;
                                sol0[id] = cos(alphap) * vr;
                                sol1[id] = sin(alphap) * vr;
                                sol2[id] = vn;
                            }
                        }
                    }
                }
            }
    }
#ifdef DUMP_DBG
    char msg[512];
    // write the source term and the solution
    sprintf(msg, "rhs_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0][0], mybc[0][1][0], mybc[1][0][0], mybc[1][1][0], mybc[2][0][0], mybc[2][1][0], nglob[0], nglob[1], nglob[2]);
    flups_hdf5_dump(topo, msg, rhs);
    sprintf(msg, "anal");
    flups_hdf5_dump(topo, msg, sol);
#endif

    //-------------------------------------------------------------------------
    /** - solve the equations */
    //-------------------------------------------------------------------------
    for (int is = 0; is < nSolve; is++) {
        flups_solve(mysolver, field, rhs, ROT);
    }


    flups_profiler_disp(prof);
    flups_profiler_free(prof);


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

    //determine the volume associated to a mesh
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
                        error[id]  = err; 
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

#ifdef DUMP_DBG
    // write the source term and the solution
    sprintf(msg, "error_%d%d%d%d%d%d_%dx%dx%d", mybc[0][0][0], mybc[0][1][0], mybc[1][0][0], mybc[1][1][0], mybc[2][0][0], mybc[2][1][0], nglob[0], nglob[1], nglob[2]);
    flups_hdf5_dump(topo, msg, error);
#endif

    char   filename[512];
    string folder = "./data";
    string ct_name = is_cell ? "CellCenter" : "NodeCenter"; 
    sprintf(filename, "%s/%s_%s_%d%d%d_typeGreen=%d.txt", folder.c_str(), __func__, ct_name.c_str(), vdir, (int) myCase.xsign, (int)myCase.ysign, typeGreen);

    if (rank == 0) {
        struct stat st = {0};
        if (stat(folder.c_str(), &st) == -1) {
            mkdir(folder.c_str(), 0770);
        }

        FILE *myfile = fopen(filename, "a+");
        if (myfile != NULL) {
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

    free(lerr2);
    free(lerri);
    free(err2);
    free(erri);
    
    flups_free(sol);
    flups_free(rhs);
    flups_free(field);
    flups_free(error);
    flups_cleanup(mysolver);
    flups_topo_free(topo);

    for (int id = 0; id < 3; id++) {
        for (int is = 0; is < 2; is++) {
            flups_free(mybc[id][is]);
        }
    }
}
