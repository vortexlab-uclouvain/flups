/**
 * @file validation.cpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-19
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "validation_3d.hpp"

/**
 * @brief computes the reference solution and the numerical one, outputs errors in a file
 * 
 * @param myCase description of the domain and initial condition
 * @param type type of solver
 * @param typeGreen type of Green function
 */
void validation_3d(const DomainDescr myCase, const SolverType type, const GreenType typeGreen) {
    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int *   nglob  = myCase.nglob;
    const int *   nproc  = myCase.nproc;
    const double *L      = myCase.L;
    const double  sigma  = myCase.sigma;
    const double *center = myCase.center;

    const double h[3] = {L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]};

    const BoundaryType mybc[3][2] = {myCase.mybc[0][0], myCase.mybc[0][1],
                                     myCase.mybc[1][0], myCase.mybc[1][1],
                                     myCase.mybc[2][0], myCase.mybc[2][1]};

    // create a real topology
    const Topology *topo = new Topology(0, nglob, nproc, false);

    //-------------------------------------------------------------------------
    /** - Initialize the solver */
    //-------------------------------------------------------------------------
    FFTW_Solver *mysolver = new FFTW_Solver(topo, mybc, h, L);
    mysolver->set_GreenType(typeGreen);
    mysolver->setup();

    //-------------------------------------------------------------------------
    /** - allocate rhs and solution */
    //-------------------------------------------------------------------------
    double *rhs = (double *)fftw_malloc(sizeof(double *) * topo->locmemsize());
    double *sol = (double *)fftw_malloc(sizeof(double *) * topo->locmemsize());

    //-------------------------------------------------------------------------
    /** - fill the rhs and the solution */
    //-------------------------------------------------------------------------
    const double oosigma  = 1.0 / (sigma);
    const double oosigma2 = 1.0 / (sigma * sigma);
    const double oosigma3 = 1.0 / (sigma * sigma * sigma);

    int istart[3];
    get_idstart_glob(istart, topo);

    /**
         * @todo change that to axis-based loops
         */
    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                double       x    = (istart[0] + i0 + 0.5) * h[0] - L[0] * center[0];
                double       y    = (istart[1] + i1 + 0.5) * h[1] - L[1] * center[1];
                double       z    = (istart[2] + i2 + 0.5) * h[2] - L[2] * center[2];
                double       rho2 = (x * x + y * y + z * z) * oosigma2;
                double       rho  = sqrt(rho2);
                const size_t id   = localindex_xyz(i0, i1, i2, topo);
                // Gaussian
                rhs[id] = -c_1o4pi * oosigma3 * sqrt(2.0 / M_PI) * exp(-rho2 * 0.5);
                sol[id] = +c_1o4pi * oosigma * 1.0 / rho * erf(rho * c_1osqrt2);
            }
        }
    }

    // add contribution from symmetry conditions
    for (int dir = 0; dir < 3; dir++) {
        for (int lr = 0; lr < 2; lr++) {
            double centerPos[3] = {L[0] * center[0], L[1] * center[1], L[2] * center[2]};
            double sign         = 1.0;
            if (mybc[dir][lr] == PER) {
                centerPos[dir] = L[dir] * center[dir] + (2 * lr - 1) * L[dir];
            } else if (mybc[dir][lr] == EVEN) {
                centerPos[dir] = -L[dir] * center[dir] + 2 * L[dir] * lr;
            } else if (mybc[dir][lr] == ODD) {
                centerPos[dir] = -L[dir] * center[dir] + 2 * L[dir] * lr;
                sign           = -1.0;
            } else {
                continue;
            }

            for (int i2 = 0; i2 < topo->nloc(2); i2++) {
                for (int i1 = 0; i1 < topo->nloc(1); i1++) {
                    for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                        double       x    = (istart[0] + i0 + 0.5) * h[0] - centerPos[0];
                        double       y    = (istart[1] + i1 + 0.5) * h[1] - centerPos[1];
                        double       z    = (istart[2] + i2 + 0.5) * h[2] - centerPos[2];
                        double       rho2 = (x * x + y * y + z * z) * oosigma2;
                        double       rho  = sqrt(rho2);
                        const size_t id   = localindex_xyz(i0, i1, i2, topo);

                        // Gaussian
                        rhs[id] -= sign * c_1o4pi * oosigma3 * sqrt(2.0 / M_PI) * exp(-rho2 * 0.5);
                        sol[id] += sign * c_1o4pi * oosigma * 1.0 / rho * erf(rho * c_1osqrt2);
                    }
                }
            }
        }
    }

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
    mysolver->solve(topo, rhs, rhs, UP_SRHS);

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

    /**
         * @todo change that to axis-based loops
         */
    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                const size_t id  = localindex_xyz(i0, i1, i2, topo);
                const double err = sol[id] - rhs[id];

                lerri = max(lerri, abs(err));
                lerr2 += (err * err) * h[0] * h[1] * h[2];
            }
        }
    }
    double erri = 0.0;
    double err2 = 0.0;
    MPI_Allreduce(&lerr2, &err2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&lerri, &erri, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    err2 = sqrt(err2);

    char filename[512];
    sprintf(filename, "data/%s_%d%d%d%d%d%d_typeGreen=%d.err",__func__, mybc[0][0], mybc[0][1], mybc[1][0], mybc[1][1], mybc[2][0], mybc[2][1],typeGreen);

    if (rank == 0) {
        FILE *myfile = fopen(filename, "a+");
        if (myfile != NULL) {
            fprintf(myfile, "%d %12.12e %12.12e\n", nglob[0], err2, erri);
            fclose(myfile);
        } else {
            UP_CHECK1(false, "unable to open file %s", filename);
        }
    }

    fftw_free(sol);
    fftw_free(rhs);
    delete (mysolver);
    delete (topo);
}