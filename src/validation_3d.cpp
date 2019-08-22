/**
 * @file validation.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-19
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "validation_3d.hpp"

/**
 * @brief 
 * 
 * @param type 
 * @param orderdiff 
 */
void validation_3d_UU_UU(const int nsample, const int *size, const SolverType type, const OrderDiff orderdiff)
{

    // int size[2] = {1024,2048};
    // int size[4] = {64,128,256,512};
    // int size[4] = {64,128,256,512};
    // int size[1] = {64};

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char filename[512];
    sprintf(filename, "./data/%s-type=%d-orderdiff=%d.error", __func__, type, orderdiff);

    if (rank == 0)
    {

        FILE *myfile = fopen(filename, "w+");
        if (myfile != NULL)
        {
            fclose(myfile);
        }
        else
        {
            UP_CHECK1(false, "unable to open file %s", filename);
        }
    }

    for (int is = 0; is < nsample; is++)
    {

        int rank, comm_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

        hdf5_dumptest();
        
        const int nglob[3]  = {size[is],size[is],size[is]};
        const int nproc[3]  = {2,2,1};
        const double L[3]   = {1.0,1.0,1.0};
        const double h[3]   = {L[0]/nglob[0],L[1]/nglob[1],L[2]/nglob[2]};
        
        // create a real topology
        const Topology* topo = new Topology(0,nglob,nproc,false);

        const BoundaryType mybc[3][2] = {{UNB, UNB}, {UNB, UNB}, {UNB, UNB}};

        //-------------------------------------------------------------------------
        /** - Initialize the solver */
        //-------------------------------------------------------------------------
        FFTW_Solver *mysolver = new FFTW_Solver(topo, mybc,h,L);
        mysolver->set_GreenType(orderdiff);
        mysolver->setup();

        //-------------------------------------------------------------------------
        /** - allocate rhs and solution */
        //-------------------------------------------------------------------------
        double *rhs = (double *)fftw_malloc(sizeof(double *) * topo->locmemsize());
        double *sol = (double *)fftw_malloc(sizeof(double *) * topo->locmemsize());

        //-------------------------------------------------------------------------
        /** - fill the rhs and the solution */
        //-------------------------------------------------------------------------
        const double sigma = 0.05;
        const double oosigma = 1.0 / (sigma);
        const double oosigma2 = 1.0 / (sigma * sigma);
        const double oosigma3 = 1.0 / (sigma * sigma * sigma);
        const double center[3] = {0.5, 0.5, 0.5};

        int istart[3];
        get_idstart_glob(istart,topo);

        /**
         * @todo change that to axis-based loops
         */
        for (int i2 = 0; i2 < topo->nloc(2); i2++)
        {
            for (int i1 = 0; i1 < topo->nloc(1); i1++)
            {
                for (int i0 = 0; i0 < topo->nloc(0); i0++)
                {
                    double x = (istart[0] + i0 + 0.5) * h[0] - L[0] * center[0];
                    double y = (istart[1] + i1 + 0.5) * h[1] - L[1] * center[1];
                    double z = (istart[2] + i2 + 0.5) * h[2] - L[2] * center[2];
                    double rho2 = (x * x + y * y + z * z) * oosigma2;
                    double rho = sqrt(rho2);
                    const size_t id = localindex_xyz(i0,i1,i2,topo);
                    // Gaussian
                    rhs[id] = +c_1o4pi * oosigma3 * sqrt(2.0 / M_PI) * exp(-rho2 * 0.5);
                    sol[id] = +c_1o4pi * oosigma * 1.0 / rho * erf(rho * c_1osqrt2);
                }
            }
        }
        // read the source term and the solution
        xmf_write(topo,"rhs","data");
        hdf5_write(topo,"rhs","data",rhs);

        //-------------------------------------------------------------------------
        /** - solve the equations */
        //-------------------------------------------------------------------------
        mysolver->solve(topo,rhs,rhs,UP_SRHS);

        //-------------------------------------------------------------------------
        /** - compute the error */
        //-------------------------------------------------------------------------
        double lerr2 = 0.0;
        double lerri = 0.0;
        double gap = 0.0;

        /**
         * @todo change that to axis-based loops
         */
        for (int i2 = 0; i2 < topo->nloc(2); i2++)
        {
            for (int i1 = 0; i1 < topo->nloc(1); i1++)
            {
                for (int i0 = 0; i0 < topo->nloc(0); i0++)
                {
                    const size_t id = localindex_xyz(i0,i1,i2,topo);
                    const double err = sol[id] - rhs[id];

                    lerri = max(lerri, abs(err));
                    lerr2 += (err * err) * h[0] * h[1];
                }
            }
        }
        double erri = 0.0;
        double err2 = 0.0;
        MPI_Allreduce(&lerr2, &err2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&lerri, &erri, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        err2 = sqrt(err2);

        if (rank == 0)
        {
            FILE *myfile = fopen(filename, "a+");
            if (myfile != NULL)
            {
                fprintf(myfile, "%d %12.12e %12.12e\n", nglob[0], err2, erri);
                fclose(myfile);
            }
            else
            {
                UP_CHECK1(false, "unable to open file %s", filename);
            }
        }

        fftw_free(sol);
        fftw_free(rhs);
        // delete (mysolver);
        delete (topo);
    }
}
