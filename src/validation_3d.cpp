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

#include "validation_2d.hpp"

/**
 * @brief 
 * 
 * @param type 
 * @param orderdiff 
 */
void validation_3d_UU_UU(const int nsample, const int* size, const SolverType type, const OrderDiff orderdiff)
{

    // int size[2] = {1024,2048};
    // int size[4] = {64,128,256,512};
    // int size[4] = {64,128,256,512};
    // int size[1] = {64};

    char filename[512];
    sprintf(filename,"./data/%s-type=%d-orderdiff=%d.error",__func__,type,orderdiff);
    FILE* myfile = fopen(filename,"w+");
    if(myfile != NULL){
        fclose(myfile);
    }
    else{
        UP_CHECK1(false,"unable to open file %s",filename);
    }

    for(int is=0; is<nsample; is++){
        const int n[3] = {size[is],size[is],size[is]};

        const double L[3] = {2.0,2.0,2.0};
        const double h[3] = {L[0]/n[0],L[1]/n[1],L[2]/n[2]};

        const BoundaryType mybc[3][2] = {{UNB,UNB},{UNB,UNB},{UNB,UNB}};

        //-------------------------------------------------------------------------
        /** - Initialize the solver */
        //-------------------------------------------------------------------------
        FFTW_Solver*  mysolver = new FFTW_Solver(n,h,L,mybc);
        mysolver->set_GreenType(orderdiff);
        mysolver->setup(type);

        //-------------------------------------------------------------------------
        /** - allocate rhs and solution */
        //-------------------------------------------------------------------------
        double** rhs = (double**) fftw_malloc(sizeof(double*)*n[1]);
        double** sol = (double**) fftw_malloc(sizeof(double*)*n[1]);

        rhs[0] = (double*) fftw_malloc(sizeof(double)*n[0]*n[1]);
        sol[0] = (double*) fftw_malloc(sizeof(double)*n[0]*n[1]);
        
        for(int iy=1; iy<n[1] ; iy++){
            rhs[iy] = rhs[iy-1] + n[0];
            sol[iy] = sol[iy-1] + n[0];
        }

        //-------------------------------------------------------------------------
        /** - fill the rhs and the solution */
        //-------------------------------------------------------------------------
        const double sigma      = 0.05;
        const double oosigma    = 1.0/(sigma);
        const double oosigma2   = 1.0/(sigma*sigma);
        const double oosigma3   = 1.0/(sigma*sigma*sigma);
        const double center[3]  = {0.5,0.5,0.5};

        for(int iz=0; iz<n[2]; iz++){
            for(int iy=0; iy<n[1]; iy++){
                for(int ix=0; ix<n[0]; ix++){
                    rhs[iy][ix] = 0.0;
                    sol[iy][ix] = 0.0;

                    double x    = (ix+0.5)*h[0]-L[0]*center[0];
                    double y    = (iy+0.5)*h[1]-L[1]*center[1];
                    double z    = (iy+0.5)*h[2]-L[2]*center[2];
                    double rho2 = (x*x+y*y+z*z)*oosigma2;
                    double rho  = sqrt(rho2);

                    // Gaussian
                    rhs[iy][ix] += + c_1o4pi * oosigma3 * sqrt(2.0/M_PI) * exp(-rho2*0.5);
                    sol[iy][ix] += + c_1o4pi * oosigma  * 1.0/rho * erf(rho*c_1osqrt2);
                }
            }
        }

        //-------------------------------------------------------------------------
        /** - solve the equations */
        //-------------------------------------------------------------------------
        mysolver->solve(rhs[0],rhs[0]);

        //-------------------------------------------------------------------------
        /** - compute the error */
        //-------------------------------------------------------------------------
        double err2 = 0.0;
        double erri = 0.0;
        double gap  = 0.0;
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                const double err = sol[iy][ix]-rhs[iy][ix];

                erri  = max(erri,abs(err));
                err2 += (err*err) * h[0]*h[1];
            }
        }
        err2 = sqrt(err2);

        
        FILE* myfile = fopen(filename,"a+");
        if(myfile != NULL){
            fprintf(myfile,"%d %12.12e %12.12e\n",n[0],err2,erri);
            fclose(myfile);
        }
        else{
            UP_CHECK1(false,"unable to open file %s",filename);
        }


        fftw_free(sol[0]);
        fftw_free(sol);
        fftw_free(rhs[0]);
        fftw_free(rhs);
        delete(mysolver);
    }
}
