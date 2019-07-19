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


int main(int argc, char* argv[])
{
    /*
    This function solves the poisson equation nabla^2 u = f
    -  We assume a cell-centered layout
    */

    const int    n[2] = {16,16};
    const double L[2] = {1.0,1.0};
    const double h[2] = {L[0]/n[0] , L[1]/n[1]};

    // BoundaryType mybc[DIM][2] = {{EVEN,ODD},{UNB,EVEN}};
    BoundaryType mybc[DIM][2] = {{UNB,UNB},{EVEN,UNB}};
    // BoundaryType mybc[DIM][2] = {{UNB,EVEN},{UNB,UNB}};
    // BoundaryType mybc[DIM][2] = {{UNB,UNB},{UNB,UNB}};

    FFTW_Solver*  mysolver = new FFTW_Solver(n,h,L,mybc);
    mysolver->set_GreenType(HEJ_2);
    mysolver->setup(UP_SRHS);


    //-------------------------------------------------------------------------------
    // allocate the memory
    //-------------------------------------------------------------------------------
    
    // real arrays
    double** freal = (double**) fftw_malloc(sizeof(double*)*n[1]);
    double** ureal = (double**) fftw_malloc(sizeof(double*)*n[1]);

    freal[0] = (double*) fftw_malloc(sizeof(double)*n[0]*n[1]);
    ureal[0] = (double*) fftw_malloc(sizeof(double)*n[0]*n[1]);
    
    for(int iy=1; iy<n[1] ; iy++){
        freal[iy] = freal[iy-1] + n[0];
        ureal[iy] = ureal[iy-1] + n[0];
    }

    // extended arrays
    //sizes - mix in X, unbounded in Y
    const int n_green[2] = {2*n[0]+1,2*n[1]};
    double** greal_ext = (double**) fftw_malloc(sizeof(double*)*n_green[1]); // Green's function
    greal_ext[0] = (double*) fftw_malloc(sizeof(double)*n_green[0]*n_green[1]);
    for(int iy=1; iy<n_green[1]; iy++)
    {
        greal_ext[iy] = greal_ext[iy-1] + n_green[0];
    }

    const int n_ext  [2] = {2*n[0]  ,2*n[1]};
    double** freal_ext = (double**) fftw_malloc(sizeof(double*)*n_ext[1]); // only one extended field is needed
    freal_ext[0] = (double*) fftw_malloc(sizeof(double)*n_ext[0]*n_ext[1]);
    for(int iy=1; iy<n_ext[1] ; iy++){
        freal_ext[iy] = freal_ext[iy-1] + n_ext[0];
        
    }
    
    //Fourier arrays (kx,y)
    double* g_kx_y = (double*) fftw_malloc(sizeof(double)*n_green[0]*n_green[1]); // Green's function (kx,ky)
    double* f_kx_y = (double*) fftw_malloc(sizeof(double)*n_ext[0]*n_ext[1]); // field Fourier transfrom

    // Fourier arrays (kx,ky)
    const int n_ghat[2] = {n_green[0],n_green[1]/2+1};
    fftw_complex** ghat = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*)*n_ghat[1]); // Green's function (kx,ky)
    ghat[0] =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n_ghat[0]*n_ghat[1]); // Green's function (kx,ky)
    for(int iy=1; iy<n_ghat[1] ; iy++){
        ghat[iy] = ghat[iy-1] + n_ghat[0];
    }

    const int n_fhat[2] = {n_ext[0],n_ext[1]/2+1};
    fftw_complex** fhat = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*)*n_fhat[1]); // Green's function (kx,ky)
    fhat[0] =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n_fhat[0]*n_fhat[1]); // Green's function (kx,ky)
    for(int iy=1; iy<n_fhat[1]; iy++){
        fhat[iy] = fhat[iy-1] + n_fhat[0];
    }

    //-------------------------------------------------------------------------------
    // DCT/DST plans
    //-------------------------------------------------------------------------------
    fftw_plan green_r2f[2];
    fftw_plan field_r2f[2];
    fftw_plan field_f2r[2];
    fftw_plan green_f2r[2];

    double fact_green = 1.0;
    double fact_field = 1.0;

    // X - direction even - unbounded
    int rank    = 1; // since we do 1D transforms
    int howmany = n_green[1]; // we do n_ext[1] transforms
    int size    = n_green[0]; // size of each of the howmany transforms
    int idist   = n_green[0]; // input of the kth transform is at in+k*idist
    int odist   = n_green[0];
    int istride = 1; // jth element of the kth transfrom is located at k*idist+j*istride
    int ostride = 1;
    fftw_r2r_kind kind = FFTW_REDFT00;
    // Green needs a DCT type I
    fact_green *= 2.0*(size-1.0);
    green_r2f[0] = fftw_plan_many_r2r(rank,&size,howmany,greal_ext[0],NULL,istride,idist,g_kx_y,NULL,ostride,odist,&kind,FFTW_MEASURE);
    green_f2r[0] = fftw_plan_many_r2r(rank,&size,howmany,g_kx_y,NULL,ostride,odist,greal_ext[0],NULL,istride,idist,&kind,FFTW_MEASURE);

    // Field needs a DCT type II 
    howmany = n_ext[1]; // we do n_ext[1] transforms
    size    = n_ext[0]; // size of each of the howmany transforms
    idist   = n_ext[0]; // input of the kth transform is at in+k*idist
    odist   = n_ext[0];
    kind    = FFTW_REDFT10;
    fact_field *= 2.0*size;
    field_r2f[0] = fftw_plan_many_r2r(rank,&size,howmany,freal_ext[0],NULL,istride,idist,f_kx_y,NULL,ostride,odist,&kind,FFTW_MEASURE);
    // inverse of type II is type III
    kind = FFTW_REDFT01;
    field_f2r[0] = fftw_plan_many_r2r(rank,&size,howmany,f_kx_y,NULL,ostride,odist,freal_ext[0],NULL,istride,idist,&kind,FFTW_MEASURE);

    // Y - direction unbounded - unbounded
    rank    = 1; // since we do 1D transforms
    idist   = 1;        // input of the kth transform is at in+k*idist
    odist   = 1;
    howmany = n_green[0]; // we do n_ext[0] transforms
    size    = n_green[1]; // size of each of the howmany transforms
    istride = n_green[0]; // jth element of the kth transfrom is located at k*idist+j*istride
    ostride = n_ghat[0];
    // Green needs a real 2 complex transform
    fact_green *= size;
    green_r2f[1] = fftw_plan_many_dft_r2c(rank,&size,howmany,g_kx_y,NULL,istride,idist,ghat[0],NULL,ostride,odist,FFTW_MEASURE);
    green_f2r[1] = fftw_plan_many_dft_c2r(rank,&size,howmany,ghat[0],NULL,ostride,odist,g_kx_y,NULL,istride,idist,FFTW_MEASURE);

    rank    = 1; // since we do 1D transforms
    idist   = 1;        // input of the kth transform is at in+k*idist
    odist   = 1;
    howmany = n_ext[0]; // we do n_ext[0] transforms
    size    = n_ext[1]; // size of each of the howmany transforms
    istride = n_ext[0]; // jth element of the kth transfrom is located at k*idist+j*istride
    ostride = n_fhat[0];
    fact_field *= size;
    field_r2f[1] = fftw_plan_many_dft_r2c(rank,&size,howmany,f_kx_y,NULL,istride,idist,fhat[0],NULL,ostride,odist,FFTW_MEASURE);
    field_f2r[1] = fftw_plan_many_dft_c2r(rank,&size,howmany,fhat[0],NULL,ostride,odist,f_kx_y,NULL,istride,idist,FFTW_MEASURE);

    //-------------------------------------------------------------------------------
    // Compute some usefull factors
    //-------------------------------------------------------------------------------
    const double oo2pi   = 1.0/(2.0*M_PI);
    const double oo4pi   = 1.0/(4.0*M_PI);
    const double eps     = 2.0*h[0];
    const double oo2eps2 = 1.0/(2.0*eps*eps);

    //-------------------------------------------------------------------------------
    // Compute the Green's function nabla^2 G = + dirac
    //-------------------------------------------------------------------------------
    for(int iy=0 ; iy < (n[1]+1) ; iy++){
        for(int ix=0 ; ix < n_green[0] ; ix++){
            double x  = ix*h[0];
            double y  = iy*h[1];
            double r2 = x*x+y*y;
            
            greal_ext[iy][ix] = oo4pi * (log(r2) + expint_ei(r2*oo2eps2));
        }
    }
    greal_ext[0][0] = - oo2pi*(GAMMA*0.5 - log(sqrt(2.0)*eps));
    // do the symmetry in Y
    for(int iy=(n[1]+1); iy<n_green[1]; iy++ ){
        for(int ix=0 ; ix<n_green[0] ; ix++){
            // we have to take the symmetry around n: n - (iy - n)
            greal_ext[iy][ix] = greal_ext[2*n[1]-iy][ix];
        }
    }
    
    // for(int iy=0; iy<n_green[1]; iy++ ){
    //     for(int ix=0 ; ix<n_green[0] ; ix++){
    //         double x  = ix*h[0];
    //         double y  = iy*h[1];
    //         double r2 = x*x+y*y;
    //         greal_ext[iy][ix] = cos(10.0 * 2.0*M_PI/(2.0*L[0]) * y) * cos(3.0 * 2.0*M_PI/(2.0*L[0]) * x);
    //     }
    // }

    //-------------------------------------------------------------------------------
    // Compute the source term
    //-------------------------------------------------------------------------------
    for(int iy=0; iy<n_ext[1]; iy++){
        for(int ix=0; ix<n_ext[0]; ix++){
            freal_ext[iy][ix] = 0.0;
        }
    }
   const double sigma    = 0.1;
   const double oosigma2 = 1.0/(sigma*sigma);
    for(int iy=0; iy<n[1]; iy++){
        for(int ix=0; ix<n[0]; ix++){
            freal_ext[iy][ix] = 0.0;
            

            double x  = (ix+0.5)*h[0]-L[0]*0.5;
            double y  = (iy+0.5)*h[1]-L[1]*0.35;
            double r2 = x*x+y*y;

            freal_ext[iy][ix] += + oo2pi * oosigma2 * exp(-r2*oosigma2*0.5);


            x  = (ix+0.5)*h[0]-L[0]*0.5;
            y  = (iy+0.5)*h[1]-L[1]*0.65;
            r2 = x*x+y*y;

            freal_ext[iy][ix] += - oo2pi * oosigma2 * exp(-r2*oosigma2*0.5);

            freal[iy][ix] = freal_ext[iy][ix];

        }
    }

    mysolver->solve(freal[0],freal[0]);
    delete(mysolver);

    int mysize[2] = {n[0],n[1]};
    write_array(mysize,freal[0],"field_real");
    write_array(n_green,greal_ext[0],"green_real");
    

    //-------------------------------------------------------------------------------
    // Do the FFT FORWARD
    //-------------------------------------------------------------------------------
    // Take the Fourier transfrom
    fftw_execute(green_r2f[0]);
    fftw_execute(field_r2f[0]);
    write_array(n_green,g_kx_y,"green_kxy");
    write_array(n_ext,f_kx_y,"field_kxy");

    // printf("last value of g_kx_y = %f\n",g_kx_y[);
    fftw_execute(green_r2f[1]);
    fftw_execute(field_r2f[1]);

    write_array(n_ghat,ghat[0],"green_fourier");

   
    const double vol = h[0]*h[1];
    for(int iy=0; iy<n_fhat[1];  iy++){
        for(int ix=0; ix<n_fhat[0]; ix++){
            fhat[iy][ix][0] = fhat[iy][ix][0]/(fact_field);
            fhat[iy][ix][1] = fhat[iy][ix][1]/(fact_field);

            fhat[iy][ix][0] = vol * (ghat[iy][ix][0]*fhat[iy][ix][0] - ghat[iy][ix][1]*fhat[iy][ix][1]);
            fhat[iy][ix][1] = vol * (ghat[iy][ix][0]*fhat[iy][ix][1] + ghat[iy][ix][1]*fhat[iy][ix][0]);
            // double tempr = ghat[iy][ix][0]*fhat[iy][ix][0] - ghat[iy][ix][1]*fhat[iy][ix][1];
            // double tempi = ghat[iy][ix][0]*fhat[iy][ix][1] + ghat[iy][ix][1]*fhat[iy][ix][0];

        }
    }
    write_array(n_fhat,fhat[0],"field_fourier");
    

    // do the inverse one
    fftw_execute(field_f2r[1]);
    fftw_execute(field_f2r[0]);

    fftw_execute(green_f2r[1]);
    fftw_execute(green_f2r[0]);

    write_array(n_ext, freal_ext[0],"field_final");
    write_array(n_green, greal_ext[0],"green_final");


    // compute the error
    // setup the source term
    // for(int iy=0; iy<n_ext[1]; iy++){
    //     for(int ix=0; ix<n_ext[0]; ix++){
    //         freal_ext[iy][ix] = 0.0;
    //     }
    // }
    double err2 = 0.0;
    double erri = 0.0;
    double gap  = 0.0;
    for(int iy=0; iy<n[1]; iy++){
        for(int ix=0; ix<n[0]; ix++){
            double solution = 0.0;
            double x  = (ix+0.5)*h[0]-L[0]*0.5;
            double y  = (iy+0.5)*h[1]-L[1]*0.35;
            double r2 = x*x+y*y;
            double x_sym = x+L[0];
            double y_sym = y;
            double r2_sym = x_sym*x_sym + y_sym*y_sym;

            // solution += 1.0/(4.0*M_PI) * (log(r2    /rhs_sigma2)+ expint_ei(r2    /rhs_sigma2));
            solution += oo4pi * (log(r2    *oosigma2*0.5)+ expint_ei(r2    *oosigma2*0.5));
            solution += oo4pi * (log(r2_sym*oosigma2*0.5)+ expint_ei(r2_sym*oosigma2*0.5));


            x  = (ix+0.5)*h[0]-L[0]*0.5;
            y  = (iy+0.5)*h[1]-L[1]*0.65;
            r2 = x*x+y*y;
            x_sym = x+L[0];
            y_sym = y;
            r2_sym = x_sym*x_sym + y_sym*y_sym;

            // solution += - 1.0/(4.0*M_PI) * (log(r2    /rhs_sigma2)+ expint_ei(r2    /rhs_sigma2));
            // solution += - 1.0/(4.0*M_PI) * (log(r2_sym/rhs_sigma2)+ expint_ei(r2_sym/rhs_sigma2));
            solution += -oo4pi * (log(r2    *oosigma2*0.5)+ expint_ei(r2    *oosigma2*0.5));
            solution += -oo4pi * (log(r2_sym*oosigma2*0.5)+ expint_ei(r2_sym*oosigma2*0.5));

            
            erri = max(erri,abs(solution-freal_ext[iy][ix]));

            freal_ext[iy][ix] = solution;
        }
    }

    printf("erri = %e\n",erri);

    write_array(n_ext, freal_ext[0],"field_sol");

    //-------------------------------------------------------------------------------
    // free the memory
    //-------------------------------------------------------------------------------
    fftw_free(freal[0]);
    fftw_free(ureal[0]);
    fftw_free(greal_ext[0]);
    fftw_free(freal_ext[0]);

    fftw_free(freal);
    fftw_free(ureal);
    fftw_free(greal_ext);
    fftw_free(freal_ext);

    fftw_destroy_plan(green_r2f[0]);
    fftw_destroy_plan(green_r2f[1]);
    fftw_destroy_plan(field_r2f[0]);
    fftw_destroy_plan(field_r2f[1]);
    fftw_destroy_plan(field_f2r[0]);
    fftw_destroy_plan(field_f2r[1]);

    fftw_cleanup();
}