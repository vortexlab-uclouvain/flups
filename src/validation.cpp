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

#include "validation.hpp"

/**
 * @brief 
 * 
 * @param type 
 * @param orderdiff 
 */
void validation_2d_UU_UU(const SolverType type, const OrderDiff orderdiff){

    // int size[2] = {1024,2048};
    int size[4] = {64,128,256,512};
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

    for(int is=0; is<4; is++){
        const int n[2] = {size[is],size[is]};

        // printf("doing size = ")

        const double L[2] = {2.0,2.0};
        const double h[2] = {L[0]/n[0],L[1]/n[1]};

        const BoundaryType mybc[2][2] = {{UNB,UNB},{UNB,UNB}};

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
        const double sigma    = 0.05;
        const double oosigma2 = 1.0/(sigma*sigma);
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                rhs[iy][ix] = 0.0;
                sol[iy][ix] = 0.0;

                double x    = (ix+0.5)*h[0]-L[0]*0.5;
                double y    = (iy+0.5)*h[1]-L[1]*0.35;
                double rho2 = (x*x+y*y)*oosigma2;

                // High order algabraic
                // rhs[iy][ix] += + c_1o2pi * oosigma2 * 4.0/pow(rho2+1.0,3);
                // sol[iy][ix] += + c_1o4pi * (log(rho2+1.0)+ rho2/(rho2+1.0));
                // Gaussian
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));

                x    = (ix+0.5)*h[0]-L[0]*0.5;
                y    = (iy+0.5)*h[1]-L[1]*0.65;
                rho2 = (x*x+y*y)*oosigma2;

                // High order algabraic
                // rhs[iy][ix] += - c_1o2pi * oosigma2 * 4.0/pow(rho2+1.0,3);
                // sol[iy][ix] += - c_1o4pi * (log(rho2+1.0)+ rho2/(rho2+1.0));
                // Gaussian
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
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

/**
 * @brief 
 * 
 * @param type 
 * @param orderdiff 
 */
void validation_2d_UU_UE(const SolverType type, const OrderDiff orderdiff){

    // int size[2] = {1024,2048};
    int size[4] = {64,128,256,512};
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

    for(int is=0; is<4; is++){
        const int n[2] = {size[is],2*size[is]};

        // printf("doing size = ")

        const double L[2] = {1.0,2.0};
        const double h[2] = {L[0]/n[0],L[1]/n[1]};

        const BoundaryType mybc[2][2] = {{UNB,UNB},{UNB,EVEN}};

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
        const double sigma    = 0.05;
        const double oosigma2 = 1.0/(sigma*sigma);
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                rhs[iy][ix] = 0.0;
                sol[iy][ix] = 0.0;

                double center[2];
                //----------------------------------------------
                center[0] = 0.5; center[1]=0.35;
                double x    = (ix+0.5)*h[0]-L[0]*center[0];
                double y    = (iy+0.5)*h[1]-L[1]*center[1];
                double rho2 = (x*x+y*y)*oosigma2;
                double x_sym = (ix+0.5)*h[0]-L[0]*center[0];
                double y_sym = (iy+0.5)*h[1]-L[1]*(2.0-center[1]);
                double rho2_sym = (x_sym*x_sym + y_sym*y_sym)*oosigma2;
                // Gaussian
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2_sym*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2_sym*0.5)+ expint_ei(rho2_sym*0.5));

                //----------------------------------------------
                center[0] = 0.5; center[1]=0.65;
                x    = (ix+0.5)*h[0]-L[0]*center[0];
                y    = (iy+0.5)*h[1]-L[1]*center[1];
                rho2 = (x*x+y*y)*oosigma2;
                x_sym = (ix+0.5)*h[0]-L[0]*center[0];
                y_sym = (iy+0.5)*h[1]-L[1]*(2.0-center[1]);
                rho2_sym = (x_sym*x_sym + y_sym*y_sym)*oosigma2;
                // Gaussian
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2_sym*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2_sym*0.5)+ expint_ei(rho2_sym*0.5));
            }
        }

        //-------------------------------------------------------------------------
        /** - solve the equations */
        //-------------------------------------------------------------------------
        mysolver->solve(rhs[0],rhs[0]);

        //-------------------------------------------------------------------------
        /** - dump the array */
        //-------------------------------------------------------------------------

        //-------------------------------------------------------------------------
        /** - compute the error */
        //-------------------------------------------------------------------------
        double err2 = 0.0;
        double erri = 0.0;
        double gap  = 0.0;
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                erri  = max(erri,fabs(sol[iy][ix]-rhs[iy][ix]));
                err2 += pow(sol[iy][ix]-rhs[iy][ix],2) * h[0]*h[1];
            }
        }
        err2 = sqrt(err2);

        
        FILE* myfile = fopen(filename,"a+");
        if(myfile != NULL){
            fprintf(myfile,"%d %e %e\n",n[0],err2,erri);
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


/**
 * @brief 
 * 
 * @param type 
 * @param orderdiff 
 */
void validation_2d_UU_UO(const SolverType type, const OrderDiff orderdiff){

    // int size[2] = {1024,2048};
    // int size[4] = {64,128,256,512}; int nsample = 4;
    int size[2] = {128,256}; int nsample = 2;
    // int size[4] = {64,128,256,512}; int nsample = 4;
    // int size[1] = {64}; int nsample = 1;

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
        const int n[2] = {size[is],2*size[is]};

        // printf("doing size = ")

        const double L[2] = {1.0,2.0};
        const double h[2] = {L[0]/n[0],L[1]/n[1]};

        const BoundaryType mybc[2][2] = {{UNB,UNB},{UNB,ODD}};

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
        const double sigma    = 0.05;
        const double oosigma2 = 1.0/(sigma*sigma);
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                rhs[iy][ix] = 0.0;
                sol[iy][ix] = 0.0;

                double center[2];
                //----------------------------------------------
                center[0] = 0.5; center[1]=0.35;
                double x    = (ix+0.5)*h[0]-L[0]*center[0];
                double y    = (iy+0.5)*h[1]-L[1]*center[1];
                double rho2 = (x*x+y*y)*oosigma2;
                double x_sym = (ix+0.5)*h[0]-L[0]*center[0];
                double y_sym = (iy+0.5)*h[1]-L[1]*(2.0-center[1]);
                double rho2_sym = (x_sym*x_sym + y_sym*y_sym)*oosigma2;
                // Gaussian
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2_sym*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2_sym*0.5)+ expint_ei(rho2_sym*0.5));

                //----------------------------------------------
                center[0] = 0.5; center[1]=0.65;
                x    = (ix+0.5)*h[0]-L[0]*center[0];
                y    = (iy+0.5)*h[1]-L[1]*center[1];
                rho2 = (x*x+y*y)*oosigma2;
                x_sym = (ix+0.5)*h[0]-L[0]*center[0];
                y_sym = (iy+0.5)*h[1]-L[1]*(2.0-center[1]);
                rho2_sym = (x_sym*x_sym + y_sym*y_sym)*oosigma2;
                // Gaussian
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2_sym*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2_sym*0.5)+ expint_ei(rho2_sym*0.5));
            }
        }

        //-------------------------------------------------------------------------
        /** - solve the equations */
        //-------------------------------------------------------------------------
        mysolver->solve(rhs[0],rhs[0]);

        //-------------------------------------------------------------------------
        /** - dump the array */
        //-------------------------------------------------------------------------
        write_array(n,rhs[0],"sol");
        write_array(n,sol[0],"anal");

        //-------------------------------------------------------------------------
        /** - compute the error */
        //-------------------------------------------------------------------------
        double err2 = 0.0;
        double erri = 0.0;
        double gap  = 0.0;
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                erri  = max(erri,fabs(sol[iy][ix]-rhs[iy][ix]));
                err2 += pow(sol[iy][ix]-rhs[iy][ix],2) * h[0]*h[1];
            }
        }
        err2 = sqrt(err2);

        
        FILE* myfile = fopen(filename,"a+");
        if(myfile != NULL){
            fprintf(myfile,"%d %e %e\n",n[0],err2,erri);
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

/**
 * @brief 
 * 
 * @param type 
 * @param orderdiff 
 */
void validation_2d_UU_OU(const SolverType type, const OrderDiff orderdiff){

    int size[2] = {128,256}; int nsample = 2;
    // int size[4] = {64,128,256,512}; int nsample = 4;
    // int size[4] = {64,128,256,512}; int nsample = 4;
    // int size[1] = {64}; int nsample = 1;

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
        const int n[2] = {size[is],2*size[is]};

        // printf("doing size = ")

        const double L[2] = {1.0,2.0};
        const double h[2] = {L[0]/n[0],L[1]/n[1]};

        const BoundaryType mybc[2][2] = {{UNB,UNB},{ODD,UNB}};

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
        const double sigma    = 0.05;
        const double oosigma2 = 1.0/(sigma*sigma);
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                rhs[iy][ix] = 0.0;
                sol[iy][ix] = 0.0;

                double center[2];
                //----------------------------------------------
                center[0] = 0.5; center[1]=0.35;
                double x    = (ix+0.5)*h[0]-L[0]*center[0];
                double y    = (iy+0.5)*h[1]-L[1]*center[1];
                double rho2 = (x*x+y*y)*oosigma2;
                double x_sym = (ix+0.5)*h[0]-L[0]*center[0];
                double y_sym = (iy+0.5)*h[1]-L[1]*(-center[1]);
                double rho2_sym = (x_sym*x_sym + y_sym*y_sym)*oosigma2;
                // Gaussian
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2_sym*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2_sym*0.5)+ expint_ei(rho2_sym*0.5));

                //----------------------------------------------
                center[0] = 0.5; center[1]=0.65;
                x    = (ix+0.5)*h[0]-L[0]*center[0];
                y    = (iy+0.5)*h[1]-L[1]*center[1];
                rho2 = (x*x+y*y)*oosigma2;
                x_sym = (ix+0.5)*h[0]-L[0]*center[0];
                y_sym = (iy+0.5)*h[1]-L[1]*(-center[1]);
                rho2_sym = (x_sym*x_sym + y_sym*y_sym)*oosigma2;
                // Gaussian
                rhs[iy][ix] += - c_1o2pi * oosigma2 * exp(-rho2*0.5);
                sol[iy][ix] += - c_1o4pi * (log(rho2*0.5)+ expint_ei(rho2*0.5));
                rhs[iy][ix] += + c_1o2pi * oosigma2 * exp(-rho2_sym*0.5);
                sol[iy][ix] += + c_1o4pi * (log(rho2_sym*0.5)+ expint_ei(rho2_sym*0.5));
            }
        }

        //-------------------------------------------------------------------------
        /** - solve the equations */
        //-------------------------------------------------------------------------
        mysolver->solve(rhs[0],rhs[0]);

        //-------------------------------------------------------------------------
        /** - dump the array */
        //-------------------------------------------------------------------------
        write_array(n,rhs[0],"sol");
        write_array(n,sol[0],"anal");

        //-------------------------------------------------------------------------
        /** - compute the error */
        //-------------------------------------------------------------------------
        double err2 = 0.0;
        double erri = 0.0;
        double gap  = 0.0;
        for(int iy=0; iy<n[1]; iy++){
            for(int ix=0; ix<n[0]; ix++){
                erri  = max(erri,fabs(sol[iy][ix]-rhs[iy][ix]));
                err2 += pow(sol[iy][ix]-rhs[iy][ix],2) * h[0]*h[1];
            }
        }
        err2 = sqrt(err2);

        
        FILE* myfile = fopen(filename,"a+");
        if(myfile != NULL){
            fprintf(myfile,"%d %e %e\n",n[0],err2,erri);
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