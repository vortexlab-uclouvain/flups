
#include "tools.hpp"

using namespace std;


void write_array(const int n[2], const double* array,const string array_name){
    write_info     (n[0],n[1],array_name);
    write_data_real(n[0],n[1],array,array_name);
}
void write_array(const int n[2], const fftw_complex* array,const string array_name){
    // we transform the array into a doubled-size array
    write_info        (n[0],n[1],array_name);
    write_data_complex(n[0],n[1],array,array_name);
}

void write_info(const int nx,const int ny, const string array_name)
{
    FILE* infofile = NULL;
    string info_name = "./data/" + array_name + ".info";
    infofile = fopen(info_name.c_str(),"w+");
    if(infofile != NULL)
    {
        fprintf(infofile,"%d %d %ld",nx,ny,sizeof(double));
        fclose(infofile);
    }
    else
    {
        printf("WARNING couldn't open file %s!!",info_name.c_str());
    }
}

void write_data_real(const int nx,const int ny, const double* array, const string array_name)
{
    FILE* datafile = NULL;
    string data_name = "./data/" + array_name + ".data";
    datafile = fopen(data_name.c_str(),"wb+");
    if(datafile != NULL)
    {
        size_t mysize = nx*ny;
        fwrite(array,sizeof(double),mysize,datafile);
        fclose(datafile);
    }
    else
    {
        printf("WARNING couldn't open file %s!!",data_name.c_str());
    }
}
void write_data_complex(const int nx,const int ny, const fftw_complex* array, const string array_name)
{
    FILE* datafile_real = NULL;
    FILE* datafile_imag = NULL;
    string data_name_real = "./data/" + array_name + "_real.data";
    string data_name_imag = "./data/" + array_name + "_imag.data";
    datafile_real = fopen(data_name_real.c_str(),"wb+");
    datafile_imag = fopen(data_name_imag.c_str(),"wb+");

    if(datafile_real != NULL && datafile_imag != NULL)
    {
        size_t mysize = nx*ny;
        for(int i=0;i<mysize; i++)
        {
            fwrite(&(array[i][0]),sizeof(double),1,datafile_real);
            fwrite(&(array[i][1]),sizeof(double),1,datafile_imag);
        }
        fclose(datafile_real);
        fclose(datafile_imag);
    }
    else
    {
        printf("WARNING couldn't open file %s or %s!!",data_name_real.c_str(),data_name_imag.c_str());
    }
}