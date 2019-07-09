#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include "fftw3.h"

#ifndef TOOLS_HPP
#define TOOLS_HPP

using namespace std;

// write the binary for the file
void write_array(const int n[2], const double* array,       const string array_name);
void write_array(const int n[2], const fftw_complex* array, const string array_name);

// write the content
void write_info        (const int nx,const int ny, const string array_name);
void write_data_real   (const int nx,const int ny, const double* array, const string array_name);
void write_data_complex(const int nx,const int ny, const fftw_complex* array, const string array_name);

#endif