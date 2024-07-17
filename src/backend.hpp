#ifndef BACKEND_HPP
#define BACKEND_HPP

#if defined (HAVE_CUDA)
#define M_HAVE_GPU 1
#include "backend_cuda.hpp"
#else
#define M_HAVE_GPU 0
#include "backend_fftw.hpp"
#endif

#endif  // BACKEND_HPP