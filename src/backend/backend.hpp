#ifndef BACKEND_HPP
#define BACKEND_HPP

#if defined(FLUPS_FFTW_BACKEND)
#include "backend_fftw.hpp"
#elif defined(FLUPS_CUDA_BACKEND)
#include "backend_cuda.hpp"
#else
#error "Unknown/no backend defined. Please specify a correct backend !"
#endif

#endif