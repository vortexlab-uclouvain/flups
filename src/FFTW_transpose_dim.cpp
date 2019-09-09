// /**
//  * @file FFTW_transpose_dim.cpp
//  * @author Thomas Gillis
//  * @brief 
//  * @version
//  * @date 2019-07-23
//  * 
//  * @copyright Copyright © UCLouvain 2019
//  * 
//  */


// #include "FFTW_transpose_dim.hpp"


// FFTW_transpose_dim::FFTW_transpose_dim(int global_size[3],const int nproc)
// {
//     BEGIN_FUNC;

//     size_t rank;
//     MPI_Comm_rank(MPI_COMM_WORLD,(int*) &rank);

//     // compute pencil start and pencil end for each dim
//     for(int id=0; id<3; id++){
//         // data is aligned in the id direction
//         _size[id][id] = global_size[id];

//         size_t sizemult=1;
//         for(int sid=0;      sid<id; sid++) sizemult *= global_size[sid];
//         for(int sid=id+1;   sid<3;sid++) sizemult *= global_size[sid];

//         // divide the sizemult by nproc - we prefer that one proc waits for everybody
//         // if I can remove one proc and obtain the correct number, I overload the last proc
//         // if not, I underload it
//         int npencil = (sizemult%(nproc-1) == 0) ? floor(sizemult/nproc) : ceil(sizemult/nproc);
//         _pencil_start[id] = (rank*npencil);
//         _pencil_end  [id] = ((rank+1) == nproc)? sizemult : (rank+1)*npencil ;


//     }

//     // compute how to go to the next one
//     for(int id=0; id<3; id++){
//         const int next = (id+1)%3;
//         // each element in a pencil in the next dimension corresponds to global_size[id] memory
//         size_t mem_min = _pencil_start[next]*global_size[next] * global_size[id];
//         size_t mem_max = _pencil_end  [next]*global_size[next] * global_size[id];

//         size_t mymem = mem_max - mem_min;
        

//     }

//     // compute the interaction list
//     // from O to 1
//     int size_self = 1;
//     int size_self = 1;

//     while()
    

// }