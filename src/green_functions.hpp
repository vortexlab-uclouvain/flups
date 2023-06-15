/**
 * @file green_functions.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#include "defines.hpp"
#include "Topology.hpp"
#include "bessel.hpp"
#include "expint.hpp"

// define macros to strigyfy, both are required!
#define STR(a) ZSTR(a)
#define ZSTR(a) #a


void cmpt_Green_3dirunbounded(const Topology *topo, const double hfact[3],                                                 const double symstart[3], double *green, GreenType typeGreen, const double length);
void cmpt_Green_2dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length);
void cmpt_Green_1dirunbounded(const Topology *topo, const double hfact[3], const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length);
void cmpt_Green_0dirunbounded(const Topology *topo, const double hgrid   , const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length);
void cmpt_Green_0dirunbounded(const Topology *topo, const double hgrid   , const double kfact[3], const double koffset[3], const double symstart[3], double *green, GreenType typeGreen, const double length, const int istart_custom[3], const int iend_custom[3]);

/**
 * @brief read the LGF file in the KERNEL_PATH folder
 * 
 * @param [in] greendim the dimension of the Green function to use, 2D or 3D
 * @param [out] N the size above which we switch to the approximation, i.e. the size of the pre-stored kernel is N^3
 * @param [out] data the data where we store the 
 */
static void lgf_readfile_(GreenType typeGreen, const int greendim, int* N, double** data) {
    BEGIN_FUNC;

    // get LGF order
    int order = 0;
    if (LGF_2 == typeGreen) {
        order = 2;
    } else if (LGF_4 == typeGreen || MEHR_4L == typeGreen || MEHR_4F == typeGreen) {
        order = 4;
    } else if (LGF_6 == typeGreen || MEHR_6L == typeGreen || MEHR_6F == typeGreen) {
        order = 6;
    } else if (LGF_8 == typeGreen) {
        order = 8;
    } else {
        FLUPS_CHECK(false, "LGF Kernel type not recognized");
    }

    // some defined parameters:
    char lgfname[512];
    char path[] = STR(KERNEL_PATH);
    int datasize = 0;
    if (greendim == 3) {
        (*N) = 32;
        datasize = (*N) * (*N) * (*N);
        if (MEHR_4L == typeGreen || MEHR_6L == typeGreen) {
            sprintf(lgfname, "%s/MEHR_%dL_3d_%d.ker", path, order, (*N));
        } else if (MEHR_4F == typeGreen || MEHR_6F == typeGreen) {
            sprintf(lgfname, "%s/MEHR_%dF_3d_%d.ker", path, order, (*N));
        } else {
            sprintf(lgfname, "%s/LGF_%d_3d_%d.ker", path, order, (*N));
        }
    } else if (greendim == 2) {
        (*N) = 32;
        datasize = (*N) * (*N);
        if (MEHR_4L == typeGreen || MEHR_6L == typeGreen || MEHR_4F == typeGreen || MEHR_6F == typeGreen) {
            FLUPS_CHECK(false, "MEHR kernels not implemented in 2D");
        } else {
            sprintf(lgfname, "%s/LGF_%d_2d_%d.ker", path, order, (*N));
        }
    } else {
        FLUPS_CHECK(false, "Greendim = %d is not available in this version", greendim);
    }

    // open the file
    FILE *lgf_file = fopen(lgfname, "r");
    // display the information to the user
    FLUPS_INFO_1("loading the LGF kernel function %s", lgfname);

    (*data) = NULL;
    // start to read the file
    if (lgf_file != NULL) {
        // allocate the data
        (*data) = (double *)m_calloc(sizeof(double) * datasize);
        fread((*data), sizeof(double), datasize, lgf_file);
        // close the file
        fclose(lgf_file);
    } else {
        FLUPS_CHECK(false, "unable to read file %s", lgfname);
    }
    END_FUNC;
}