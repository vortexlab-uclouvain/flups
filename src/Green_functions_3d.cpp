/**
 * @file Green_functions_3d.cpp
 * @author Denis-Gabriel Caprace, Thomas Gillis
 * @brief 
 * @version
 * @date 2019-07-22
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * -------------------------------
 * **Symmetry computation:**
 * 
 * We have to take the symmetry around symstart. e.g. in X direction: `symstart[0] - (ix - symstart[0]) = 2 symstart[0] - ix`
 * 
 * In some cases when we have an R2C transform, it ask for 2 additional doubles.
 * The value is meaningless but we would like to avoid segfault and nan's.
 * To do so, we use 2 tricks:
 * - The `abs` is used to stay on the positivie side and hence avoid negative memory access
 * - The `max` is used to prevent the computation of the value in 0, which is never used in the symmetry.
 * 
 * As an example, the final formula is then ( in the X direction):
 * `max( abs(2 symstart[0] - ix) , 1)`
 * 
 */

#include "Green_functions_3d.hpp"

void Green_3D_3dirunbounded_0dirspectral_chat2(const Topology *topo, const double hfact[3], double *green)
{
    /**
     * @todo get the 3dir unbounded green function
     */

    UP_CHECK0(false, "non implemented function");
}


/**
 * @brief generic type for Green kernel, takes a table of parameters that can be used depending on the kernel
 * 
 */
typedef double (*GreenKernel)(const void* );

//notice that these function will likely not be inlined as we have a pointer to them...
static inline double _hej_2(const void* params) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    return c_1o4pi / r * (erf(r / eps * c_1osqrt2));
}
static inline double _hej_4(const void* params) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (rho)*exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double _hej_6(const void* params) {
    double r   = ((double*)params) [0];
    double eps = ((double*)params) [1];
    double rho = r / eps;
    return c_1o4pi / r * (c_1osqrt2 * c_1osqrtpi * (c_7o4 * rho - c_1o4 * pow(rho, 3)) * exp(-rho * rho * .5 ) + erf(rho * c_1osqrt2));
}
static inline double _chat_2(const void* params) {
    double r   = ((double*)params) [0];
    return c_1o4pi / r ;
}

/**
 * @brief Compute the Green function for 3dirunbounded
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors
 * @param symstart the symmetry plan asked
 * @param green the Green function array
 * @param typeGreen the type of Green function to be 
 * @param alpha the smoothing parameter (only used for HEJ kernels)
 * 
 */
void Green_3D_3dirunbounded_0dirspectral(const Topology *topo, const double hfact[3], const int symstart[3], double *green, GreenType typeGreen, const double alpha)
{
    BEGIN_FUNC

    UP_CHECK0(!(topo->isComplex()),"Green topology cannot been complex");

    // assert that the green spacing is not 0.0 everywhere
    UP_CHECK0(hfact[0] != 0.0, "grid spacing cannot be 0");
    UP_CHECK0(hfact[1] != 0.0, "grid spacing cannot be 0");
    UP_CHECK0(hfact[2] != 0.0, "grid spacing cannot be 0");

    int ax0 = topo->axis();
    int ax1 = (ax0 + 1) % 3;
    int ax2 = (ax0 + 2) % 3;

    const double eps     = alpha * hfact[0];
    
    double      G0;  //value of G in 0
    GreenKernel G;

    switch (typeGreen) {
        case HEJ_2:
            G  = &_hej_2;
            G0 =       M_SQRT2 / (4.0 * eps * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_4:
            G  = &_hej_4;
            G0 = 3.0 * M_SQRT2 / (8.0 * eps * sqrt(M_PI * M_PI * M_PI));
            break;
        case HEJ_6:
            G  = &_hej_6;
            G0 = 15.0 * M_SQRT2 / (32.0 * eps * sqrt(M_PI * M_PI * M_PI));
            break;
        case CHAT_2:
            G  = &_chat_2;
            G0 = .5 * pow(1.5 * c_1o2pi * hfact[0] * hfact[1] * hfact[2], 2. / 3.);
            break;
        case LGF_2:
            UP_ERROR("Lattice Green Function not implemented yet.");
            //please add the parameters you need to params
            break;
        default:
            UP_ERROR("Green Function type unknow.");
    }

    int istart[3];
    get_idstart_glob(istart,topo);

    for (int i2 = 0; i2 < topo->nloc(ax2); i2++)
    {
        for (int i1 = 0; i1 < topo->nloc(ax1); i1++)
        {
            for (int i0 = 0; i0 < topo->nloc(ax0); i0++)
            {
                // exact indexes
                const int ie0 = (istart[ax0] + i0);
                const int ie1 = (istart[ax1] + i1);
                const int ie2 = (istart[ax2] + i2);

                // symmetrize indexes
                const int is0 = (symstart[ax0]==0 || ie0 <= symstart[ax0]) ? ie0 : std::max(abs(2*(int)symstart[ax0]-ie0),1);
                const int is1 = (symstart[ax1]==0 || ie1 <= symstart[ax1]) ? ie1 : std::max(abs(2*(int)symstart[ax1]-ie1),1);
                const int is2 = (symstart[ax2]==0 || ie2 <= symstart[ax2]) ? ie2 : std::max(abs(2*(int)symstart[ax2]-ie2),1);

                // symmetrized position
                const double x0 = (is0) * hfact[ax0];
                const double x1 = (is1) * hfact[ax1];
                const double x2 = (is2) * hfact[ax2];

                // green function value
                const double r2  = x0 * x0 + x1 * x1 + x2 * x2;
                const double r   = sqrt(r2);
                const size_t id  = i0 + topo->nloc(ax0) * (i1 + i2 * topo->nloc(ax1));
                
                const double tmp[2] = {r,eps};
                green[id] = - G( tmp );
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0)
    {
        green[0] = -G0;
    }
}