/**
 * @file green_functions_3d.cpp
 * @author Thomas Gillis
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

#include "green_functions_3d.hpp"

void Green_3D_3dirunbounded_0dirspectral_chat2(const Topology *topo, const double hfact[3], double *green)
{
    /**
     * @todo get the 3dir unbounded green function
     */

    UP_CHECK0(false, "non implemented function");
}

/**
 * @brief Compute the Green function of Hejlesen order 2
 * 
 * @param topo the topology associated to the Green's function
 * @param hfact the h multiplication factors
 * @param symstart the symmetry plan asked
 * @param green the Green function array
 * @param alpha the smoothing parameter
 * 
 */
void Green_3D_3dirunbounded_0dirspectral_hej2(const Topology *topo, const double hfact[3], const int symstart[3], double *green, const double alpha)
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

    const double eps = alpha * hfact[0];
    const double oo2eps2 = 1.0 / (2.0 * eps * eps);

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
                const double r = sqrt(x0 * x0 + x1 * x1 + x2 * x2);
                const size_t id = i0 + topo->nloc(ax0) * (i1 + i2 * topo->nloc(ax1));
                green[id] = -c_1o4pi * 1.0 / r * (erf(r * c_1osqrt2 / eps));
            }
        }
    }
    // reset the value in 0.0
    if (istart[ax0] == 0 && istart[ax1] == 0 && istart[ax2] == 0)
    {
        green[0] = -M_SQRT2 / (4.0 * eps * sqrt(M_PI * M_PI * M_PI));
    }
}

void Green_3D_3dirunbounded_0dirspectral_hej4(const Topology *topo, const double hfact[3], double *green, const double alpha)
{
    BEGIN_FUNC

    /**
     * @todo get the 3dir unbounded green function
     */

    UP_CHECK0(false, "non implemented function");
}