/**
 * @file dothemagic_rhs.ipp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief contains the domagic functions 
 * @version
 * 
 * @copyright Copyright © UCLouvain 2019
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 *  
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 *  
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

#if (KIND == 0)
/**
 * @brief perform the convolution for real to real cases
 * 
 */
void Solver::dothemagic_rhs_real(double *data) {
#elif (KIND == 1)
/**
 * @brief perform the convolution for complex to complex cases
 * 
 */
void Solver::dothemagic_rhs_complex(double *data) {
#endif

    BEGIN_FUNC;
    int cdim = _ndim - 1;  // get current dim
#if (KIND == 0)
    FLUPS_CHECK(_topo_hat[cdim]->nf() == 1, "The topo_hat[2] (field) has to be complex", LOCATION);
#else
    FLUPS_CHECK(_topo_hat[cdim]->nf() == 2, "The topo_hat[2] (field) has to be complex", LOCATION);
#endif
    // get the axis
    const int nf  = _topo_hat[cdim]->nf();
    const int ax0 = _topo_hat[cdim]->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    
    // get the norm factor
    const double         normfact = _normfact;

    // get the adresses
    opt_double_ptr       mydata   = data;
    const opt_double_ptr mygreen  = _green;
    FLUPS_CHECK(FLUPS_ISALIGNED(mygreen) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
    FLUPS_CHECK(FLUPS_ISALIGNED(mydata) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
    FLUPS_ASSUME_ALIGNED(mydata, FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mygreen, FLUPS_ALIGNMENT);

    // get the number of pencils for the field and green
    const size_t onmax = _topo_hat[cdim]->nloc(ax1) * _topo_hat[cdim]->nloc(ax2) * _topo_hat[cdim]->lda();
    const size_t ondim = _topo_hat[cdim]->nloc(ax1) * _topo_hat[cdim]->nloc(ax2);
    const size_t inmax = _topo_hat[cdim]->nloc(ax0);
    // get the memory details
    const size_t memdim  = _topo_hat[cdim]->memdim();
    const int    nmem[3] = {_topo_hat[cdim]->nmem(0), _topo_hat[cdim]->nmem(1), _topo_hat[cdim]->nmem(2)};
    
    // do the loop
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, ondim, inmax, memdim, nmem, mydata, mygreen, normfact, ax0, nf)
    for (size_t id = 0; id < onmax; id++) {
        // get the lia and the io index
        const size_t lia = id / ondim;
        const size_t io  = id % ondim;

        // get the starting pointer
        opt_double_ptr greenloc = mygreen + collapsedIndex(ax0, 0, io, nmem, nf);  //lda of Green is only 1
        opt_double_ptr dataloc  = mydata + lia * memdim + collapsedIndex(ax0, 0, io, nmem, nf);

        FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);
        FLUPS_ASSUME_ALIGNED(greenloc, FLUPS_ALIGNMENT);

        // do the actual convolution
        for (size_t ii = 0; ii < inmax; ii++) {
#if (KIND == 0)
            dataloc[ii] *= normfact * greenloc[ii];
#elif (KIND == 1)
            const double a = dataloc[ii * 2 + 0];
            const double b = dataloc[ii * 2 + 1];
            const double c = greenloc[ii * 2 + 0];
            const double d = greenloc[ii * 2 + 1];
            // update the values
            dataloc[ii * 2 + 0] = normfact * (a * c - b * d);
            dataloc[ii * 2 + 1] = normfact * (a * d + b * c);
#endif
        }
    }

    END_FUNC;
}
