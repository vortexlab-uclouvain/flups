/**
 * @file dothemagic_div.ipp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief contains the domagic functions for the rotational case
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
    void Solver::dothemagic_div_real_p1(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (KIND == -1)
    void Solver::dothemagic_div_real_m1(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (KIND == 1)
    void Solver::dothemagic_div_complex_p1(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (KIND == 2)
    void Solver::dothemagic_div_complex_m1(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (KIND == 3)
    void Solver::dothemagic_div_complex_pi(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (KIND == 4)
    void Solver::dothemagic_div_complex_mi(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#endif

        BEGIN_FUNC;
        int cdim = _ndim - 1;  // get current dim
#if (KIND == 0)
        FLUPS_CHECK(_topo_hat[cdim]->nf() == 1, "The topo_hat[2] (field) has to be complex", LOCATION);
#else
        FLUPS_CHECK(_topo_hat[cdim]->nf() == 2, "The topo_hat[2] (field) has to be complex", LOCATION);
#endif
        // get the axis
        const int ax0 = _topo_hat[cdim]->axis();
        const int ax1 = (ax0 + 1) % 3;
        const int ax2 = (ax0 + 2) % 3;
        const int nf  = _topo_hat[cdim]->nf();

        int istart[3];
        _topo_hat[cdim]->get_istart_glob(istart);

        // get the factors
        const double         normfact = _normfact;
        opt_double_ptr       mydata   = data;
        const opt_double_ptr mygreen  = _green;

        FLUPS_ASSUME_ALIGNED(mydata, FLUPS_ALIGNMENT);
        FLUPS_ASSUME_ALIGNED(mygreen, FLUPS_ALIGNMENT);
        const size_t onmax   = _topo_hat[cdim]->nloc(ax1) * _topo_hat[cdim]->nloc(ax2);
        const size_t inmax   = _topo_hat[cdim]->nloc(ax0);
        const int    nmem[3] = {_topo_hat[cdim]->nmem(0), _topo_hat[cdim]->nmem(1), _topo_hat[cdim]->nmem(2)};

        FLUPS_CHECK(FLUPS_ISALIGNED(mygreen) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
        FLUPS_CHECK(FLUPS_ISALIGNED(mydata) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);

        // get the memory size between two dimension
        const size_t memdim = (size_t)nmem[0] * (size_t)nmem[1] * (size_t)nmem[2];

        // do the loop
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, inmax, nmem, memdim, mydata, mygreen, normfact, ax0, ax1, ax2, nf, kfact, koffset)
        for (int io = 0; io < onmax; io++) {
            //local indexes start
            opt_double_ptr greenloc  = mygreen + collapsedIndex(ax0, 0, io, nmem, nf);
            opt_double_ptr dataloc_0 = mydata + 0 * memdim + collapsedIndex(ax0, 0, io, nmem, nf);
            opt_double_ptr dataloc_1 = mydata + 1 * memdim + collapsedIndex(ax0, 0, io, nmem, nf);
            opt_double_ptr dataloc_2 = mydata + 2 * memdim + collapsedIndex(ax0, 0, io, nmem, nf);
            // check the alignmeno
            FLUPS_ASSUME_ALIGNED(greenloc, FLUPS_ALIGNMENT);
            FLUPS_ASSUME_ALIGNED(dataloc_0, FLUPS_ALIGNMENT);
            FLUPS_ASSUME_ALIGNED(dataloc_1, FLUPS_ALIGNMENT);
            FLUPS_ASSUME_ALIGNED(dataloc_2, FLUPS_ALIGNMENT);

            // compute the k mode
            const double k1 = ((io % _topo_hat[cdim]->nloc(ax1)) + koffset[ax1]) * kfact[ax1];
            const double k2 = ((io / _topo_hat[cdim]->nloc(ax1)) + koffset[ax2]) * kfact[ax2];

            for (int ii = 0; ii < inmax; i0++) {
                const double k0 = (ii + koffset[ax0]) * kfact[ax0];
#if (KIND <= 0)
                // copy the values before the updtes
                const double gr = greenloc[ii];
                // real part field
                const double f0r = dataloc_0[ii];
                const double f1r = dataloc_1[ii];
                const double f2r = dataloc_2[ii];
                // compute the rotational
                const double divr = k0 * f0r + k1 * f1r + k2 * f2r;
#if (KIND == 0)
                dataloc_0[ii]     = normfact * divr * gr;
#else
                dataloc_0[ii]     = -normfact * divr * gr;
#endif

#else
                // copy the values before the updtes
                const double gr = greenloc[ii * 2 + 0];
                const double gc = greenloc[ii * 2 + 1];
                // real part field
                const double f0r = dataloc_0[ii * 2 + 0];
                const double f1r = dataloc_1[ii * 2 + 0];
                const double f2r = dataloc_2[ii * 2 + 0];
                // imag part field
                const double f0c = dataloc_0[ii * 2 + 1];
                const double f1c = dataloc_1[ii * 2 + 1];
                const double f2c = dataloc_2[ii * 2 + 1];
                // compute the rotational
                const double divr = k0 * f0r + k1 * f1r + k2 * f2r;
                const double divr = k0 * f0c + k1 * f1c + k2 * f2c;
                // update the values
#if (KIND == 1)
                dataloc_0[ii * 2 + 0] = normfact * (divr * gr - divc * gc);
                dataloc_0[ii * 2 + 1] = normfact * (divr * gc + divc * gr);
#elif (KIND == 2)
                dataloc_0[ii * 2 + 0] = -normfact * (divr * gr - divc * gc);
                dataloc_0[ii * 2 + 1] = -normfact * (divr * gc + divc * gr);
#elif (KIND == 3)
                dataloc_0[ii * 2 + 0] = -normfact * (divr * gc + divc * gr);
                dataloc_0[ii * 2 + 1] = normfact * (divr * gr - divc * gc);
#elif (KIND == 4)
                dataloc_0[ii * 2 + 0] = normfact * (divr * gc + divc * gr);
                dataloc_0[ii * 2 + 1] = -normfact * (divr * gr - divc * gc);
#endif
#endif
            }
        }
    }
    END_FUNC;
}
