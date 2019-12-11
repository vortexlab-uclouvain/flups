/**
 * @file dothemagic_rot.ipp
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

#if (KIND == 00)
void Solver::dothemagic_rot_real_p1(double *data, const double kfact[3], const double koffset[3], const double symstart[3], const int orderdiff, const bool* shift) {
#elif (KIND == 01)
void Solver::dothemagic_rot_real_m1(double *data, const double kfact[3], const double koffset[3], const double symstart[3], const int orderdiff, const bool* shift) {
#elif (KIND == 10)
void Solver::dothemagic_rot_complex_p1(double *data, const double kfact[3], const double koffset[3], const double symstart[3], const int orderdiff, const bool* shift) {
#elif (KIND == 11)
void Solver::dothemagic_rot_complex_m1(double *data, const double kfact[3], const double koffset[3], const double symstart[3], const int orderdiff, const bool* shift) {
#elif (KIND == 12)
void Solver::dothemagic_rot_complex_pi(double *data, const double kfact[3], const double koffset[3], const double symstart[3], const int orderdiff, const bool* shift) {
#elif (KIND == 13)
void Solver::dothemagic_rot_complex_mi(double *data, const double kfact[3], const double koffset[3], const double symstart[3], const int orderdiff, const bool* shift) {
#endif
    BEGIN_FUNC;
    
    // get the current topology
    Topology* topo = _topo_hat[_ndim - 1];

    // sanity checks
#if (KIND < 9)
    FLUPS_CHECK(topo->nf() == 1, "The topo_hat[2] (field) has to be complex", LOCATION);
#else
    FLUPS_CHECK(topo->nf() == 2, "The topo_hat[2] (field) has to be complex", LOCATION);
#endif

    // get the axis
    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nf  = topo->nf();

    //compute the memory offset due to 

    // get the starting indexes of the current block
    int istart[3];
    topo->get_istart_glob(istart);

    // get some factors and memory info
    const double normfact = _normfact;
    const size_t onmax    = topo->nloc(ax1) * topo->nloc(ax2);
    const size_t inmax    = topo->nloc(ax0);
    const int    nmem[3]  = {topo->nmem(0), topo->nmem(1), topo->nmem(2)};

    // get the pointers
    const opt_double_ptr mygreen = _green;
    const opt_double_ptr mydata0 = data;                       // this is the current dimension
    opt_double_ptr       mydata1 = data + 1 * topo->memdim();  // this will get d/dx2 data
    opt_double_ptr       mydata2 = data + 2 * topo->memdim();  // this will get - d/dx1 data

    FLUPS_CHECK(FLUPS_ISALIGNED(mygreen) && (nmem[ax0] * topo->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
    FLUPS_CHECK(FLUPS_ISALIGNED(mydata0) && (nmem[ax0] * topo->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
    FLUPS_CHECK(FLUPS_ISALIGNED(mydata1) && (nmem[ax0] * topo->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
    FLUPS_CHECK(FLUPS_ISALIGNED(mydata2) && (nmem[ax0] * topo->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);

    FLUPS_ASSUME_ALIGNED(mygreen, FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mydata0, FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mydata1, FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mydata2, FLUPS_ALIGNMENT);

    // get the memory size between two dimension
    const size_t nloc_ax1 = topo->nloc(ax1);

    // do the loop
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, inmax, nmem, memdim, mydata, myres1, myres2, mygreen, normfact, ax0, ax1, ax2, nf, kfact, koffset, istart, symstart, nloc_ax1)
    for (size_t io = 0; io < onmax; io++) {
        //local indexes start
        opt_double_ptr greenloc = mygreen + collapsedIndex(ax0, 0, io, nmem, nf);
        opt_double_ptr dataloc0  = mydata0 + collapsedIndex(ax0, 0, io, nmem, nf);
        opt_double_ptr dataloc1  = mydata1 + collapsedIndex(ax0, 0, io, nmem, nf);
        opt_double_ptr dataloc2  = mydata2 + collapsedIndex(ax0, 0, io, nmem, nf);
        // check the alignment
        FLUPS_ASSUME_ALIGNED(greenloc, FLUPS_ALIGNMENT);
        FLUPS_ASSUME_ALIGNED(dataloc0, FLUPS_ALIGNMENT);
        FLUPS_ASSUME_ALIGNED(dataloc1, FLUPS_ALIGNMENT);
        FLUPS_ASSUME_ALIGNED(dataloc2, FLUPS_ALIGNMENT);

        

        for (int ii = 0; ii < inmax; ii++) {
            int is[3];
            cmpt_symID(ax0, ii, io % nloc_ax1, io / nloc_ax1, istart, symstart, 0, is);

#if (KIND < 9)
            //----------------------------------------------------
            // REAL
            //----------------------------------------------------
            // copy the real values before the updtes
            const double gr = greenloc[ii];
            // real part field
            const double f0r = dataloc_0[ii];
            const double f1r = dataloc_1[ii];
            const double f2r = dataloc_2[ii];
            // compute the rotational
            // ...this case only happens when we do 3dirspectral:
            // combining rephasing (which will always be i or -i)
            // and the derivative (which involves i), what we
            // actually do is computing (k) x (f)
            // and then choosing the correct sign
            const double rot0r = k1 * f2r - k2 * f1r;
            const double rot1r = k2 * f0r - k0 * f2r;
            const double rot2r = k0 * f1r - k1 * f0r;
#if (KIND == 00)
            dataloc_0[ii] = normfact * rot0r * gr;
            dataloc_1[ii] = normfact * rot1r * gr;
            dataloc_2[ii] = normfact * rot2r * gr;
#elif (KIND == 01)
            dataloc_0[ii] = -normfact * rot0r * gr;
            dataloc_1[ii] = -normfact * rot1r * gr;
            dataloc_2[ii] = -normfact * rot2r * gr;
#endif

#else
            //----------------------------------------------------
            // COMPLEX
            //----------------------------------------------------
            // copy the values before the updates
            const double gr  = greenloc[ii * 2 + 0];
            const double gc  = greenloc[ii * 2 + 1];
            const double f0r = dataloc_0[ii * 2 + 0];
            const double f1r = dataloc_1[ii * 2 + 0];
            const double f2r = dataloc_2[ii * 2 + 0];
            const double f0c = dataloc_0[ii * 2 + 1];
            const double f1c = dataloc_1[ii * 2 + 1];
            const double f2c = dataloc_2[ii * 2 + 1];

            // (symmetrized) wave number
            // ki = k in the ith direction and NOT the k aligned in the axis!!
            // the offset is different for each component of the field
            // ki_cj is the k in the direction i for the component j
            const double k0_c1 = ((double)is[0] + koffset[0][1]) * kfact[0];
            const double k0_c2 = ((double)is[0] + koffset[0][2]) * kfact[0];
            const double k1_c0 = ((double)is[1] + koffset[1][0]) * kfact[1];
            const double k1_c2 = ((double)is[1] + koffset[1][2]) * kfact[1];
            const double k2_c0 = ((double)is[2] + koffset[2][0]) * kfact[2];
            const double k2_c1 = ((double)is[2] + koffset[2][1]) * kfact[2];



            // // compute the rotational:  (ik) x (f)
            const double rot0r = -k1 * f2c + k2 * f1c;
            const double rot0c = +k1 * f2r - k2 * f1r;
            const double rot1r = -k2 * f0c + k0 * f2c;
            const double rot1c = +k2 * f0r - k0 * f2r;
            const double rot2r = -k0 * f1c + k1 * f0c;
            const double rot2c = +k0 * f1r - k1 * f0r;

            // update the values
            // accounting for the rephasing:
#if (KIND == 10)
            dataloc_0[ii * 2 + 0] = normfact * (rot0r * gr - rot0c * gc);
            dataloc_0[ii * 2 + 1] = normfact * (rot0r * gc + rot0c * gr);
            dataloc_1[ii * 2 + 0] = normfact * (rot1r * gr - rot1c * gc);
            dataloc_1[ii * 2 + 1] = normfact * (rot1r * gc + rot1c * gr);
            dataloc_2[ii * 2 + 0] = normfact * (rot2r * gr - rot2c * gc);
            dataloc_2[ii * 2 + 1] = normfact * (rot2r * gc + rot2c * gr);
#elif (KIND == 11)
            dataloc_0[ii * 2 + 0] = -normfact * (rot0r * gr - rot0c * gc);
            dataloc_0[ii * 2 + 1] = -normfact * (rot0r * gc + rot0c * gr);
            dataloc_1[ii * 2 + 0] = -normfact * (rot1r * gr - rot1c * gc);
            dataloc_1[ii * 2 + 1] = -normfact * (rot1r * gc + rot1c * gr);
            dataloc_2[ii * 2 + 0] = -normfact * (rot2r * gr - rot2c * gc);
            dataloc_2[ii * 2 + 1] = -normfact * (rot2r * gc + rot2c * gr);
#elif (KIND == 12)
            dataloc_0[ii * 2 + 0] = -normfact * (rot0r * gc + rot0c * gr);
            dataloc_0[ii * 2 + 1] = +normfact * (rot0r * gr - rot0c * gc);
            dataloc_1[ii * 2 + 0] = -normfact * (rot1r * gc + rot1c * gr);
            dataloc_1[ii * 2 + 1] = +normfact * (rot1r * gr - rot1c * gc);
            dataloc_2[ii * 2 + 0] = -normfact * (rot2r * gc + rot2c * gr);
            dataloc_2[ii * 2 + 1] = +normfact * (rot2r * gr - rot2c * gc);

#elif (KIND == 13)
            dataloc_0[ii * 2 + 0] = +normfact * (rot0r * gc + rot0c * gr);
            dataloc_0[ii * 2 + 1] = -normfact * (rot0r * gr - rot0c * gc);
            dataloc_1[ii * 2 + 0] = +normfact * (rot1r * gc + rot1c * gr);
            dataloc_1[ii * 2 + 1] = -normfact * (rot1r * gr - rot1c * gc);
            dataloc_2[ii * 2 + 0] = +normfact * (rot2r * gc + rot2c * gr);
            dataloc_2[ii * 2 + 1] = -normfact * (rot2r * gr - rot2c * gc);
#endif
#endif
        }
    }
    END_FUNC;
}
